Vararg_Tuple{T} = Tuple{Vararg{T}}
compute_idx(index::Vararg_Tuple{T}, t::T, ::OneConfig) where {T} = index
compute_idx(index::Vararg_Tuple{T}, t::T, ::ReconfigAllowed) where {T} = (t, index...)
# ---------------------------------------------------------------------------- #
#                         Voltage Reference Constraints                        #
# ---------------------------------------------------------------------------- #
# 1 pu reference at each substation if its is built 
# don't forget to think also to add the fact that some substations can already
# be there
function _add_RefVoltages!(model::JuMP.AbstractModel)::Nothing
    network_data = model[:network_data]
    T  = model[:time_steps]
    Ns = network_data.nb_substations
   
    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], model[:V_sqr][t, i] - 1 <= 
                                        (network_data.buses[i].V_limits.V_max^2 - 1) *
                                        (1 - model[:Beta][i])
                        [t=1:T, i=1:Ns], model[:V_sqr][t, i] - 1 >= 
                                        (network_data.buses[i].V_limits.V_min^2 - 1) * 
                                        (1 - model[:Beta][i])
                    end)
    return 
end

# ---------------------------------------------------------------------------- #
#                        Load Oversatisfaction constraints                     #
# ---------------------------------------------------------------------------- #
# ----------------------------------- NoDG ----------------------------------- #
function _add_LoadOverSatisfaction!(model::JuMP.AbstractModel, ::NoDG)::Nothing
    network_data = model[:network_data]
    T  = model[:time_steps]
    Ns = network_data.nb_substations
    N  = Ns + network_data.nb_loads
    buses = network_data.buses

    JuMP.@constraints(  model, begin
                        [t=1:T], sum(model[:P_sub][t, i] for i in 1:Ns) >= 
                                 sum(buses[i].load_profile.time_serie[t] .*
                                     buses[i].cos_phi for i in (Ns+1):N)
                        
                        [t=1:T], sum(model[:Q_sub][t, i] for i in 1:Ns) >= 
                                 sum(buses[i].load_profile.time_serie[t] .* 
                                 sin(acos(buses[i].cos_phi)) for i in (Ns+1):N)
                    end)
    return 
end
# ------------------------------------ DG ------------------------------------ #
function _add_LoadOverSatisfaction!(model::JuMP.AbstractModel, ::DG)::Nothing
    network_data = model[:network_data]
    T  = model[:time_steps]
    Ns = network_data.nb_substations
    Nu = network_data.nb_loads
    N  = Ns + Nu
    buses = network_data.buses

    JuMP.@constraints(  model, begin
                        [t=1:T], sum(model[:P_sub][t, i] for i in 1:Ns) + 
                                 sum(model[:p_pv][t, i] for i in 1:Nu) >= 
                                 sum(buses[i].load_profile.time_serie[t] .*
                                 buses[i].cos_phi for i in (Ns+1):N)
                        
                        [t=1:T], sum(model[:Q_sub][t, i] for i in 1:Ns) >= 
                                 sum(buses[i].load_profile.time_serie[t] .* 
                                 sin(acos(buses[i].cos_phi)) for i in (Ns+1):N)
                    end)
    return 
end

# ---------------------------------------------------------------------------- #
#                               Substation constraints                         #
# ---------------------------------------------------------------------------- #
# Think about initial substations !!!!
function _add_SubstationConstraints!(model::JuMP.AbstractModel, ::NonConvex)::Nothing

    network_data = model[:network_data]
    T  = model[:time_steps]
    Ns = network_data.nb_substations
    buses = network_data.buses

    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], model[:S_sub][t, i]^2 == model[:P_sub][t, i]^2 + model[:Q_sub][t, i]^2
                        [t=1:T, i=1:Ns], model[:S_sub][t, i] <= model[:S_sub_capa][i]
                        [i=1:Ns], model[:S_sub_capa][i] <= model[:Beta][i] * buses[i].S_rating_max
                    end)
    return 
end

function _add_SubstationConstraints!(model::JuMP.AbstractModel, ::Convex)::Nothing

    network_data = model[:network_data]
    T  = model[:time_steps]
    Ns = network_data.nb_substations
    buses = network_data.buses

    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], [model[:S_sub][t, i]; model[:P_sub][t, i]; model[:Q_sub][t, i]] in JuMP.SecondOrderCone()
                        [t=1:T, i=1:Ns], model[:S_sub][t, i] <= model[:S_sub_capa][i]
                        [i=1:Ns], model[:S_sub_capa][i] <= model[:Beta][i] * buses[i].S_rating_max
                    end)
    return 
end
# ---------------------------------------------------------------------------- #
#                     Voltage Operational constraints                          #
# ---------------------------------------------------------------------------- #
# Voltage operational constraints 
# Add the possibility to be relaxed or not
function _add_VoltageOpConstraints!(model::JuMP.AbstractModel, ::StrongVoltages)::Nothing

    network_data = model[:network_data]
    T  = model[:time_steps]
    N  = network_data.nb_substations + network_data.nb_loads
    buses = network_data.buses

    for t in 1:T, i in 1:N
        JuMP.set_lower_bound(model[:V_sqr][t, i], buses[i].V_limits.V_min^2)
        JuMP.set_upper_bound(model[:V_sqr][t, i], buses[i].V_limits.V_max^2)
    end
    return 
end

################## TO DO ##############################
#######################################################
function _add_VoltageOpConstraints!(model::JuMP.AbstractModel, ::RelaxedVoltages)::Nothing

    network_data = model[:network_data]
    T  = model[:time_steps]
    N  = network_data.nb_substations + network_data.nb_loads
    buses = network_data.buses

    JuMP.@variables( model, begin 
                    V_violation_up[1:T, 1:L, 1:K], (binary=true)
                    V_violation_low[1:T, 1:L, 1:K], (binary=true)
                    V_slack_up[1:T, 1:L, 1:K]
                    V_slack_low[1:T, 1:L, 1:K]
                    end)


    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:N], model[:V_slack_low][t, i] <= 
                                        model[:V_violation_low][t, i] * buses[i].V_limits.V_min^2
                        [t=1:T, i=1:N], model[:V_slack_up][t, i] <= 
                                        model[:V_violation_up][t, i] * buses[i].V_limits.V_max^2
                        [t=1:T, i=1:N], model[:V_sqr][t, i] + model[:V_slack_low][t, i] >= 
                                        buses[i].V_limits.V_min^2
                        [t=1:T, i=1:N], model[:V_sqr][t, i] - model[:V_slack_up][t, i] <= 
                                        buses[i].V_limits.V_max^2
                    end)
    return 
end


# ---------------------------------------------------------------------------- #
#                     Current Operational constraints                          #
# ---------------------------------------------------------------------------- #
function _add_CurrentOpConstraints!(model::JuMP.AbstractModel, 
                                    topology_choice::TopologyChoiceFormulation, 
                                    ::StrongCurrents
                                )::Nothing

    network_data = model[:network_data]
    T  = model[:time_steps]
    L  = network_data.nb_lines
    K  = network_data.nb_conductors
    conductors = network_data.conductors

    cond_choice = isa(topology_choice, OneConfig) ? model[:Alpha] : model[:Gamma]
    
    JuMP.@constraint(model, 
                    [t=1:T, l=1:L, k=1:K],
                    model[:I_sqr_k][t, l, k] <= conductors[k].max_i^2 * 
                    cond_choice[compute_idx((l, k), t, topology_choice)...])

    return 
end

function _add_CurrentOpConstraints!(model::JuMP.AbstractModel, 
                                    topology_choice::TopologyChoiceFormulation, 
                                    ::RelaxedCurrents
                                )::Nothing

    network_data = model[:network_data]
    T = model[:time_steps]
    L = network_data.nb_lines
    K = network_data.nb_conductors
    conductors = network_data.conductors

    cond_choice = isa(topology_choice, OneConfig) ? model[:Alpha] : model[:Gamma]

    JuMP.@variables( model, begin 
                    I_violation[1:T, 1:L, 1:K], (binary=true)
                    I_slack[1:T, 1:L, 1:K]
                    end)

    JuMP.@constraint(model, 
                    [t=1:T, l=1:L, k=1:K],  
                    model[:I_slack][t, l, k] <= conductors[k].max_i^2 * model[:I_violation][t, l, k])


    JuMP.@constraint(model,
                    [t=1:T, l=1:L, k=1:K],  
                    model[:I_sqr_k][t, l, k] - model[:I_slack][t, l, k] <= conductors[k].max_i^2 * 
                    cond_choice[compute_idx((l, k), t, topology_choice)...])


    return 
end

# ---------------------------------------------------------------------------- #
#                               PowerBalance constraints                       #
# ---------------------------------------------------------------------------- #
function _add_PowerBalanceConstraints!( model::JuMP.AbstractModel, 
                                        prod_type::TypeofProdFormulation, 
                                        ::BFM
                                    )::Nothing
   
    Omega_sending = model[:network_topology].sending_lines
    Omega_receiving = model[:network_topology].receiving_lines
    T  = model[:time_steps]
    network_data = model[:network_data]
    K = network_data.nb_conductors
    Ns = network_data.nb_substations
    Nu = network_data.nb_loads
    N = Ns + Nu
    L = network_data.nb_lines
    conductors = network_data.conductors
    buses = network_data.buses
    lines = network_data.lines
    P_consumed = [buses[i].load_profile.time_serie .* buses[i].cos_phi for i in (Ns+1):N]
    Q_consumed = [buses[i].load_profile.time_serie .* sin(acos(buses[i].cos_phi)) for i in (Ns+1):N]
    p_pv = isa(prod_type, DG) ? model[:p_pv] : zeros(Float64, T, Nu)

    #println(Q_consumed)

    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], -  model[:P_sub][t, i] ==  
                                            sum((model[:P_ij_k][t, l, k] - 
                                            conductors[k].r * lines[l].length * model[:I_sqr_k][t, l, k])
                                            for l in Omega_receiving[i], k in 1:K) -
                                            sum(model[:P_ij_k][t, l, k] for l in Omega_sending[i], k in 1:K)

                        [t=1:T, i=1:Ns], -  model[:Q_sub][t, i] ==  
                                            sum((model[:Q_ij_k][t, l, k] - 
                                            conductors[k].x * lines[l].length * model[:I_sqr_k][t, l, k])
                                            for l in Omega_receiving[i], k in 1:K) -
                                            sum(model[:Q_ij_k][t, l, k] for l in Omega_sending[i], k in 1:K)

                        [t=1:T, i=1:Nu], -  p_pv[t, i] + P_consumed[i][t] ==  
                                            sum(model[:P_ij_k][t, l, k] - 
                                            conductors[k].r * lines[l].length * model[:I_sqr_k][t, l, k]
                                            for l in Omega_receiving[Ns + i], k in 1:K) -
                                            sum(model[:P_ij_k][t, l, k] for l in Omega_sending[Ns + i], k in 1:K)

                        [t=1:T, i=1:Nu],    Q_consumed[i][t] ==  
                                            sum(model[:Q_ij_k][t, l, k] - 
                                            conductors[k].x * lines[l].length * model[:I_sqr_k][t, l, k]
                                            for l in Omega_receiving[Ns + i], k in 1:K) -
                                            sum(model[:Q_ij_k][t, l, k] for l in Omega_sending[Ns + i], k in 1:K)
                        
                        # Maybe move this at another place in the code but temporary to test
                        [t=1:T, l=1:L, k=1:K], model[:P_ij_k][t, l, k] <= conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max * model[:Alpha][l, k] # indispensable
                        [t=1:T, l=1:L, k=1:K], model[:P_ij_k][t, l, k] >= -conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max * model[:Alpha][l, k] # indispensable
                        [t=1:T, l=1:L, k=1:K], model[:Q_ij_k][t, l, k] <= conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max * model[:Alpha][l, k] # indispensable
                        [t=1:T, l=1:L, k=1:K], model[:Q_ij_k][t, l, k] >= -conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max * model[:Alpha][l, k] # indispensable
                        end
                    )
    return
end

function _add_PowerBalanceConstraints!( model::JuMP.AbstractModel, 
                                        prod_type::TypeofProdFormulation, 
                                        ::BIM
                                    )::Nothing


    Omega_sending = model[:network_topology].sending_lines
    Omega_receiving = model[:network_topology].receiving_lines
    network_data = model[:network_data]
    T  = model[:time_steps]
    K = network_data.nb_conductors
    Ns = network_data.nb_substations
    Nu = network_data.nb_loads
    buses = network_data.buses
    P_consumed = [buses[Ns + i].load_profile.time_serie .* buses[Ns + i].cos_phi for i in 1:Nu]
    Q_consumed = [buses[Ns + i].load_profile.time_serie .* sin(acos(buses[Ns + i].cos_phi)) for i in 1:Nu]
    p_pv = isa(prod_type, DG) ? model[:p_pv] : zeros(Float64, T, Nu)


    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], -  model[:P_sub][t, i] ==  
                                            sum(model[:P_ji_k][t, l, k] for l in Omega_receiving[i], k in 1:K) -
                                            sum(model[:P_ij_k][t, l, k] for l in Omega_sending[i], k in 1:K)

                        [t=1:T, i=1:Ns], -  model[:Q_sub][t, i] ==  
                                            sum(model[:Q_ji_k][t, l, k] for l in Omega_receiving[i], k in 1:K) -
                                            sum(model[:Q_ij_k][t, l, k] for l in Omega_sending[i], k in 1:K)

                        [t=1:T, i=1:Nu], -  p_pv[t, i] + P_consumed[i][t] ==  
                                            sum(model[:P_ji_k][t, l, k] for l in Omega_receiving[Ns + i], k in 1:K) -
                                            sum(model[:P_ij_k][t, l, k] for l in Omega_sending[Ns + i], k in 1:K)

                        [t=1:T, i=1:Nu],    Q_consumed[i][t] ==  
                                            sum(model[:Q_ji_k][t, l, k] for l in Omega_receiving[Ns + i], k in 1:K) -
                                            sum(model[:Q_ij_k][t, l, k] for l in Omega_sending[Ns + i], k in 1:K)
                    end)

                   
    return 
end

# ---------------------------------------------------------------------------- #
#                               PowerFlow constraints                          #
# ---------------------------------------------------------------------------- #
function _add_RotatedConicConstraints!( model::JuMP.AbstractModel, 
                                        ::BFM, ::Convex
                                    )::Nothing

    network_data = model[:network_data]
    T  = model[:time_steps]
    L = network_data.nb_lines
    K = network_data.nb_conductors
    lines = network_data.lines


    JuMP.@constraint(model, 
                    [t=1:T, l=1:L],
                    [model[:V_sqr][t, lines[l].edge.from_node.id] / 2; 
                    sum(model[:I_sqr_k][t, l, k] for k in 1:K); 
                    sum(model[:P_ij_k][t, l, k] for k in 1:K); 
                    sum(model[:Q_ij_k][t, l, k] for k in 1:K)] in 
                    JuMP.RotatedSecondOrderCone())
    return 
end

function _add_RotatedConicConstraints!(model::JuMP.AbstractModel, 
                                        ::BFM, ::NonConvex
                                    )::Nothing
    network_data = model[:network_data]
    T  = model[:time_steps]
    L = network_data.nb_lines
    K = network_data.nb_conductors
    lines = network_data.lines

    JuMP.@constraint(model,
                    [t=1:T, l=1:L],
                    model[:V_sqr][t, lines[l].edge.from_node.id] * 
                    sum(model[:I_sqr_k][t, l, k] for k in 1:K) == 
                    sum(model[:P_ij_k][t, l, k] for k in 1:K)^2 + 
                    sum(model[:Q_ij_k][t, l, k] for k in 1:K)^2)
    return 
end

function _add_RotatedConicConstraints!( model::JuMP.AbstractModel, 
                                        ::BIM, ::Convex
                                    )::Nothing
    network_data = model[:network_data]
    T  = model[:time_steps]
    L = network_data.nb_lines
    K = network_data.nb_conductors
    lines = network_data.lines
   
    JuMP.@constraint(model, 
                    [t=1:T, l=1:L],
                    [sum(model[:X_ij_k_i][t, l, k, lines[l].edge.from_node.id] for k in 1:K) / 2; 
                    sum(model[:X_ij_k_i][t, l, k, lines[l].edge.to_node.id] for k in 1:K); 
                    sum(model[:X_ij_k_re][t, l, k] for k in 1:K); 
                    sum(model[:X_ij_k_im][t, l, k] for k in 1:K)] in 
                    JuMP.RotatedSecondOrderCone())
    return 
end

function _add_RotatedConicConstraints!( model::JuMP.AbstractModel, 
                                        ::BIM, ::NonConvex
                                    )::Nothing
    network_data = model[:network_data]
    T  = model[:time_steps]
    L = network_data.nb_lines
    K = network_data.nb_conductors
    lines = network_data.lines

    JuMP.@constraint(model, 
                    [t=1:T, l=1:L],
                    sum(model[:X_ij_k_i][t, l, k, lines[l].edge.from_node.id] for k in 1:K) * 
                    sum(model[:X_ij_k_i][t, l, k, lines[l].edge.to_node.id] for k in 1:K) == 
                    sum(model[:X_ij_k_re][t, l, k] for k in 1:K)^2 + 
                    sum(model[:X_ij_k_im][t, l, k] for k in 1:K)^2)
    return 
end

function _add_PowerFlowConstraints!(model::JuMP.AbstractModel, 
                                    topology_choice::TopologyChoiceFormulation, 
                                    graph::TypeOfGraph,
                                    ::BIM)::Nothing

    # THINK ABOUT WAYS TO GET THE INDICES OF ALPHA 
    network_data = model[:network_data]
    T  = model[:time_steps]
    L = network_data.nb_lines
    K = network_data.nb_conductors
    lines = network_data.lines
    conductors = network_data.conductors
    buses = network_data.buses
    cond_choice = isa(topology_choice, OneConfig) ? model[:Alpha] : model[:Gamma]

    for l in 1:L 
        ifrom = lines[l].edge.from_node.id
        ito   = lines[l].edge.to_node.id
        for k in 1:K
            admittance = 1/(conductors[k].r * lines[l].length + im * (conductors[k].x * lines[l].length))
            G = real(admittance) 
            B = abs(imag(admittance))
            for t in 1:T
                index = compute_idx((l, k), t, topology_choice)

                JuMP.@constraints(model, begin
                                
                    model[:P_ij_k][t, l, k] ==  G * (model[:X_ij_k_i][t, l, k, ifrom] - 
                                                model[:X_ij_k_re][t, l, k]) + 
                                                B * model[:X_ij_k_im][t, l, k]

                    model[:P_ji_k][t, l, k] ==  G * (model[:X_ij_k_i][t, l, k, ito] - 
                                                model[:X_ij_k_re][t, l, k]) - 
                                                B * model[:X_ij_k_im][t, l, k]

                    model[:Q_ij_k][t, l, k] ==  B * (model[:X_ij_k_i][t, l, k, ifrom] - 
                                                model[:X_ij_k_re][t, l, k]) - 
                                                G * model[:X_ij_k_im][t, l, k]

                    model[:Q_ji_k][t, l, k] ==  B * (model[:X_ij_k_i][t, l, k, ifrom] - 
                                                model[:X_ij_k_re][t, l, k]) + 
                                                G * model[:X_ij_k_im][t, l, k]

                    model[:I_sqr_k][t, l, k] == (G^2 + B^2) * (model[:X_ij_k_i][t, l, k, ifrom] + 
                                                model[:X_ij_k_i][t, l, k, ito] - 
                                                2 * model[:X_ij_k_re][t, l, k])

                    model[:X_ij_k_i][t, l, k, ifrom] >= buses[ifrom].V_limits.V_min^2 * 
                                                        cond_choice[index...]

                    model[:X_ij_k_i][t, l, k, ifrom] <= buses[ifrom].V_limits.V_max^2 * 
                                                        cond_choice[index...]

                    model[:X_ij_k_i][t, l, k, ito] >=   buses[ito].V_limits.V_min^2 * 
                                                        cond_choice[index...]

                    model[:X_ij_k_i][t, l, k, ito] <=   buses[ito].V_limits.V_max^2 * 
                                                        cond_choice[index...]

                    model[:X_ij_k_re][t, l, k] <=   buses[ifrom].V_limits.V_max * 
                                                    buses[ito].V_limits.V_max * 
                                                    cond_choice[index...]

                    model[:X_ij_k_im][t, l, k] <=   buses[ifrom].V_limits.V_max * 
                                                    buses[ito].V_limits.V_max * 
                                                    cond_choice[index...]

                    model[:X_ij_k_im][t, l, k] >=   - buses[ifrom].V_limits.V_max * 
                                                    buses[ito].V_limits.V_max * 
                                                    cond_choice[index...]

                    model[:V_sqr][t, ifrom] - model[:X_ij_k_i][t, l, k, ifrom] >= 
                        buses[ifrom].V_limits.V_min^2 * (1 - cond_choice[index...])
    
                    model[:V_sqr][t, ifrom] - model[:X_ij_k_i][t, l, k, ifrom] <= 
                        buses[ifrom].V_limits.V_max^2 * (1 - cond_choice[index...])
        
                    model[:V_sqr][t, ito] - model[:X_ij_k_i][t, l, k, ito] >= 
                        buses[ito].V_limits.V_min^2 * (1 - cond_choice[index...])
        
                    model[:V_sqr][t, ito] - model[:X_ij_k_i][t, l, k, ito] <= 
                            buses[ito].V_limits.V_max^2 * (1 - cond_choice[index...])

                    end)
            end
        end
    end
    
    return 
end


function _add_PowerFlowConstraints!(model::JuMP.AbstractModel, 
                                    topology_choice::TopologyChoiceFormulation,
                                    graph::TypeOfGraph,
                                    ::BFM
                                )::Nothing

    # Redo this function and check several times !!!!!

    network_data = model[:network_data]
    T  = model[:time_steps]
    L = network_data.nb_lines
    K = network_data.nb_conductors
    conductors = network_data.conductors
    lines = network_data.lines
    buses = network_data.buses

    # [l = L, t = T], V_sqr[line_ends[l][2], t] - V_sqr[line_ends[l][1], t] <=
    #                 sum(-2 * (R[l, k] * P_ij_k[l, k, t] + X[l, k] * Q_ij_k[l, k, t]) +
    #                      (R[l, k]^2 + X[l, k]^2) * I_sqr_k[l, k, t] for k in K) +
    #                 M * (1 - Y[l])
    # [l = L, t = T], V_sqr[line_ends[l][2], t] - V_sqr[line_ends[l][1], t] >=
    #                 sum(-2 * (R[l, k] * P_ij_k[l, k, t] + X[l, k] * Q_ij_k[l, k, t]) +
    #                      (R[l, k]^2 + X[l, k]^2) * I_sqr_k[l, k, t] for k in K) -
    #                 M * (1 - Y[l])


    voltage_expr = JuMP.@expression(
        model, 
        [t=1:T, l=1:L], sum( -2 * (conductors[k].r * lines[l].length * model[:P_ij_k][t, l, k] + 
        conductors[k].x * lines[l].length * model[:Q_ij_k][t, l, k]) +
        ((conductors[k].r * lines[l].length)^2 + (conductors[k].x * lines[l].length)^2) * 
        model[:I_sqr_k][t, l, k] for k in 1:K)) 

    if isa(graph, Undirected)
        Y = JuMP.@expression(  model, 
                                [t=1:T, l=1:L], 
                                (1 - model[:Y][compute_idx((l,), t, topology_choice)...]) * 
                                (buses[lines[l].edge.to_node.id].V_limits.V_max^2 - 
                                buses[lines[l].edge.from_node.id].V_limits.V_min^2)
                            )

    elseif isa(graph, Directed)
        Y = JuMP.@expression(  model, 
                                [t=1:T, l=1:L], 
                                (1 - (model[:Y_send][compute_idx((l,), t, topology_choice)...] + 
                                model[:Y_rec][compute_idx((l,), t, topology_choice)...])) *
                                (buses[lines[l].edge.to_node.id].V_limits.V_max^2 - 
                                buses[lines[l].edge.from_node.id].V_limits.V_min^2)
                            )
    end

    for l in 1:L 
        ifrom = lines[l].edge.from_node.id
        ito   = lines[l].edge.to_node.id
        for t in 1:T
            JuMP.@constraints(  model, begin
                model[:V_sqr][t, ito] - model[:V_sqr][t, ifrom] <= voltage_expr[t, l] + 
                                            Y[compute_idx((l,), t, topology_choice)...]

                model[:V_sqr][t, ito] - model[:V_sqr][t, ifrom] >= voltage_expr[t, l] - 
                                            Y[compute_idx((l,), t, topology_choice)...]
                
            end)
        end
    end
    return 
end

# ---------------------------------------------------------------------------- #
#                               PV Operation constraints                       #
# ---------------------------------------------------------------------------- #
# I AM HERE : ERROR WITH PV_PROD
function _add_PVOperationConstraints!(model::JuMP.AbstractModel)::Nothing
    network_data = model[:network_data]
    T  = model[:time_steps]
    Nu = network_data.nb_loads
    Ns = network_data.nb_substations
    PV_prod = [network_data.buses[Ns + i].PV_installation.profile.time_serie for i in 1:Nu]

    JuMP.@constraints(model, begin
                    [t=1:T, i=1:Nu], model[:p_pv][t, i] <= model[:s_conv_pv][i] # The PV power output at time t is always bounded by the capacity of the converter
                    [t=1:T, i=1:Nu], model[:p_pv][t, i] <= PV_prod[i][t] * model[:p_pv_max][i] # The PV power output at time t is also bounded by the available sun power PV_PRODUCTION is in [W/wp]
    end)
    
    return
end

# ---------------------------------------------------------------------------- #
#                               Radiality constraints                          #
# ---------------------------------------------------------------------------- #
#=
function _add_RadialityConstraints!(model::JuMP.AbstractModel, 
                                    topology_choice::TopologyChoiceFormulation,
                                    graph_type::TypeOfGraph,
                                    ::RadialityFormulation)::Nothing

    
    return
end
=#