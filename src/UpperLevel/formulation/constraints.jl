# ---------------------------------------------------------------------------- #
#                         Voltage Reference Constraints                        #
# ---------------------------------------------------------------------------- #
# 1 pu reference at each substation if its is built 
# don't forget to think also to add the fact that some substations can already
# be there
function _add_RefVoltages!(model::JuMP.AbstractModel)::Nothing
    network_data = model[:network_data]
    T  = model[:time_steps]
    Ns = get_nb_substations(network_data)
   
    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], model[:V_sqr][t, i] - 1 <= 
                                        (network.buses[i].V_limits.V_max^2 - 1) *
                                        (1 - model[:Beta][i])
                        [t=1:T, i=1:Ns], model[:V_sqr][t, i] - 1 >= 
                                        (network.buses[i].V_limits.V_min^2 - 1) * 
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
    Ns = get_nb_substations(network_data)
    N  = get_nb_buses(network_data)
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
    Nu = get_nb_load_bus(network_data)
    Ns = get_nb_sub_bus(network_data)
    buses = network_data.buses

    JuMP.@constraints(  model, begin
                        [t=1:T], sum(model[:P_sub][t, i] for i in 1:Ns) + 
                                 sum(model[:p_pv][t, i] for i in (Ns+1):N) >= 
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
    Ns = get_nb_substations(network_data)
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
    Ns = get_nb_substations(network_data)
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
    N  = get_nb_buses(network_data)
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
    N  = get_nb_buses(network_data)
    buses = network_data.buses

    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], model[:V_sqr][t, i] >= sub_buses[i].V_limits.V_min^2
                        [t=1:T, i=1:Ns], model[:V_sqr][t, i] <= sub_buses[i].V_limits.V_max^2
                        [t=1:T, i=1:Nu], model[:V_sqr][t, Ns + i] >= load_buses[i].V_limits.V_min^2
                        [t=1:T, i=1:Nu], model[:V_sqr][t, Ns + i] <= load_buses[i].V_limits.V_max^2
                    end)
    return 
end


# ---------------------------------------------------------------------------- #
#                     Current Operational constraints                          #
# ---------------------------------------------------------------------------- #
function _add_CurrentOpConstraints!(model::JuMP.AbstractModel, 
                                    choice_topology::TopologyChoiceFormulation, 
                                    ::StrongCurrents
                                )::Nothing

    network_data = model[:network_data]
    T  = model[:time_steps]
    L  = get_nb_lines(network_data)
    K  = get_nb_conductors(network_data)
    conductors = network_data.conductors

    JuMP.@constraint(model, 
                    [t=1:T, l=1:L, k=1:K],
                    model[:I_sqr_k][t, l, k] <= conductors[k].max_i^2 * 
                    model[:Alpha][compute_index((l, k), t, choice_topology)...])
    return 
end

function _add_CurrentOpConstraints!(model::JuMP.AbstractModel, 
                                    choice_topology::TopologyChoiceFormulation, 
                                    ::RelaxedCurrents
                                )::Nothing

    network_data = model[:network_data]
    T = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)
    conductors = network_data.conductors

    JuMP.@variables( model, begin 
                    I_violation[1:T, 1:L, 1:K], (binary=true)
                    I_slack[1:T, 1:L, 1:K]
                    end)

    JuMP.@constraints(  model, begin
                        [t=1:T, l=1:L, k=1:K],  model[:I_slack][t, l, k] <= 
                                                conductors[k].max_i^2 * model[:I_violation][t, l, k]
                        [t=1:T, l=1:L, k=1:K],  model[:I_sqr_k][t, l, k] - model[:I_slack][t, l, k] <= 
                                                conductors[k].max_i^2 * 
                                                model[:Alpha][compute_index((l, k), t, choice_topology)...] # indispensable mais relaxée
                    end)
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
    K = get_nb_conductors(network_data)
    Ns = get_nb_substations(network_data)
    Nu = get_nb_loads(network_data)
    conductors = network_data.conductors
    lines = network_data.lines

    p_pv = isa(prod_type, DG) ? model[:p_pv] : zeros(Float64, T, Nu)

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

                        [t=1:T, i=1:Nu],    p_pv[t, i] ==  
                                            sum(model[:P_ij_k][t, l, k] - 
                                            conductors[k].r * lines[l].length * model[:I_sqr_k][t, l, k]
                                            for l in Omega_receiving[Ns + i], k in 1:K) -
                                            sum(model[:P_ij_k][t, l, k] for l in Omega_sending[Ns + i], k in 1:K)

                        [t=1:T, i=1:Nu],    0 ==  
                                            sum(model[:Q_ij_k][t, l, k] - 
                                            conductors[k].x * lines[l].length * model[:I_sqr_k][t, l, k]
                                            for l in Omega_receiving[Ns + i], k in 1:K) -
                                            sum(model[:Q_ij_k][t, l, k] for l in Omega_sending[Ns + i], k in 1:K)
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
    K = get_nb_conductors(network_data)
    Ns = get_nb_sub_bus(network_data)
    Nu = get_nb_load_bus(network_data)
 
    p_pv = isa(prod_type, DG) ? model[:p_pv] : zeros(Float64, T, Nu)

    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], -  model[:P_sub][t, i] ==  
                                            sum(model[:P_ji_k][t, l, k] for l in Omega_receiving[i], k in 1:K) -
                                            sum(model[:P_ij_k][t, l, k] for l in Omega_sending[i], k in 1:K)

                        [t=1:T, i=1:Ns], -  model[:Q_sub][t, i] ==  
                                            sum(model[:Q_ji_k][t, l, k] for l in Omega_receiving[i], k in 1:K) -
                                            sum(model[:Q_ij_k][t, l, k] for l in Omega_sending[i], k in 1:K)

                        [t=1:T, i=1:Nu],    p_pv[t, i] ==  
                                            sum(model[:P_ji_k][t, l, k] for l in Omega_receiving[Ns + i], k in 1:K) -
                                            sum(model[:P_ij_k][t, l, k] for l in Omega_sending[Ns + i], k in 1:K)

                        [t=1:T, i=1:Nu],    0 ==  
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

    T  = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)
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

    T  = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)
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

    T  = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)
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
                                        ::BIM, ::Convex
                                    )::Nothing

    T  = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)
    lines = network_data.lines

    JuMP.@constraint(model, 
                    [t=1:T, l=1:L],
                    sum(model[:X_ij_k_i][t, l, k, lines[l].edge.from_node.id] for k in 1:K) * 
                    sum(model[:X_ij_k_i][t, l, k, lines[l].edge.to_node.id] for k in 1:K) == 
                    sum(model[:X_ij_k_re][t, l, k] for k in 1:K)^2 + 
                    sum(model[:X_ij_k_im][t, l, k] for k in 1:K)^2)
    return 
end

function _add_PowerFlowConstraints!(model::JuMP.AbstractModel, topology_choice::TopologyChoiceFormulation, ::BIM)::Nothing

    # THINK ABOUT WAYS TO GET THE INDICES OF ALPHA 
    T  = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)
    lines = network_data.lines
    conductors = network_data.conductors
    buses = network_data.conductors

    for l in 1:L 
        ifrom = lines[l].edge.from_node.id
        ito   = lines[l].edge.to_node.id
        for k in 1:K
            Y = 1/(conductors[k].r * lines[l].length + im * (conductors[k].x * lines[l].length))
            G = real(y) 
            B = abs(imag(y))
            for t in 1:T
                alpha_index = compute_index((l, k), t, topology_choice)

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
                                                        model[:Alpha][alpha_index...]

                    model[:X_ij_k_i][t, l, k, ifrom] <= buses[ifrom].V_limits.V_max^2 * 
                                                        model[:Alpha][alpha_index...]

                    model[:X_ij_k_i][t, l, k, ito] >=   buses[ito].V_limits.V_min^2 * 
                                                        model[:Alpha][alpha_index...]

                    model[:X_ij_k_i][t, l, k, ito] <=   buses[ito].V_limits.V_max^2 * 
                                                        model[:Alpha][alpha_index...]

                    model[:X_ij_k_re][t, l, k] <=   buses[ifrom].V_limits.V_max * 
                                                    buses[ito].V_limits.V_max * 
                                                    model[:Alpha][alpha_index...]

                    model[:X_ij_k_im][t, l, k] <=   buses[ifrom].V_limits.V_max * 
                                                    buses[ito].V_limits.V_max * 
                                                    model[:Alpha][alpha_index...]

                    model[:X_ij_k_im][t, l, k] >=   - buses[ifrom].V_limits.V_max * 
                                                    buses[ito].V_limits.V_max * 
                                                    model[:Alpha][alpha_index...]

                    model[:V_sqr][t, ifrom] - model[:X_ij_k_i][t, l, k, ifrom] >= 
                        buses[ifrom].V_limits.V_min^2 * (1 - Alpha[alpha_index...])
    
                    model[:V_sqr][t, ifrom] - model[:X_ij_k_i][t, l, k, ifrom] <= 
                        buses[ifrom].V_limits.V_max^2 * (1 - Alpha[alpha_index...])
        
                    model[:V_sqr][t, ito] - model[:X_ij_k_i][t, l, k, ito] >= 
                        buses[ito].V_limits.V_min^2 * (1 - Alpha[alpha_index...])
        
                    model[:V_sqr][t, ito] - model[:X_ij_k_i][t, l, k, ito] <= 
                            buses[ito].V_limits.V_max^2 * (1 - Alpha[alpha_index...])

                    end)
            end
        end
    end
    
    return 
end


function _add_PowerFlowConstraints!(model::JuMP.AbstractModel, 
                                    topology_choice::TopologyChoiceFormulation,
                                    graph::TypeOfGraph,
                                )::Nothing

    network_data = model[:network_data]
    T  = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)
    conductors = network_data.conductors
    lines = network_data.lines

    voltage_expr = JuMP.@expression(
        model, 
        [t=1:T, l=1:L],
        - sum(2 * (conductors[k].r * len * model[:P_ij_k][t, l, k] + 
        conductors[k].x * len * model[:Q_ij_k][t, l, k]) +
        ((conductors[k].r * len)^2 + (conductors[k].x * len)^2) * 
        model[:I_sqr_k][t, l, k] for k in K)) 

    if isa(graph, Undirected)
        Y = JuMP.@expression(  model, 
                                [t=1:T, l=1:L], 
                                (1 - model[:Y][compute_index((l), t, topology_choice)...]) * 
                                (buses[lines[l].edge.to_node.id].V_limits.V_max^2 - 
                                buses[lines[l].edge.from_node.id].V_limits.V_min^2)
                            )
    elseif isa(graph, Directed)
        Y = JuMP.@expression(  model, 
                                [t=1:T, l=1:L], 
                                (1 - (model[:Y_send][compute_index((l), t, topology_choice)...] + 
                                model[:Y_rec][compute_index((l), t, topology_choice)...])) *
                                (buses[lines[l].edge.to_node.id].V_limits.V_max^2 - 
                                buses[lines[l].edge.from_node.id].V_limits.V_min^2)
                            )
    end

    for l in 1:L 
        ifrom = lines[l].edge.from_node.id
        ito   = lines[l].edge.to_node.id
        len = lines[l].length
        for t in 1:T
            JuMP.@constraints(  model, begin
                model[:V_sqr][t, ito] - model[:V_sqr][t, ifrom] <= voltage_expr[t, l] + 
                                            Y[compute_index((l), t, topology_choice)...]

                model[:V_sqr][t, ito] - model[:V_sqr][t, ifrom] >= voltage_expr[t, l] + 
                                            Y[compute_index((l), t, topology_choice)...]
                
            end)
        end
    end
    return 
end

# ---------------------------------------------------------------------------- #
#                        Conductor choice constraints                          #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#                               Radiality constraints                          #
# ---------------------------------------------------------------------------- #
function _add_RadialityConstraints!(model::JuMP.AbstractModel, 
                                    topology_choice::TopologyChoiceFormulation,
                                    graph_type::TypeOfGraph,
                                    ::RadialityFormulation)::Nothing

    
    return
end