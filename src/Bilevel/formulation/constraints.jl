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
    Ns_init = network_data.nb_init_subs
    # Nsinit and Nsnotinit do this 
    #Ns_init = # nb of buses with S_rating
    V_sqr = model[:V_sqr]
    Beta  = model[:Beta]
    
    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns_init], V_sqr[t, i] == 1
                        [t=1:T, i=(Ns_init+1):Ns], V_sqr[t, i] - 1 <= 
                                        (network_data.buses[i].V_limits.V_max^2 - 1) *
                                        (1 - Beta[i])
                        [t=1:T, i=(Ns_init+1):Ns], V_sqr[t, i] - 1 >= 
                                        (network_data.buses[i].V_limits.V_min^2 - 1) * 
                                        (1 - Beta[i])
                    end)
    return 
end

# ---------------------------------------------------------------------------- #
#                        Load Oversatisfaction constraints                     #
# ---------------------------------------------------------------------------- #

# ------------------------------------ DG ------------------------------------ #
function _add_LoadOverSatisfaction!(model::JuMP.AbstractModel)::Nothing
    network_data = model[:network_data]
    T  = model[:time_steps]
    Ns = network_data.nb_substations
    Nu = network_data.nb_loads
    N  = Ns + Nu
    buses = network_data.buses
    P_sub = model[:P_sub]
    Q_sub = model[:Q_sub]
    p_pv = model[:p_pv]

    JuMP.@constraints(  model, begin
                        [t=1:T], sum(P_sub[t, i] for i in 1:Ns) + 
                                 sum(p_pv[t, i] for i in 1:Nu) >= 
                                 sum(buses[i].load_profile.time_serie[t] .*
                                 buses[i].cos_phi for i in (Ns+1):N)
                        
                        [t=1:T], sum(Q_sub[t, i] for i in 1:Ns) >= 
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

    S_sub = model[:S_sub]
    P_sub = model[:P_sub]
    Q_sub = model[:Q_sub]
    S_sub_capa = model[:S_sub_capa]
    Beta = model[:Beta]

    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], S_sub[t, i]^2 == P_sub[t, i]^2 + Q_sub[t, i]^2
                        [t=1:T, i=1:Ns], S_sub[t, i] <= buses[i].S_rating + S_sub_capa[i] 
                        [i=1:Ns], S_sub_capa[i] <= Beta[i] * (buses[i].S_rating_max - buses[i].S_rating)
                    end)
    return 
end

function _add_SubstationConstraints!(model::JuMP.AbstractModel, ::Convex)::Nothing

    network_data = model[:network_data]
    T  = model[:time_steps]
    Ns = network_data.nb_substations
    buses = network_data.buses

    S_sub = model[:S_sub]
    P_sub = model[:P_sub]
    Q_sub = model[:Q_sub]
    S_sub_capa = model[:S_sub_capa]
    Beta = model[:Beta]

    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], [S_sub[t, i]; P_sub[t, i]; Q_sub[t, i]] in JuMP.SecondOrderCone()
                        [t=1:T, i=1:Ns], S_sub[t, i] <= buses[i].S_rating + S_sub_capa[i]
                        [i=1:Ns], S_sub_capa[i] <= Beta[i] * (buses[i].S_rating_max - buses[i].S_rating)
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

    V_sqr = model[:V_sqr]

    for i in 1:N, t in 1:T
        JuMP.set_lower_bound(V_sqr[t, i], buses[i].V_limits.V_min^2)
    end

    for i in 1:N, t in 1:T
        JuMP.set_upper_bound(V_sqr[t, i], buses[i].V_limits.V_max^2)
    end

    
    return 
end


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

    V_sqr = model[:V_sqr]
    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:N], V_slack_low[t, i] <= 
                                        V_violation_low[t, i] * buses[i].V_limits.V_min^2
                        [t=1:T, i=1:N], V_slack_up[t, i] <= 
                                        V_violation_up[t, i] * buses[i].V_limits.V_max^2
                        [t=1:T, i=1:N], V_sqr[t, i] + V_slack_low[t, i] >= 
                                        buses[i].V_limits.V_min^2
                        [t=1:T, i=1:N], V_sqr[t, i] - V_slack_up[t, i] <= 
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
    I_sqr_k = model[:I_sqr_k]
    
    JuMP.@constraint(model, 
                    [t=1:T, l=1:L, k=1:K],
                    I_sqr_k[t, l, k] <= conductors[k].max_i^2 * 
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
                    I_slack[t, l, k] <= conductors[k].max_i^2 * I_violation[t, l, k])

    I_sqr_k = model[:I_sqr_k]
    JuMP.@constraint(model,
                    [t=1:T, l=1:L, k=1:K],  
                    I_sqr_k[t, l, k] - I_slack[t, l, k] <= conductors[k].max_i^2 * 
                    cond_choice[compute_idx((l, k), t, topology_choice)...])

    return 
end

# ---------------------------------------------------------------------------- #
#                               PowerBalance constraints                       #
# ---------------------------------------------------------------------------- #
function _add_PowerBalanceConstraints!( model::JuMP.AbstractModel,
                                        topology_choice::TopologyChoiceFormulation, 
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

    P_sub = model[:P_sub]
    Q_sub = model[:Q_sub]
    P_ij_k = model[:P_ij_k]
    Q_ij_k = model[:Q_ij_k]
    I_sqr_k = model[:I_sqr_k]
    Alpha = model[:Alpha]
    p_imp = model[:p_imp]
    p_exp = model[:p_exp]
    q_imp = model[:q_imp]
    q_exp = model[:q_exp]

    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], -  P_sub[t, i] ==  
                                            sum((P_ij_k[t, l, k] - 
                                            conductors[k].r * lines[l].length * I_sqr_k[t, l, k])
                                            for l in Omega_receiving[i], k in 1:K) -
                                            sum(P_ij_k[t, l, k] for l in Omega_sending[i], k in 1:K)

                        [t=1:T, i=1:Ns], -  Q_sub[t, i] ==  
                                            sum((Q_ij_k[t, l, k] - 
                                            conductors[k].x * lines[l].length * I_sqr_k[t, l, k])
                                            for l in Omega_receiving[i], k in 1:K) -
                                            sum(Q_ij_k[t, l, k] for l in Omega_sending[i], k in 1:K)

                        [t=1:T, i=1:Nu],    p_imp[t, i] - p_exp[t, i] ==  
                                            sum(P_ij_k[t, l, k] - 
                                            conductors[k].r * lines[l].length * I_sqr_k[t, l, k]
                                            for l in Omega_receiving[Ns + i], k in 1:K) -
                                            sum(P_ij_k[t, l, k] for l in Omega_sending[Ns + i], k in 1:K)

                        [t=1:T, i=1:Nu],    q_imp[t, i] - q_exp[t, i] ==  
                                            sum(Q_ij_k[t, l, k] - 
                                            conductors[k].x * lines[l].length * I_sqr_k[t, l, k]
                                            for l in Omega_receiving[Ns + i], k in 1:K) -
                                            sum(Q_ij_k[t, l, k] for l in Omega_sending[Ns + i], k in 1:K)
                        
                        end
                    )

    if isa(topology_choice, OneConfig)
        Alpha = model[:Alpha]
        JuMP.@constraints(model, begin
            [t=1:T, l=1:L, k=1:K], P_ij_k[t, l, k] <= 2*(conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max) * Alpha[l, k] # indispensable
            [t=1:T, l=1:L, k=1:K], P_ij_k[t, l, k] >= -2*(conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max) * Alpha[l, k] # indispensable
            [t=1:T, l=1:L, k=1:K], Q_ij_k[t, l, k] <= 2*(conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max) * Alpha[l, k] # indispensable
            [t=1:T, l=1:L, k=1:K], Q_ij_k[t, l, k] >= -2*(conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max * Alpha[l, k]) # indispensable
        end)
    elseif isa(topology_choice, ReconfigAllowed)
        Gamma = model[:Gamma]
        JuMP.@constraints(model, begin
            [t=1:T, l=1:L, k=1:K], P_ij_k[t, l, k] <= 2* conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max * Gamma[t, l, k] # indispensable
            [t=1:T, l=1:L, k=1:K], P_ij_k[t, l, k] >= - 2 * conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max * Gamma[t, l, k]  # indispensable
            [t=1:T, l=1:L, k=1:K], Q_ij_k[t, l, k] <= 2 * conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max * Gamma[t, l, k]  # indispensable
            [t=1:T, l=1:L, k=1:K], Q_ij_k[t, l, k] >= - 2* conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max * Gamma[t, l, k]  # indispensable
        end)
    end
    
    return
end

function _add_PowerBalanceConstraints!( model::JuMP.AbstractModel, 
                                        topology_choice::TopologyChoiceFormulation, 
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

    P_sub = model[:P_sub]
    Q_sub = model[:Q_sub]
    P_ij_k = model[:P_ij_k]
    Q_ij_k = model[:Q_ij_k]

    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], -  P_sub[t, i] ==  
                                            sum(P_ji_k[t, l, k] for l in Omega_receiving[i], k in 1:K) -
                                            sum(P_ij_k[t, l, k] for l in Omega_sending[i], k in 1:K)

                        [t=1:T, i=1:Ns], -  Q_sub[t, i] ==  
                                            sum(Q_ji_k[t, l, k] for l in Omega_receiving[i], k in 1:K) -
                                            sum(Q_ij_k[t, l, k] for l in Omega_sending[i], k in 1:K)

                        [t=1:T, i=1:Nu],    p_imp[t, i] - p_exp[t, i] ==  
                                            sum(P_ji_k[t, l, k] for l in Omega_receiving[Ns + i], k in 1:K) -
                                            sum(P_ij_k[t, l, k] for l in Omega_sending[Ns + i], k in 1:K)

                        [t=1:T, i=1:Nu],    q_imp[t, i] - q_exp[t, i] ==  
                                            sum(Q_ji_k[t, l, k] for l in Omega_receiving[Ns + i], k in 1:K) -
                                            sum(Q_ij_k[t, l, k] for l in Omega_sending[Ns + i], k in 1:K)
                    end)

    if isa(topology_choice, OneConfig)
        Alpha = model[:Alpha]
        JuMP.@constraints(model, begin
            [t=1:T, l=1:L, k=1:K], P_ij_k[t, l, k] <= conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max * Alpha[l, k] # indispensable
            [t=1:T, l=1:L, k=1:K], P_ji_k[t, l, k] <= conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max * Alpha[l, k] # indispensable
            [t=1:T, l=1:L, k=1:K], Q_ij_k[t, l, k] <= conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max * Alpha[l, k] # indispensable
            [t=1:T, l=1:L, k=1:K], Q_ji_k[t, l, k] <= conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max * Alpha[l, k] # indispensable
        end)
    elseif isa(topology_choice, ReconfigAllowed)
        Gamma = model[:Gamma]
        JuMP.@constraints(model, begin
            [t=1:T, l=1:L, k=1:K], P_ij_k[t, l, k] <= conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max * Gamma[t, l, k] # indispensable
            [t=1:T, l=1:L, k=1:K], P_ij_k[t, l, k] <= conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max * Gamma[t, l, k]  # indispensable
            [t=1:T, l=1:L, k=1:K], Q_ij_k[t, l, k] <= conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max * Gamma[t, l, k]  # indispensable
            [t=1:T, l=1:L, k=1:K], Q_ij_k[t, l, k] <= conductors[k].max_i * buses[lines[l].edge.from_node.id].V_limits.V_max * Gamma[t, l, k]  # indispensable
        end)
    end          
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

    V_sqr = model[:V_sqr]
    I_sqr = model[:I_sqr]
    P_ij = model[:P_ij]
    Q_ij = model[:Q_ij]

    JuMP.@constraint(model, 
                    [t=1:T, l=1:L],
                    [V_sqr[t, lines[l].edge.from_node.id] / 2; 
                    I_sqr[t, l]; 
                    P_ij[t, l]; 
                    Q_ij[t, l]] in 
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

    V_sqr = model[:V_sqr]
    I_sqr_k = model[:I_sqr_k]
    P_ij_k = model[:P_ij_k]
    Q_ij_k = model[:Q_ij_k]

    JuMP.@constraint(model,
                    [t=1:T, l=1:L],
                    V_sqr[t, lines[l].edge.from_node.id] * 
                    I_sqr[t, l] == 
                    P_ij[t, l]^2 + 
                    Q_ij[t, l]^2)
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

    X_ij_k_i = model[:X_ij_k_i]
    X_ij_k_re = model[:X_ij_k_re]
    X_ij_k_im = model[:X_ij_k_im]
  
    JuMP.@constraint(model, 
                    [t=1:T, l=1:L],
                    [sum(X_ij_k_i[t, l, k, lines[l].edge.from_node.id] for k in 1:K) / 2; 
                    sum(X_ij_k_i[t, l, k, lines[l].edge.to_node.id] for k in 1:K); 
                    sum(X_ij_k_re[t, l, k] for k in 1:K); 
                    sum(X_ij_k_im[t, l, k] for k in 1:K)] in 
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

    X_ij_k_i = model[:X_ij_k_i]
    X_ij_k_re = model[:X_ij_k_re]
    X_ij_k_im = model[:X_ij_k_im]

    JuMP.@constraint(model, 
                    [t=1:T, l=1:L],
                    sum(X_ij_k_i[t, l, k, lines[l].edge.from_node.id] for k in 1:K) * 
                    sum(X_ij_k_i[t, l, k, lines[l].edge.to_node.id] for k in 1:K) == 
                    sum(X_ij_k_re[t, l, k] for k in 1:K)^2 + 
                    sum(X_ij_k_im[t, l, k] for k in 1:K)^2)
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

    P_ij_k = model[:P_ij_k]
    P_ji_k = model[:P_ji_k]
    Q_ij_k = model[:Q_ij_k]
    Q_ji_k = model[:Q_ji_k]
    I_sqr_k = model[:I_sqr_k]
    X_ij_k_i = model[:X_ij_k_i]
    X_ij_k_re = model[:X_ij_k_re]
    X_ij_k_im = model[:X_ij_k_im]
    V_sqr = model[:V_sqr]

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
                                
                    P_ij_k[t, l, k] ==  G * (X_ij_k_i[t, l, k, ifrom] - 
                                                X_ij_k_re[t, l, k]) + 
                                                B * X_ij_k_im[t, l, k]

                    P_ji_k[t, l, k] ==  G * (X_ij_k_i[t, l, k, ito] - 
                                                X_ij_k_re[t, l, k]) - 
                                                B * X_ij_k_im[t, l, k]

                    Q_ij_k[t, l, k] ==  B * (X_ij_k_i[t, l, k, ifrom] - 
                                                X_ij_k_re[t, l, k]) - 
                                                G * X_ij_k_im[t, l, k]

                    Q_ji_k[t, l, k] ==  B * (X_ij_k_i[t, l, k, ifrom] - 
                                                X_ij_k_re[t, l, k]) + 
                                                G * X_ij_k_im[t, l, k]

                    I_sqr_k[t, l, k] == (G^2 + B^2) * (X_ij_k_i[t, l, k, ifrom] + 
                                                X_ij_k_i[t, l, k, ito] - 
                                                2 * X_ij_k_re[t, l, k])

                    X_ij_k_i[t, l, k, ifrom] >= buses[ifrom].V_limits.V_min^2 * 
                                                        cond_choice[index...]

                    X_ij_k_i[t, l, k, ifrom] <= buses[ifrom].V_limits.V_max^2 * 
                                                        cond_choice[index...]

                    X_ij_k_i[t, l, k, ito] >=   buses[ito].V_limits.V_min^2 * 
                                                        cond_choice[index...]

                    X_ij_k_i[t, l, k, ito] <=   buses[ito].V_limits.V_max^2 * 
                                                        cond_choice[index...]

                    X_ij_k_re[t, l, k] <=   buses[ifrom].V_limits.V_max * 
                                                    buses[ito].V_limits.V_max * 
                                                    cond_choice[index...]

                    X_ij_k_im[t, l, k] <=   buses[ifrom].V_limits.V_max * 
                                                    buses[ito].V_limits.V_max * 
                                                    cond_choice[index...]

                    X_ij_k_im[t, l, k] >=   - buses[ifrom].V_limits.V_max * 
                                                    buses[ito].V_limits.V_max * 
                                                    cond_choice[index...]

                    V_sqr[t, ifrom] - X_ij_k_i[t, l, k, ifrom] >= 
                        buses[ifrom].V_limits.V_min^2 * (1 - cond_choice[index...])
    
                    V_sqr[t, ifrom] - X_ij_k_i[t, l, k, ifrom] <= 
                        buses[ifrom].V_limits.V_max^2 * (1 - cond_choice[index...])
        
                    V_sqr[t, ito] - X_ij_k_i[t, l, k, ito] >= 
                        buses[ito].V_limits.V_min^2 * (1 - cond_choice[index...])
        
                    V_sqr[t, ito] - X_ij_k_i[t, l, k, ito] <= 
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

    P_ij_k = model[:P_ij_k]
    Q_ij_k = model[:Q_ij_k]
    I_sqr_k = model[:I_sqr_k]
    voltage_expr = JuMP.@expression(
        model, 
        [t=1:T, l=1:L], sum( -2 * (conductors[k].r * lines[l].length * P_ij_k[t, l, k] + 
        conductors[k].x * lines[l].length * Q_ij_k[t, l, k]) +
        ((conductors[k].r * lines[l].length)^2 + (conductors[k].x * lines[l].length)^2) * 
        I_sqr_k[t, l, k] for k in 1:K)) 

    if isa(graph, Undirected)
        Y = model[:Y]
        Y_expr = JuMP.@expression(  model, 
                                [t=1:T, l=1:L], 
                                (1 - Y[compute_idx((l,), t, topology_choice)...]) * 
                                (buses[lines[l].edge.to_node.id].V_limits.V_max^2 - 
                                buses[lines[l].edge.from_node.id].V_limits.V_min^2)
                            )

    elseif isa(graph, Directed)
        Y_send = model[:Y_send]
        Y_rec = model[:Y_rec]

        Y_expr = JuMP.@expression(  model, 
                                [t=1:T, l=1:L], 
                                (1 - (Y_send[compute_idx((l,), t, topology_choice)...] + 
                                Y_rec[compute_idx((l,), t, topology_choice)...])) *
                                (buses[lines[l].edge.to_node.id].V_limits.V_max^2 - 
                                buses[lines[l].edge.from_node.id].V_limits.V_min^2)
                            )
    end

    V_sqr = model[:V_sqr]
    for l in 1:L 
        ifrom = lines[l].edge.from_node.id
        ito   = lines[l].edge.to_node.id
        for t in 1:T
            JuMP.@constraints(  model, begin
                V_sqr[t, ito] - V_sqr[t, ifrom] <= voltage_expr[t, l] + 
                                            Y_expr[compute_idx((l,), t, topology_choice)...]

                V_sqr[t, ito] - V_sqr[t, ifrom] >= voltage_expr[t, l] - 
                                            Y_expr[compute_idx((l,), t, topology_choice)...]
                
            end)
        end
    end
    return 
end

# ---------------------------------------------------------------------------- #
#                               Lower Level constraints                        #
# ---------------------------------------------------------------------------- #
# I AM HERE : ERROR WITH PV_PROD
function _add_LowerConstraints!(model::JuMP.AbstractModel)::Nothing
    network_data = model[:network_data]
    T  = model[:time_steps]
    Nu = network_data.nb_loads
    Ns = network_data.nb_substations
    N = Ns + Nu
    buses = network_data.buses 

    PV_prod = [(isnothing(buses[Ns + i].PV_installation) ? 0.0 : buses[Ns + i].PV_installation.profile.time_serie[t]) for t in 1:T, i in 1:Nu]

    P_consumed = [buses[i].load_profile.time_serie[t] * buses[i].cos_phi for t in 1:T, i in (Ns+1):N]
    Q_consumed = [buses[i].load_profile.time_serie[t] * sin(acos(buses[i].cos_phi)) for t in 1:T, i in (Ns+1):N]

    p_pv = model[:p_pv]
    s_conv_pv = model[:s_conv_pv]
    p_pv_max = model[:p_pv_max]
    p_imp = model[:p_imp]
    p_exp = model[:p_exp]
    q_imp = model[:q_imp]
    q_exp = model[:q_exp]
    s_grid_max = model[:s_grid_max]

    # Power balance
    power_balance = JuMP.@expression(model, [t=1:T, i=1:Nu], P_consumed[t, i] - p_pv[t, i])

    if model[:storage]
        # Fetching the required data
        storage_capacity = model[:storage_capacity]
        storage_state = model[:storage_state]
        p_storage_charge = model[:p_storage_charge]
        p_storage_discharge = model[:p_storage_discharge]
        NB_PROFILES = model[:nb_sign_days]
        STORAGE_EFF = [buses[i].storage.efficiency for i in Ns+1:N]
        DELTA_T = model[:delta_t]


        # Storage constraints
        JuMP.@constraint(model, [t=1:T, i=1:Nu], p_imp[t, i] - p_exp[t, i] == P_consumed[t, i] - p_pv[t, i] + p_storage_charge[t, i] - p_storage_discharge[t, i])

        @constraint(model,  [t=1:T, i=1:Nu], storage_state[t, i] <= storage_capacity[i])

        for season in 1:NB_PROFILES
            idx_first = Int((season - 1) * T / NB_PROFILES + 1)
            idx_last = Int(season * T  / NB_PROFILES)

            @constraints(model, begin
                [t=idx_first+1:idx_last, i=1:Nu], storage_state[t, i] == storage_state[t-1, i] + (p_storage_charge[t, i] * STORAGE_EFF[i] - (1/(STORAGE_EFF[i]- 0.05)) * p_storage_discharge[t, i]) * DELTA_T/60
                # Boundary effects
                [i=1:Nu], storage_state[idx_first, i] == 0.1*storage_capacity[i]
                [i=1:Nu], storage_state[idx_first, i] == storage_state[idx_last, i] + (p_storage_charge[idx_first, i] * STORAGE_EFF[i] - (1/(STORAGE_EFF[i]- 0.05)) * p_storage_discharge[idx_first, i]) * DELTA_T/60
            end)
        end
    else
        JuMP.@constraint(model, [t=1:T, i=1:Nu], p_imp[t, i] - p_exp[t, i] == P_consumed[t, i] - p_pv[t, i])
    end

    
    @constraints(model, begin
        [t=1:T, i=1:Nu], q_imp[t, i] - q_exp[t, i] == Q_consumed[t, i]
        [t=1:T, i=1:Nu], p_imp[t, i] <= s_grid_max[i]
        [t=1:T, i=1:Nu], q_imp[t, i] <= s_grid_max[i]
        [t=1:T, i=1:Nu], p_exp[t, i] <= s_grid_max[i]
        [t=1:T, i=1:Nu], q_exp[t, i] <= s_grid_max[i]
        [t=1:T, i=1:Nu], p_pv[t, i] <= s_conv_pv[i] 
        [t=1:T, i=1:Nu], p_pv[t, i] <= PV_prod[t, i] * p_pv_max[i] # T
        [t=1:T, i=1:Nu], p_exp[t, i] <= p_pv[t, i]
    end)
    return
end

# ---------------------------------------------------------------------------- #
#                               Radiality constraints                          #
# ---------------------------------------------------------------------------- #

# -------------------------- Single Commodity Flow --------------------------- #

function _add_RadialityConstraints!(model::JuMP.AbstractModel, 
                                    graph_type::TypeOfGraph,
                                    ::OneConfig,
                                    ::SingleCommodityFlow)::Nothing

    Omega_sending = model[:network_topology].sending_lines
    Omega_receiving = model[:network_topology].receiving_lines
    network = model[:network_data]
    L = network.nb_lines
    Ns_init = network.nb_init_subs
    Ns = network.nb_substations
    Nu = network.nb_loads
    N = Ns + Nu
    Beta = model[:Beta]
    k_ij = model[:k_ij]

    # idea : sum of all the flows sent by substation = Nu
    # flow in the two directions ?
    JuMP.@constraints(model, begin
                    [i=(Ns+1):N], - sum(k_ij[l] for l in Omega_receiving[i]) + 
                                sum(k_ij[l] for l in Omega_sending[i]) == -1
                    [i=1:Ns],   - sum(k_ij[l] for l in Omega_receiving[i]) + 
                                sum(k_ij[l] for l in Omega_sending[i]) >= 0
                    [i=1:Ns_init],  - sum(k_ij[l] for l in Omega_receiving[i]) +
                                    sum(k_ij[l] for l in Omega_sending[i]) <= Nu
                    [i=(Ns_init+1):Ns], - sum(k_ij[l] for l in Omega_receiving[i]) +
                                        sum(k_ij[l] for l in Omega_sending[i]) <= Nu * Beta[i]
                        
    end)

    if isa(graph_type, Undirected)
        Y = model[:Y]
        JuMP.@constraints(model, begin
            [l=1:L], k_ij[l] <= Nu * Y[l]               
            [l=1:L], k_ij[l] >= - Nu * Y[l]
        end)

    elseif isa(graph_type, Directed)
        Y_send = model[:Y_send]
        Y_rec = model[:Y_rec]
        JuMP.@constraints(model, begin
            [l=1:L], k_ij[l] <= Nu * (Y_send[l] + Y_rec[l])               
            [l=1:L], k_ij[l] >= - Nu * (Y_send[l] + Y_rec[l])     
        end)
    end

    return
end

function _add_RadialityConstraints!(model::JuMP.AbstractModel, 
                                    graph_type::TypeOfGraph,
                                    ::ReconfigAllowed,
                                    ::SingleCommodityFlow)::Nothing

    Omega_sending = model[:network_topology].sending_lines
    Omega_receiving = model[:network_topology].receiving_lines
    network = model[:network_data]
    T = model[:time_steps]
    L = network.nb_lines
    Ns_init = network.nb_init_subs
    Ns = network.nb_substations
    Nu = network.nb_loads
    N = Ns + Nu
    Beta = model[:Beta]
    k_ij = model[:k_ij]

    JuMP.@constraints(model, begin
                    [t=1:T, i=(Ns+1):N], - sum(k_ij[t, l] for l in Omega_receiving[i]) +
                                        sum(k_ij[t, l] for l in Omega_sending[i]) == -1
                    [t=1:T, i=1:Ns],   - sum(k_ij[t, l] for l in Omega_receiving[i]) +
                                sum(k_ij[t, l] for l in Omega_sending[i]) >= 0
                    [t=1:T, i=1:Ns_init],  - sum(k_ij[t, l] for l in Omega_receiving[i]) +
                                     sum(k_ij[t, l] for l in Omega_sending[i]) <= Nu
                    [t=1:T, i=(Ns_init+1):Ns], - sum(k_ij[t, l] for l in Omega_receiving[i]) +
                                        sum(k_ij[t, l] for l in Omega_sending[i]) <= Nu * Beta[i]
                        
    end)

    if isa(graph_type, Undirected)
        Y = model[:Y]
        JuMP.@constraints(model, begin
            [t=1:T, l=1:L], k_ij[t, l] <= Nu * (Y[t, l])               
            [t=1:T, l=1:L], k_ij[t, l] >= - Nu * (Y[t, l])
        end)
        
    elseif isa(graph_type, Directed)
        Y_send = model[:Y_send]
        Y_rec = model[:Y_rec]
        JuMP.@constraints(model, begin
            [t=1:T, l=1:L], k_ij[t, l] <= Nu * (Y_send[t, l] + Y_rec[t, l])               
            [t=1:T, l=1:L], k_ij[t, l] >= - Nu * (Y_send[t, l] + Y_rec[t, l])     
        end)
    end

    return
end

# -------------------------- Multi Commodity Flow --------------------------- #
function _add_RadialityConstraints!(model::JuMP.AbstractModel, 
                                    graph_type::TypeOfGraph,
                                    ::OneConfig,
                                    ::MultiCommodityFlow)::Nothing

    Omega_sending = model[:network_topology].sending_lines
    Omega_receiving = model[:network_topology].receiving_lines
    network = model[:network_data]
    T = model[:time_steps]
    L = network.nb_lines
    Ns_init = network.nb_init_subs
    Ns = network.nb_substations
    Nu = network.nb_loads
    N = Ns + Nu
    Beta = model[:Beta]
    k_ij = model[:k_ij]

    
    JuMP.@constraints(model, begin
                    [i=1:Ns, w=(Ns+1):N], - sum(k_ij[l, w] for l in Omega_receiving[i]) + 
                                            sum(k_ij[l, w] for l in Omega_sending[i]) >= 0
                    [i=1:Ns_init, w=(Ns+1):N], - sum(k_ij[l, w] for l in Omega_receiving[i]) +
                                                sum(k_ij[l, w] for l in Omega_sending[i]) <= 1
                    [i=(Ns_init+1):Ns, w=(Ns+1):N],  -  sum(k_ij[l, w] for l in Omega_receiving[i]) +
                                                        sum(k_ij[l, w] for l in Omega_sending[i]) <= Beta[i]
                    [i=(Ns+1):N], - sum(k_ij[l, i] for l in Omega_receiving[i]) + 
                                    sum(k_ij[l, i] for l in Omega_sending[i]) == -1
                    [i=(Ns+1):N, w=(Ns+1):N; i != w], - sum(k_ij[l, w] for l in Omega_receiving[i]) +
                                                        sum(k_ij[l, w] for l in Omega_sending[i]) == 0
    end)
   

    if isa(graph_type, Undirected)
        Y = model[:Y]
        JuMP.@constraints(model, begin
            [l=1:L, w=(Ns+1):N], k_ij[l, w] <= Y[l]
            [l=1:L, w=(Ns+1):N], k_ij[l, w] >= -Y[l]
        end)

    elseif isa(graph_type, Directed)
        Y_send = model[:Y_send]
        Y_rec = model[:Y_rec]
        JuMP.@constraints(model, begin
            [l=1:L, w=(Ns+1):N], k_ij[l, w] <= Y_send[l]
            [l=1:L, w=(Ns+1):N], k_ij[l, w] >= -Y_rec[l] 
        end)
    end

    return
end

function _add_RadialityConstraints!(model::JuMP.AbstractModel, 
                                    graph_type::TypeOfGraph,
                                    ::ReconfigAllowed,
                                    ::MultiCommodityFlow)::Nothing

    Omega_sending = model[:network_topology].sending_lines
    Omega_receiving = model[:network_topology].receiving_lines
    network = model[:network_data]
    T = model[:time_steps]
    L = network.nb_lines
    Ns_init = network.nb_init_subs
    Ns = network.nb_substations
    Nu = network.nb_loads
    N = Ns + Nu
    Beta = model[:Beta]
    k_ij = model[:k_ij]


    JuMP.@constraints(model, begin
        [t=1:T, i=1:Ns, w=(Ns+1):N], - sum(k_ij[t, l, w] for l in Omega_receiving[i]) + 
                    sum(k_ij[t, l, w] for l in Omega_sending[i]) >= 0
        [t=1:T, i=1:Ns_init, w=(Ns+1):N], - sum(k_ij[t, l, w] for l in Omega_receiving[i]) +
                        sum(k_ij[t, l, w] for l in Omega_sending[i]) <= 1
        [t=1:T, i=(Ns_init+1):Ns, w=(Ns+1):N],  -  sum(k_ij[t, l, w] for l in Omega_receiving[i]) +
                                sum(k_ij[t, l, w] for l in Omega_sending[i]) <= Beta[i]
        [t=1:T, i=(Ns+1):N], - sum(k_ij[t, l, i] for l in Omega_receiving[i]) + 
            sum(k_ij[t, l, i] for l in Omega_sending[i]) == -1
        [t=1:T, i=(Ns+1):N, w=(Ns+1):N; i != w], - sum(k_ij[t, l, w] for l in Omega_receiving[i]) +
                                sum(k_ij[t, l, w] for l in Omega_sending[i]) == 0
    end)


    if isa(graph_type, Undirected)
        Y = model[:Y]
        JuMP.@constraints(model, begin
            [t=1:T, l=1:L, w=(Ns+1):N], k_ij[t, l, w] <= Y[t, l]
            [t=1:T, l=1:L, w=(Ns+1):N], k_ij[t, l, w] >= -Y[t, l]
    end)

    elseif isa(graph_type, Directed)
        Y_send = model[:Y_send]
        Y_rec = model[:Y_rec]
        JuMP.@constraints(model, begin
            [t=1:T, l=1:L, w=(Ns+1):N], k_ij[t, l, w] <= Y_send[t, l]
            [t=1:T, l=1:L, w=(Ns+1):N], k_ij[t, l, w] >= -Y_rec[t, l] 
        end)
    end

    return
end