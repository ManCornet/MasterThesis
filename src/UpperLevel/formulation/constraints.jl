# ---------------------------------------------------------------------------- #
#                         Voltage Reference Constraints                        #
# ---------------------------------------------------------------------------- #
# 1 pu reference at each substation if its is built 
# don't forget to think also to add the fact that some substations can already
# be there
function _add_RefVoltages!(model::JuMP.AbstractModel)::Nothing
    network_data = model[:network_data]
    T  = model[:time_steps]
    Ns = get_nb_sub_bus(network_data)
    println("helloeeeee")
    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], model[:V_sqr][t, i] - 1 <= 
                                        (network.sub_buses[i].V_limits.V_max^2 - 1) *
                                        (1 - model[:Beta][i])
                        [t=1:T, i=1:Ns], model[:V_sqr][t, i] - 1 >= 
                                        (network.sub_buses[i].V_limits.V_min^2 - 1) * 
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
    Nu = get_nb_load_bus(network_data)
    Ns = get_nb_sub_bus(network_data)
    load_buses = network_data.load_buses

    JuMP.@constraints(  model, begin
                        [t=1:T], sum(model[:P_sub][t, i] for i in 1:Ns) >= 
                                 sum(load_buses[i].load_profile.time_serie[t] .*
                                 load_buses[i].cos_phi for i in 1:Nu)
                        
                        [t=1:T], sum(model[:Q_sub][t, i] for i in 1:Ns) >= 
                                 sum(load_buses[i].load_profile.time_serie[t] .* 
                                 sin(acos(load_buses[i].cos_phi)) for i in 1:Nu)
                    end)
    return 
end
# ------------------------------------ DG ------------------------------------ #
function _add_LoadOverSatisfaction!(model::JuMP.AbstractModel, ::DG)::Nothing
    network_data = model[:network_data]
    T  = model[:time_steps]
    Nu = get_nb_load_bus(network_data)
    Ns = get_nb_sub_bus(network_data)
    load_buses = network_data.load_buses

    JuMP.@constraints(  model, begin
                        [t=1:T], sum(model[:P_sub][t, i] for i in 1:Ns) + 
                                 sum(model[:p_pv][t, i] for i in 1:Nu) >= 
                                 sum(load_buses[i].load_profile.time_serie[t] .*
                                 load_buses[i].cos_phi for i in 1:Nu)
                        
                        [t=1:T], sum(model[:Q_sub][t, i] for i in 1:Ns) >= 
                                 sum(load_buses[i].load_profile.time_serie[t] .* 
                                 sin(acos(load_buses[i].cos_phi)) for i in 1:Nu)
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
    Ns = get_nb_sub_bus(network_data)
    sub_buses = network_data.sub_buses

    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], model[:S_sub][t, i]^2 == model[:P_sub][t, i]^2 + model[:Q_sub][t, i]^2
                        [t=1:T, i=1:Ns], model[:S_sub][t, i] <= model[:S_sub_capa][i]
                        [i=1:Ns], model[:S_sub_capa][i] <= model[:Beta][i] * sub_buses[i].S_rating_max
                    end)
    return 
end

function _add_SubstationConstraints!(model::JuMP.AbstractModel, ::Convex)::Nothing

    network_data = model[:network_data]
    T  = model[:time_steps]
    Ns = get_nb_sub_bus(network_data)
    sub_buses = network_data.sub_buses

    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], [model[:S_sub][t, i]; model[:P_sub][t, i]; model[:Q_sub][t, i]] in JuMP.SecondOrderCone()
                        [t=1:T, i=1:Ns], model[:S_sub][t, i] <= model[:S_sub_capa][i]
                        [i=1:Ns], model[:S_sub_capa][i] <= model[:Beta][i] * sub_buses[i].S_rating_max
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
    Ns = get_nb_sub_bus(network_data)
    Nu = get_nb_load_bus(network_data)
    load_buses = network_data.load_buses
    sub_buses = network_data.sub_buses

    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], model[:V_sqr][t, i] >= sub_buses[i].V_limits.V_min^2
                        [t=1:T, i=1:Ns], model[:V_sqr][t, i] <= sub_buses[i].V_limits.V_max^2
                        [t=1:T, i=1:Nu], model[:V_sqr][t, Ns + i] >= load_buses[i].V_limits.V_min^2
                        [t=1:T, i=1:Nu], model[:V_sqr][t, Ns + i] <= load_buses[i].V_limits.V_max^2
                    end)
    return 
end

################## TO DO ##############################
#######################################################
function _add_VoltageOpConstraints!(model::JuMP.AbstractModel, ::RelaxedVoltages)::Nothing

    network_data = model[:network_data]
    T  = model[:time_steps]
    Ns = get_nb_sub_bus(network_data)
    Nu = get_nb_load_bus(network_data)
    load_buses = network_data.load_buses
    sub_buses = network_data.sub_buses

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
function _add_CurrentOpConstraints!(model::JuMP.AbstractModel, ::StrongCurrents, ::OneConfig)::Nothing

    network_data = model[:network_data]
    T  = model[:time_steps]
    L  = get_nb_lines(network_data)
    K  = get_nb_conductors(network_data)
    conductors = network_data.conductors

    JuMP.@constraint(model, [t=1:T, l=1:L, k=1:K], I_sqr_k[t, l, k] <= conductors[k].max_i^2 * Alpha[l, k])
    return 
end

function _add_CurrentOpConstraints!(model::JuMP.AbstractModel, ::StrongCurrents, ::ReconfigAllowed)::Nothing

    network_data = model[:network_data]
    T  = model[:time_steps]
    L  = get_nb_lines(network_data)
    K  = get_nb_conductors(network_data)
    conductors = network_data.conductors

    JuMP.@constraint(model, [t=1:T, l=1:L, k=1:K], I_sqr_k[t, l, k] <= conductors[k].max_i^2 * Alpha[t, l, k])
    return 
end


function _add_CurrentOpConstraints!(model::JuMP.AbstractModel, ::RelaxedCurrents, ::OneConfig)::Nothing

    network_data = model[:network_data]
    T  = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)
    conductors = network_data.conductors
    

    JuMP.@variables( model, begin 
                    I_violation[1:T, 1:L, 1:K], (binary=true)
                    I_slack[1:T, 1:L, 1:K]
                    end)

    JuMP.@constraints(  model, begin
                        [t=1:T, l=1:L, k=1:K], model[:I_slack][t, l, k] <= conductors[k].max_i^2 * model[:I_violation][t, l, k]
                        [t=1:T, l=1:L, k=1:K], model[:I_sqr_k][t, l, k] - model[:I_slack][t, l, k] <= conductors[k].max_i^2 * model[:Alpha][l, k] # indispensable mais relaxée
                    end)
    return 
end

function _add_CurrentOpConstraints!(model::JuMP.AbstractModel, ::RelaxedCurrents, ::ReconfigAllowed)::Nothing

    network_data = model[:network_data]
    T  = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)
    conductors = network_data.conductors
    

    JuMP.@variables( model, begin 
                    I_violation[1:T, 1:L, 1:K], (binary=true)
                    I_slack[1:T, 1:L, 1:K]
                    end)

    JuMP.@constraints(  model, begin
                        [t=1:T, l=1:L, k=1:K], model[:I_slack][t, l, k] <= conductors[k].max_i^2 * model[:I_violation][t, l, k]
                        [t=1:T, l=1:L, k=1:K], model[:I_sqr_k][t, l, k] - model[:I_slack][t, l, k] <= conductors[k].max_i^2 * model[:Alpha][t, l, k] # indispensable mais relaxée
                    end)
    return 
end
# ---------------------------------------------------------------------------- #
#                               PowerFlow constraints                          #
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
#                               PowerBalance constraints                       #
# ---------------------------------------------------------------------------- #
function _add_PowerBalance!(model::JuMP.AbstractModel, ::BFM, ::NoDG)::Nothing

    network_data = model[:network_data]
    T  = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)
    Ns = get_nb_sub_bus(network_data)
    Nu = get_nb_load_bus(network_data)
    conductors = network_data.conductors
    lines = network_data.lines

    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], -  model[:P_sub][t, i] ==  
                                            sum(model[:P_ij_k][t, l, k] - 
                                            conductors[k].r * lines[l].length * model[:I_sqr_k][t, l, k]
                                            for l in Omega_receiving[i], k in K) -
                                            sum(model[:P_ij_k][t, l, k] for l in Omega_sending[i], k in K)

                        [t=1:T, i=1:Ns], -  model[:Q_sub][t, i] ==  
                                            sum(model[:Q_ij_k][t, l, k] - 
                                            conductors[k].x * lines[l].length * model[:I_sqr_k][t, l, k]
                                            for l in Omega_receiving[i], k in K) -
                                            sum(model[:Q_ij_k][t, l, k] for l in Omega_sending[i], k in K)

                        [t=1:T, i=1:Nu],    0 ==  
                                            sum(model[:P_ij_k][t, l, k] - 
                                            conductors[k].r * lines[l].length * model[:I_sqr_k][t, l, k]
                                            for l in Omega_receiving[Ns + i], k in K) -
                                            sum(model[:P_ij_k][t, l, k] for l in Omega_sending[Ns + i], k in K)

                        [t=1:T, i=1:Nu],    0 ==  
                                            sum(model[:Q_ij_k][t, l, k] - 
                                            conductors[k].x * lines[l].length * model[:I_sqr_k][t, l, k]
                                            for l in Omega_receiving[Ns + i], k in K) -
                                            sum(model[:Q_ij_k][t, l, k] for l in Omega_sending[Ns + i], k in K)
                    end)
    return 
end

function _add_PowerBalance!(model::JuMP.AbstractModel, ::BFM, ::DG)::Nothing

    network_data = model[:network_data]
    T  = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)
    Ns = get_nb_sub_bus(network_data)
    Nu = get_nb_load_bus(network_data)
    conductors = network_data.conductors
    lines = network_data.lines

    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], -  model[:P_sub][t, i] ==  
                                            sum(model[:P_ij_k][t, l, k] - 
                                            conductors[k].r * lines[l].length * model[:I_sqr_k][t, l, k]
                                            for l in Omega_receiving[i], k in K) -
                                            sum(model[:P_ij_k][t, l, k] for l in Omega_sending[i], k in K)

                        [t=1:T, i=1:Ns], -  model[:Q_sub][t, i] ==  
                                            sum(model[:Q_ij_k][t, l, k] - 
                                            conductors[k].x * lines[l].length * model[:I_sqr_k][t, l, k]
                                            for l in Omega_receiving[i], k in K) -
                                            sum(model[:Q_ij_k][t, l, k] for l in Omega_sending[i], k in K)

                        [t=1:T, i=1:Nu], -  model[:p_pv][i] ==  
                                            sum(model[:P_ij_k][t, l, k] - 
                                            conductors[k].r * lines[l].length * model[:I_sqr_k][t, l, k]
                                            for l in Omega_receiving[Ns + i], k in K) -
                                            sum(model[:P_ij_k][t, l, k] for l in Omega_sending[Ns + i], k in K)

                        [t=1:T, i=1:Nu],    0 ==  
                                            sum(model[:Q_ij_k][t, l, k] - 
                                            conductors[k].x * lines[l].length * model[:I_sqr_k][t, l, k]
                                            for l in Omega_receiving[Ns + i], k in K) -
                                            sum(model[:Q_ij_k][t, l, k] for l in Omega_sending[Ns + i], k in K)
                    end)
    return 
end

function _add_PowerBalance!(model::JuMP.AbstractModel, ::BIM, ::NoDG)::Nothing

    network_data = model[:network_data]
    T  = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)
    Ns = get_nb_sub_bus(network_data)
    Nu = get_nb_load_bus(network_data)
    conductors = network_data.conductors
    lines = network_data.lines

    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], -  model[:P_sub][t, i] ==  
                                            sum(model[:P_ji_k][t, l, k] for l in Omega_receiving[i], k in K) -
                                            sum(model[:P_ij_k][t, l, k] for l in Omega_sending[i], k in K)

                        [t=1:T, i=1:Ns], -  model[:Q_sub][t, i] ==  
                                            sum(model[:Q_ji_k][t, l, k] for l in Omega_receiving[i], k in K) -
                                            sum(model[:Q_ij_k][t, l, k] for l in Omega_sending[i], k in K)

                        [t=1:T, i=1:Nu],    0 ==  
                                            sum(model[:P_ji_k][t, l, k] for l in Omega_receiving[Ns + i], k in K) -
                                            sum(model[:P_ij_k][t, l, k] for l in Omega_sending[Ns + i], k in K)

                        [t=1:T, i=1:Nu],    0 ==  
                                            sum(model[:Q_ji_k][t, l, k] for l in Omega_receiving[i], k in K) -
                                            sum(model[:Q_ij_k][t, l, k] for l in Omega_sending[i], k in K)
                    end)
    return 
end

function _add_PowerBalance!(model::JuMP.AbstractModel, ::BIM, ::DG)::Nothing

    network_data = model[:network_data]
    T  = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)
    Ns = get_nb_sub_bus(network_data)
    Nu = get_nb_load_bus(network_data)
    conductors = network_data.conductors
    lines = network_data.lines

    JuMP.@constraints(  model, begin
                        [t=1:T, i=1:Ns], -  model[:P_sub][t, i] ==  
                                            sum(model[:P_ji_k][t, l, k] for l in Omega_receiving[i], k in K) -
                                            sum(model[:P_ij_k][t, l, k] for l in Omega_sending[i], k in K)

                        [t=1:T, i=1:Ns], -  model[:Q_sub][t, i] ==  
                                            sum(model[:Q_ji_k][t, l, k] for l in Omega_receiving[i], k in K) -
                                            sum(model[:Q_ij_k][t, l, k] for l in Omega_sending[i], k in K)

                        [t=1:T, i=1:Nu], -  model[:p_pv][t, i]  ==  
                                            sum(model[:P_ji_k][t, l, k] for l in Omega_receiving[Ns + i], k in K) -
                                            sum(model[:P_ij_k][t, l, k] for l in Omega_sending[Ns + i], k in K)

                        [t=1:T, i=1:Nu],    0 ==  
                                            sum(model[:Q_ji_k][t, l, k] for l in Omega_receiving[i], k in K) -
                                            sum(model[:Q_ij_k][t, l, k] for l in Omega_sending[i], k in K)
                    end)
    return 
end