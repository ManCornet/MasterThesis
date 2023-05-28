# ---------------------------------------------------------------------------- #
#                               Bus variables                                  #
# ---------------------------------------------------------------------------- #
function _add_BusVariables!(model::JuMP.AbstractModel, ::NoDG)::Nothing 

    network_data = model[:network_data]
    T  = model[:time_steps]
    N  = get_nb_bus(network_data)
    Ns = get_nb_sub_bus(network_data)

    JuMP.@variables(model,   
                    begin 
                        V_sqr[1:T, 1:N]
                        P_sub[1:T, 1:Ns]                            
                        Q_sub[1:T, 1:Ns]                                 
                        S_sub[1:T, 1:Ns] >= 0
                        S_sub_capa[1:Ns] >= 0                                    
                        Beta[1:Ns], (binary=true)  
                    end
                    )
    return
end

function _add_BusVariables!(model::JuMP.AbstractModel, ::DG)::Nothing 

    network_data = model[:network_data]
    T  = model[:time_steps]
    N  = get_nb_bus(network_data)
    Ns = get_nb_sub_bus(network_data)
    Nu = get_nb_load_bus(network_data)

    JuMP.@variables(model,   
                    begin 
                        V_sqr[1:T, 1:N]
                        P_sub[1:T, 1:Ns]                            
                        Q_sub[1:T, 1:Ns]                                 
                        S_sub[1:T, 1:Ns] >= 0
                        S_sub_capa[1:Ns] >= 0                                    
                        p_pv[1:T, 1:Nu]  >= 0
                        s_conv_pv[1:Nu]  >= 0
                        p_pv_max[1:Nu]   >= 0
                        Beta[1:Ns], (binary=true)  
                    end
                    )

    # don't forget to add this is in operational constraints 
    # (2 config: strong constraints or relaxed version)
    
    # load_buses = get_load_buses(network_data)
    #JuMP.@constraint(model, [i = 1:Nu], p_pv_max[i] <= load_buses[i].max_pv_capa)
    return
end


# ---------------------------------------------------------------------------- #
#                             Branch variables                                 #
# ---------------------------------------------------------------------------- #
function _add_BranchVariables!(model::JuMP.AbstractModel, ::BIM)::Nothing 

    network_data = model[:network_data]
    T = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)
    N = get_nb_bus(network_data)

    JuMP.@variables(model,   
                    begin 
                        P_ij_k[1:T, 1:L, 1:K]
                        P_ji_k[1:T, 1:L, 1:K] 
                        Q_ij_k[1:T, 1:L, 1:K]                           
                        Q_ji_k[1:T, 1:L, 1:K]
                        X_ij_k_i[1:T, 1:L, 1:K, 1:N]
                        X_ij_k_re[1:T, 1:L, 1:K] >= 0
                        X_ij_k_im[1:T, 1:L, 1:K]
                        I_sqr_k[1:T, 1:L, 1:K]
                        Alpha[1:L, 1:K], (binary = true)
                        Y[1:L], (binary = true) 
                    end
                    )
    return
end

function _add_BranchVariables!(model::JuMP.AbstractModel, ::BFM)::Nothing 

    network_data = model[:network_data]
    T = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)

    JuMP.@variables(model,   
                    begin 
                        P_ij_k[1:T, 1:L, 1:K]
                        Q_ij_k[1:T, 1:L, 1:K]                           
                        I_sqr_k[1:T, 1:L, 1:K] >= 0
                    end
                    )
    return
end

# ---------------------------------------------------------------------------- #
#                          Conductor choice variables                          #
# ---------------------------------------------------------------------------- #
function _add_CondChoiceVariables!(model::JuMP.AbstractModel, ::OneConfig, ::Undirected)

    network_data = model[:network_data]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)

    JuMP.@variables(model,   
                    begin 
                        Alpha[1:L, 1:K], (binary = true)
                        Y[1:L], (binary = true) 
                    end
                    )
    return
end

function _add_CondChoiceVariables!(model::JuMP.AbstractModel, ::ReconfigAllowed, ::Undirected)

    network_data = model[:network_data]
    T = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)

    JuMP.@variables(model,   
                    begin 
                        Alpha[1:T, 1:L, 1:K], (binary = true)
                        Y[1:T, 1:L], (binary = true) 
                    end
                    )
    return
end

function _add_CondChoiceVariables!(model::JuMP.AbstractModel, ::OneConfig, ::Directed)

    network_data = model[:network_data]
    T = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)

    JuMP.@variables(model,   
                    begin 
                        Alpha[1:L, 1:K], (binary = true)
                        Y_send[1:L], (binary = true)  
                        Y_rec[1:L], (binary = true)  
                    end
                    )
    return
end

function _add_CondChoiceVariables!(model::JuMP.AbstractModel, ::ReconfigAllowed, ::Directed)

    network_data = model[:network_data]
    T = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)

    JuMP.@variables(model,   
                    begin 
                        Alpha[1:T, 1:L, 1:K], (binary = true)
                        Y_send[1:T, 1:L], (binary = true)  
                        Y_rec[1:T, 1:L], (binary = true)  
                    end
                    )
    return
end


# ---------------------------------------------------------------------------- #
#                          Radiality variables                                 #
# ---------------------------------------------------------------------------- #
# For now, only single commodity
# Here, we want to test:
# Reconfiguration vs only one configuration 
# 1° Single-Commodity constraints 
# 2° Multi-Commodity constraints 
# 3° SpanningTree constraints (called like that by Jabr what is the real name of that ?)
#    => parent child relation ship
# 4° + for each allow the version with reconfiguration etc
function _add_RadialityVariables!(model::JuMP.AbstractModel, ::OneConfig, ::SingleCommodityFlow)::Nothing 

    network_data = model[:network_data]
    T = model[:time_steps]
    L = get_nb_lines(network_data)

    JuMP.@variable(model, k_ij[1:L])

    return
end

function _add_RadialityVariables!(model::JuMP.AbstractModel, ::ReconfigAllowed, ::SingleCommodityFlow)::Nothing 

    network_data = model[:network_data]
    T = model[:time_steps]
    L = get_nb_lines(network_data)

    JuMP.@variable(model, k_ij[1:T, 1:L])

    return
end






