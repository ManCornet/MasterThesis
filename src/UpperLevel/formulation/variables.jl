# ---------------------------------------------------------------------------- #
#                               Bus variables                                  #
# ---------------------------------------------------------------------------- #
function _add_BusVariables!(model::JuMP.AbstractModel, ::NoDG)::Nothing 

    network_data = model[:network_data]
    T  = model[:time_steps]
    N  = get_nb_buses(network_data)
    Ns = get_nb_substations(network_data)

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
    N  = get_nb_buses(network_data)
    Ns = get_nb_substations(network_data)
    Nu = get_nb_loads(network_data)

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
    N = get_nb_buses(network_data)

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
Vararg_Tuple{T} = Tuple{Vararg{T}}
compute_index(index::Vararg_Tuple{T}, t::T, ::OneConfig) where {T} = index
compute_index(index::Vararg_Tuple{T}, t::T, ::ReconfigAllowed) where {T} = (t, index...)

function _add_CondChoiceVariables!(model::JuMP.AbstractModel, topology_choice::TopologyChoiceFormulation, ::Undirected)

    T = model[:time_steps]
    network_data = model[:network_data]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)
    
    JuMP.@variables(model,   
                    begin 
                        Alpha[compute_index((1:L, 1:K), 1:T, topology_choice)...], (binary = true)
                        Y[compute_index((1:L), 1:T, topology_choice)...], (binary = true) 
                    end
                    )
    return
end

function _add_CondChoiceVariables!(model::JuMP.AbstractModel, topology_choice::TopologyChoiceFormulation, ::Directed)

    network_data = model[:network_data]
    T = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)

    JuMP.@variables(model,   
                    begin 
                        Alpha[compute_index((1:L, 1:K), 1:T, topology_choice)...], (binary = true)
                        Y_send[compute_index((1:L), 1:T, topology_choice)...], (binary = true)  
                        Y_rec[compute_index((1:L), 1:T, topology_choice)...], (binary = true)  
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
function _add_RadialityVariables!(model::JuMP.AbstractModel, topology_choice::TopologyChoiceFormulation, ::SingleCommodityFlow)::Nothing 

    network_data = model[:network_data]
    T = model[:time_steps]
    L = get_nb_lines(network_data)

    JuMP.@variable(model, k_ij[compute_index((1:L), 1:T, topology_choice)...])

    return
end





