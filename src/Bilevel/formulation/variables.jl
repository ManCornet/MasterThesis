# ---------------------------------------------------------------------------- #
#                               Bus variables                                  #
# ---------------------------------------------------------------------------- #
function _add_BusVariables!(model::JuMP.AbstractModel)::Nothing 

    network_data = model[:network_data]
    T  = model[:time_steps]
    N  = get_nb_buses(network_data)
    Ns = get_nb_substations(network_data)
    Nu = get_nb_loads(network_data)
    buses = network_data.buses

    upper = model[:bilevel] ? Upper(model) : model
    lower = model[:bilevel] ? Lower(model) : model

    JuMP.@variables(upper,   
                    begin 
                        V_sqr[1:T, 1:N] >= 0
                        P_sub[1:T, 1:Ns]                            
                        Q_sub[1:T, 1:Ns]                                 
                        S_sub[1:T, 1:Ns] >= 0
                        S_sub_capa[1:Ns] >= 0                                    
                        Beta[1:Ns], Bin
                        DSO_fixed_costs >= 0
                        DSO_loss_costs >= 0
                        DSO_op_limits >= 0
                    end
                    )

    JuMP.@variables(lower,   
        begin 
            p_imp[1:T, 1:Nu] >= 0
            q_imp[1:T, 1:Nu] >= 0
            p_exp[1:T, 1:Nu] >= 0
            q_exp[1:T, 1:Nu] >= 0
            s_grid_max[1:Nu] >= 0
            # pv limit 
            p_pv[1:T, 1:Nu]  >= 0
            #q_pv[1:T, 1:Nu]  >= 0  # DEACTIVATED in a first time to reimplement and test later.
            s_conv_pv[1:Nu]  >= 0
            p_pv_max[1:Nu]   >= 0
            user_costs[1:Nu]              # Cost per user: PV_costs + energy_costs - energy_revenues + grid_costs
            PV_costs[1:Nu] >= 0           # Costs for PV investment
            energy_costs[1:Nu] >= 0       # Costs related to energy imported, for each user
            energy_revenues[1:Nu] >= 0    # Costs related to energy imported, for each user
            grid_costs[1:Nu] >= 0         # Costs related to network for each user, both capacy and energy-based      
        end
        )

    for i in 1:Nu
        JuMP.set_upper_bound(p_pv_max[i], buses[Ns + i].max_pv_capa)
    end
    
    if model[:storage] 
        JuMP.@variables(lower,   
            begin 
            p_storage_charge[1:T, 1:Nu] >= 0               # active power to storage device at time t, positive when the battery is charging
            p_storage_discharge[1:T, 1:Nu] >= 0             # active power to storage device at time t, positive when the battery is charging
            storage_state[1:T, 1:Nu] >= 0      # storage capacity
            storage_capacity[1:Nu] >= 0        # storage capacity
            storage_costs[1:Nu] >= 0           # Costs related to storage investment for each user
            end)
    end


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
                        X_ij_k_i[1:T, 1:L, 1:K, 1:N] >= 0
                        X_ij_k_re[1:T, 1:L, 1:K] >= 0
                        X_ij_k_im[1:T, 1:L, 1:K]
                        I_sqr_k[1:T, 1:L, 1:K] 
                    end
                    )

    for i in 1:N, k in 1:K, l in 1:L, t in 1:T
        ifrom = network_data.lines[l].edge.from_node.id
        ito   = network_data.lines[l].edge.to_node.id
        #println("($ifrom, $ito)")
        if i in (ifrom, ito) 
            continue
        else
            JuMP.fix(X_ij_k_i[t, l, k, i], 0.0; force=true)
        end
    end
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
                        P_ij[1:T, 1:L]
                        Q_ij_k[1:T, 1:L, 1:K]          
                        Q_ij[1:T, 1:L]
                        I_sqr_k[1:T, 1:L, 1:K] >= 0
                        I_sqr[1:T, 1:L]
                    end
                    )
    JuMP.@constraints(model, begin 
                    [t=1:T, l=1:L], P_ij[t,l] == sum(P_ij_k[t,l,k] for k in 1:K)
                    [t=1:T, l=1:L], Q_ij[t,l] == sum(Q_ij_k[t,l,k] for k in 1:K)
                    [t=1:T, l=1:L], I_sqr[t,l] == sum(I_sqr_k[t,l,k] for k in 1:K)
                end)
    return
end

# ---------------------------------------------------------------------------- #
#                          Conductor choice variables                          #
# ---------------------------------------------------------------------------- #


function _add_CondChoiceVariables!(model::JuMP.AbstractModel, topology_choice::TopologyChoiceFormulation, ::Undirected)

    T = model[:time_steps]
    network_data = model[:network_data]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)
    Nu = get_nb_loads(network_data)
    
    
    if isa(topology_choice, OneConfig) 
        JuMP.@variables(model,   
                    begin 
                        Alpha[1:L, 1:K], Bin
                        Y[1:L], Bin
                    end
                    )

        JuMP.@constraints(  model, 
                            begin
                                [l=1:L], sum(Alpha[l, k] for k in 1:K) == Y[l]
                                sum(Y[l] for l in 1:L) == Nu 
                            end
                        )
                        
    elseif isa(topology_choice, ReconfigAllowed)
        println("hello")
        JuMP.@variables(model,   
                    begin 
                        Alpha[1:L, 1:K], Bin
                        Gamma[1:T, 1:L, 1:K], Bin
                        Y[1:T, 1:L], Bin 
                    end
                    )

        JuMP.@constraints(  model, 
                            begin
                                [l=1:L], sum(Alpha[l, k] for k in 1:K) <= 1 # we select only one conductor
                                [t=1:T, l=1:L], sum(Gamma[t, l, k] for k in 1:K) == Y[t, l]
                                [t=1:T, l=1:L, k=1:K], Gamma[t, l, k] <= Alpha[l, k]
                                [t=1:T], sum(Y[t, l] for l in 1:L) == Nu # radial in operation
                            end
                        )
    end
    return
end

function _add_CondChoiceVariables!(model::JuMP.AbstractModel, topology_choice::TopologyChoiceFormulation, ::Directed)::Nothing
    
    network_data = model[:network_data]
    T = model[:time_steps]
    L = get_nb_lines(network_data)
    K = get_nb_conductors(network_data)
    Nu = get_nb_loads(network_data)       

    # If one line has been built (always select the same conductor for all time_steps)
    # Add the constraints linking alpha and other constraints

    # Add the links between
    if isa(topology_choice, OneConfig) 
        JuMP.@variables(model,   
                    begin 
                        Alpha[1:L, 1:K], Bin
                        Y_send[1:L], Bin 
                        Y_rec[1:L], Bin 
                    end
                    )

        JuMP.@constraints(  model, 
                            begin
                                [l=1:L], sum(Alpha[l, k] for k in 1:K) == Y_send[l] + Y_rec[l] 
                                sum(Y_send[l] + Y_rec[l] for l in 1:L) == Nu 
                            end
                        )

    elseif isa(topology_choice, ReconfigAllowed)
        JuMP.@variables(model,   
                    begin 
                        Alpha[1:L, 1:K], (binary = true)
                        Gamma[1:T, 1:L, 1:K], (binary = true)
                        Y_send[1:T, 1:L], (binary = true)  
                        Y_rec[1:T, 1:L], (binary = true)  
                    end
                    )
        JuMP.@constraints(  model, 
                            begin
                                [l=1:L], sum(Alpha[l, k] for k in 1:K) <= 1 # we select only one conductor
                                [t=1:T, l=1:L], sum(Gamma[t, l, k] for k in 1:K) == Y_send[t, l] + Y_rec[t, l]
                                [t=1:T, l=1:L, k=1:K], Gamma[t, l, k] <= Alpha[l, k]
                                [t=1:T], sum(Y_send[t, l] + Y_rec[t, l] for l in 1:L) == Nu # radial in operation
                            end
                        )
    end
   
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

    if isa(topology_choice, OneConfig)
        JuMP.@variable(model, k_ij[1:L])
    elseif isa(topology_choice, ReconfigAllowed)
        JuMP.@variable(model, k_ij[1:T, 1:L])
    end

    return
end

function _add_RadialityVariables!(model::JuMP.AbstractModel,topology_choice::TopologyChoiceFormulation, ::MultiCommodityFlow)::Nothing 

    network_data = model[:network_data]
    T = model[:time_steps]
    L = get_nb_lines(network_data)
    N = get_nb_loads(network_data) + get_nb_substations(network_data)
    if isa(topology_choice, OneConfig)
        JuMP.@variable(model, k_ij[1:L, 1:N])
    elseif isa(topology_choice, ReconfigAllowed)
        JuMP.@variable(model, k_ij[1:T, 1:L, 1:N])
    end

    return
end







