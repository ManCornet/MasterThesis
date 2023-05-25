

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
                        beta[1:Ns], (binary=true)  
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
                        p_pv[1:T, 1:Nu]  >= 0
                        q_pv[1:T, 1:Nu]  >= 0
                        s_conv_pv[1:Nu]  >= 0
                        p_pv_max[1:Nu]   >= 0
                        beta[1:Ns], (binary=true)  
                    end
                    )

    load_buses = get_load_buses(network_data)

    JuMP.@constraint(model, [i = 1:Nu], p_pv_max[i] <= load_buses[i].capa_pv_max)
    return
end


# ---------------------------------------------------------------------------- #
#                             Branch variables                                 #
# ---------------------------------------------------------------------------- #