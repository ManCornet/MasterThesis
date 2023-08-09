
function _set_Objective!(   model::JuMP.AbstractModel, 
                            c::CurrentConstraintsFormulation,
                            v::VoltageConstraintsFormulation)::Nothing

    # Definition of the upper and lower problems
    upper = model[:bilevel] ? Upper(model) : model
    lower = model[:bilevel] ? Lower(model) : model

    #############################################################
    #                       Upper objective                     #
    #############################################################

    # -- General network data --
    network_data = model[:network_data]
    T  = model[:time_steps]
    L  = get_nb_lines(network_data)
    K  = get_nb_conductors(network_data)
    Ns = get_nb_substations(network_data)
    Nu = get_nb_loads(network_data)

    # -- Operational limits violation term of the objective --
    DSO_costs    = model[:DSO_costs] 
    DAYS_A_YEAR  = 365
    MULTIPLIER   = DAYS_A_YEAR / model[:nb_sign_days] 

    violations = JuMP.@expression(upper, 0.0)
    if isa(c, RelaxedCurrents)
        I_violation = model[:I_violation]
        WEIGHT_I     = DSO_costs.WEIGHT_I 
        JuMP.add_to_expression!(violations, 
                                WEIGHT_I * MULTIPLIER,
                                sum(I_violation[t, l, k] 
                                for t in 1:T, l in 1:L, k in 1:K)
                            )
    end

    if isa(v, RelaxedVoltages)
        WEIGHT_V     = DSO_costs.WEIGHT_V
        V_violation_up  = model[:V_violation_up]
        V_violation_low = model[:V_violation_low]
        JuMP.add_to_expression!(violations, 
                                WEIGHT_V * MULTIPLIER,
                                sum(V_violation_up[t, l, k] + V_violation_low[t, l, k] 
                                for t in 1:T, l in 1:L, k in 1:K)
                            )
    end

    # -- Definition of the fixed investment costs and the loss costs --
    # Fetching the variables
    DSO_fixed_costs = model[:DSO_fixed_costs]
    DSO_loss_costs  = model[:DSO_loss_costs]
    DSO_op_limits   = model[:DSO_op_limits]
    Alpha           = model[:Alpha]
    S_sub_capa      = model[:S_sub_capa]
    I_sqr_k         = model[:I_sqr_k]
    grid_costs      = model[:grid_costs]

    # Fetching the data
    lines        = network_data.lines
    conductors   = network_data.conductors
    BASE_POWER   = network_data.pu_basis.base_power
    AMORTIZATION = DSO_costs.amortization
    IR           = DSO_costs.interest_rate
    LOSS_COST    = DSO_costs.loss
    SUB_COST     = DSO_costs.substation
    DELTA_T      = model[:delta_t]
   
    JuMP.@constraints(upper,
                    begin
                        DSO_fixed_costs == (sum(Alpha[l, k] * lines[l].length * conductors[k].cost for l in 1:L, k in 1:K) +
                                            sum(S_sub_capa[i] * BASE_POWER * SUB_COST for i in 1:Ns))
                        DSO_loss_costs == MULTIPLIER * LOSS_COST * DELTA_T/60 * BASE_POWER * sum(lines[l].length * conductors[k].r *
                        I_sqr_k[t, l, k] for t in 1:T, l in 1:L, k in 1:K)
                        DSO_op_limits == violations
                        DSO_fixed_costs * ((1 + IR)^AMORTIZATION) + DSO_loss_costs * AMORTIZATION <= sum(grid_costs[i] for i in 1:Nu) * AMORTIZATION
                    end)

    #############################################################
    #                       Lower objective                     #
    #############################################################
    user_costs          = model[:user_costs]
    energy_costs        = model[:energy_costs]
    energy_revenues     = model[:energy_revenues]
    PV_costs            = model[:PV_costs]
    grid_costs          = model[:grid_costs]
    s_conv_pv           = model[:s_conv_pv]
    p_pv_max            = model[:p_pv_max]
    p_imp               = model[:p_imp]
    p_exp               = model[:p_exp]
    s_grid_max          = model[:s_grid_max]

    User_costs          = model[:User_costs]
    CONV_COST           = User_costs.PV_conv
    PV_COST             = User_costs.PV
    AMORTIZATION_CONV   = User_costs.amortization_PVC
    AMORTIZATION_PV     = User_costs.amortization_PV
    EI_COST             = User_costs.EI
    EE_COST             = User_costs.EE
    DSOEI_COST          = User_costs.DSOEI
    DSOEE_COST          = User_costs.DSOEE
    GC_COST             = User_costs.GCC
   
    # # User costs expression
    user_costs_exp = JuMP.@expression(lower, [i=1:Nu], energy_costs[i] - energy_revenues[i] + PV_costs[i] + grid_costs[i])

    if model[:storage] 
        # Fetching required data
        storage_capacity     = model[:storage_capacity]
        storage_costs        = model[:storage_costs]
        STORAGE_COST         = User_costs.storage
        AMORTIZATION_STORAGE = User_costs.amortization_storage
        # Updating user costs with storage costs
        user_costs_exp = JuMP.@expression(lower, [i=1:Nu], user_costs_exp[i] + storage_costs[i])

        JuMP.@constraints(lower, begin
            [i=1:Nu], storage_costs[i] == BASE_POWER * (storage_capacity[i] * STORAGE_COST / AMORTIZATION_STORAGE)
        end)
    end

    JuMP.@constraints(lower, 
                begin
                    [i in 1:Nu], user_costs[i] == user_costs_exp[i]
                    [i in 1:Nu], PV_costs[i] == BASE_POWER * (s_conv_pv[i] * CONV_COST / AMORTIZATION_CONV + p_pv_max[i] * PV_COST / AMORTIZATION_PV)
                    [i in 1:Nu], energy_costs[i] == MULTIPLIER * BASE_POWER * DELTA_T/60 * sum(p_imp[t, i] * EI_COST for t in 1:T)
                    [i in 1:Nu], energy_revenues[i] == MULTIPLIER * BASE_POWER * DELTA_T/60 * sum(p_exp[t, i] * EE_COST for t in 1:T)  
                    [i in 1:Nu], grid_costs[i] == s_grid_max[i] * GC_COST + MULTIPLIER * BASE_POWER * DELTA_T/60 * sum((p_imp[t, i] * DSOEI_COST + p_exp[t, i] * DSOEE_COST) for t in 1:T)
                end)

    if model[:bilevel] 
        JuMP.@objective(upper, Min, DSO_fixed_costs/AMORTIZATION + DSO_loss_costs + DSO_op_limits)
        JuMP.@objective(lower, Min, sum(user_costs))
    else 
        WEIGHT_OBJ1 = DSO_costs.WEIGHT_OBJ
        WEIGHT_OBJ2 = User_costs.WEIGHT_OBJ
        JuMP.@objective(model, Min, WEIGHT_OBJ1 * (DSO_fixed_costs/AMORTIZATION + DSO_loss_costs + DSO_op_limits) + WEIGHT_OBJ2 * (sum(user_costs)))
    end

    return
end