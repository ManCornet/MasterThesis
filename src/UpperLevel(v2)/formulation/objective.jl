function PV_coeff(tau, lambda)
    return (1 - 1/(1 + tau)^lambda)/tau  
end

# Capital recovery rate formula
# tau: interest rate 
# n : number of annuity received
function CRF(tau, n)
    return (tau * (1 + tau)^n)/((1 + tau)^n - 1)  
end


# function _set_Objective!(   model::JuMP.AbstractModel, 
#                             c::CurrentConstraintsFormulation,
#                             v::VoltageConstraintsFormulation)::Nothing

    
#     network_data = model[:network_data]
#     DSO_costs = model[:DSO_costs] 
#     WEIGHT_I = model[:DSO_costs].WEIGHT_I 
#     WEIGHT_V = model[:DSO_costs].WEIGHT_V

#     T = model[:time_steps]
#     L = network_data.nb_lines
#     K = network_data.nb_conductors
#     Ns = network_data.nb_substations

#     lines = network_data.lines
#     conductors = network_data.conductors
#     DAYS_A_YEAR = 365
#     pu_basis = network_data.pu_basis
   
#     # Here it depends on what ?
#     # if alpha depends on time or not then => Gamma
#     # how to compute line cost ?
#     # SUB_INSTALL COST verify it is in DSO_COSTS ?
#     # LOSS_COST ?
#     # TIME_STEP: we can get it but in minutes i think
#     # BASE_POWER is in 

#     # Here take into account the fact that substation can exist
    
#     Alpha = model[:Alpha]
#     S_sub_capa = model[:S_sub_capa]

#     MULTIPLIER = DAYS_A_YEAR / model[:nb_sign_days] 
#     fixed_costs = JuMP.@expression(
#                             model, 
#                             sum(Alpha[l, k] * lines[l].length 
#                                 * conductors[k].cost for l in 1:L, k in 1:K) +
#                             sum(S_sub_capa[i] * pu_basis.base_power * DSO_costs.substation for i in 1:Ns))
   
#     delta_t = model[:delta_t]
#     I_sqr_k = model[:I_sqr_k]
#     loss_costs = JuMP.@expression(
#                             model,
#                             MULTIPLIER * DSO_costs.loss * delta_t/60 *  
#                             pu_basis.base_power * sum(lines[l].length * conductors[k].r * 
#                             I_sqr_k[t, l, k] for t in 1:T, l in 1:L, k in 1:K)
#     )

#     violations = 0
#     if isa(c, RelaxedCurrents)
#         I_violation = model[:I_violation]
#         violations = JuMP.@expression(
#                                     model,
#                                     violations + WEIGHT_I * MULTIPLIER * 
#                                     sum(I_violation[t, l, k] 
#                                     for t in 1:T, l in 1:L, k in 1:K)
#         )
#     end

#     if isa(v, RelaxedVoltages)
#         V_violation_up = model[:V_violation_up]
#         V_violation_low = model[:V_violation_low]
#         violations = JuMP.@expression(
#                                     model,
#                                     violations + WEIGHT_V * MULTIPLIER * 
#                                     sum(V_violation_up[t, l, k] + V_violation_low[t, l, k] 
#                                     for t in 1:T, l in 1:L, k in 1:K)
#         )
#     end

#     JuMP.@objective(model, Min, fixed_costs / DSO_costs.amortization + loss_costs + violations)
#     return
# end


function _set_Objective!(   model::JuMP.AbstractModel, 
    c::CurrentConstraintsFormulation,
    v::VoltageConstraintsFormulation)::Nothing


    network_data = model[:network_data]
    DSO_costs = model[:DSO_costs] 
    WEIGHT_I = model[:DSO_costs].WEIGHT_I 
    WEIGHT_V = model[:DSO_costs].WEIGHT_V

    T = model[:time_steps]
    L = network_data.nb_lines
    K = network_data.nb_conductors
    Ns = network_data.nb_substations

    lines = network_data.lines
    conductors = network_data.conductors
    DAYS_A_YEAR = 365
    BASE_POWER = network_data.pu_basis.base_power

  
    JuMP.@variables(model, 
                    begin
                        DSO_fixed_costs 
                        DSO_loss_costs
                        DSO_op_limits
                    end)


    Alpha = model[:Alpha]
    S_sub_capa = model[:S_sub_capa]
    DELTA_T = model[:delta_t]
    I_sqr_k = model[:I_sqr_k]

    MULTIPLIER = DAYS_A_YEAR / model[:nb_sign_days] 
    violations = JuMP.AffExpr(0.0)
    
    if isa(c, RelaxedCurrents)
        I_violation = model[:I_violation]
        JuMP.add_to_expression!(violations, 
                                WEIGHT_I * MULTIPLIER,
                                sum(I_violation[t, l, k] 
                                for t in 1:T, l in 1:L, k in 1:K)
                            )
    end

    if isa(v, RelaxedVoltages)
        V_violation_up = model[:V_violation_up]
        V_violation_low = model[:V_violation_low]
        JuMP.add_to_expression!(violations, 
                                WEIGHT_V * MULTIPLIER,
                                sum(V_violation_up[t, l, k] + V_violation_low[t, l, k] 
                                for t in 1:T, l in 1:L, k in 1:K)
                            )
    end

    JuMP.@constraints(model,
                    begin
                        DSO_fixed_costs == (sum(Alpha[l, k] * lines[l].length 
                                            * conductors[k].cost for l in 1:L, k in 1:K) +
                                            sum(S_sub_capa[i] * BASE_POWER * DSO_costs.substation for i in 1:Ns))
                        DSO_loss_costs == MULTIPLIER * DSO_costs.loss * DELTA_T/60 *
                                         BASE_POWER * sum(lines[l].length * conductors[k].r *
                                         I_sqr_k[t, l, k] for t in 1:T, l in 1:L, k in 1:K)
                        DSO_op_limits == violations
                    end)

    # JuMP.add_to_expression!(objective, 
    #                             1.0/DSO_costs.amortization,
    #                             sum(Alpha[l, k] * lines[l].length 
    #                             * conductors[k].cost for l in 1:L, k in 1:K) 
    #                             + sum(S_sub_capa[i] * BASE_POWER * DSO_costs.substation for i in 1:Ns)
    #                         )

    # JuMP.add_to_expression!(objective, 
    #                         1.0,
    #                         MULTIPLIER * DSO_costs.loss * DELTA_T/60 
    #                         * BASE_POWER * sum(lines[l].length * conductors[k].r 
    #                         * I_sqr_k[t, l, k] for t in 1:T, l in 1:L, k in 1:K)
    #                     )
  
    JuMP.@objective(model, Min, DSO_fixed_costs/DSO_costs.amortization + DSO_loss_costs + DSO_op_limits)
    return
end





