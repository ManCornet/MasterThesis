function PV_coeff(tau, lambda)
    return (1 - 1/(1 + tau)^lambda)/tau  
end

# Capital recovery rate formula
# tau: interest rate 
# n : number of annuity received
function CRF(tau, n)
    return (tau * (1 + tau)^n)/((1 + tau)^n - 1)  
end


function _set_Objective!(   model::JuMP.AbstractModel, 
                            c::CurrentConstraintsFormulation,
                            v::VoltageConstraintsFormulation;
                            WEIGHT_I=nothing,
                            WEIGHT_V=nothing)::Nothing

    network_data = model[:network_data]
    DSO_costs = model[:DSO_costs] 

    T = model[:time_steps]
    L = network_data.nb_lines
    K = network_data.nb_conductors
    Ns = network_data.nb_substations

    lines = network_data.lines
    conductors = network_data.conductors
    DAYS_A_YEAR = 365
    pu_basis = network_data.pu_basis
   
    # Here it depends on what ?
    # if alpha depends on time or not then => Gamma
    # how to compute line cost ?
    # SUB_INSTALL COST verify it is in DSO_COSTS ?
    # LOSS_COST ?
    # TIME_STEP: we can get it but in minutes i think
    # BASE_POWER is in 
    MULTIPLIER = DAYS_A_YEAR / model[:nb_sign_days] 
    fixed_costs = JuMP.@expression(
                            model, 
                            sum(model[:Alpha][l, k] * lines[l].length 
                                * conductors[k].cost for l in 1:L, k in 1:K) +
                            sum(model[:S_sub_capa][i] * pu_basis.base_power * DSO_costs.substation for i in 1:Ns))
   

    loss_costs = JuMP.@expression(
                            model,
                            MULTIPLIER * DSO_costs.loss * model[:delta_t]/60 *  
                            pu_basis.base_power * sum(lines[l].length * conductors[k].r * 
                            model[:I_sqr_k][t, l, k] for t in 1:T, l in 1:L, k in 1:K)
    )

    violations = 0
    if isa(c, RelaxedCurrents)
        violations = JuMP.@expression(
                                    model,
                                    violations + WEIGHT_I * MULTIPLIER * 
                                    sum(model[:I_violation][t, l, k] 
                                    for t in 1:T, l in 1:L, k in 1:K)
        )
    end

    if isa(v, RelaxedVoltages)
        violations = JuMP.@expression(
                                    model,
                                    violations + WEIGHT_V * MULTIPLIER * 
                                    sum(model[:V_violation_up][t, l, k] + model[:V_violation_low][t, l, k] 
                                    for t in 1:T, l in 1:L, k in 1:K)
        )
    end

    JuMP.@objective(model, Min, fixed_costs / DSO_costs.amortization + loss_costs + violations)
    return
end