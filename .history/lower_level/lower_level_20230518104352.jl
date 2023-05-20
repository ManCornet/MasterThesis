using JuMP, BilevelJuMP, Gurobi
using GraphRecipes, Graphs, Plots, Dates, DataFrames


model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "TimeLimit", 200)
set_optimizer_attribute(model, "Presolve", 0)


@objective(LOWER_MODEL, Min, sum(user_costs))


@variables(LOWER_MODEL, begin
    p_imp[N, T] >= 0, (container = Array) # active power imported at time t
    q_imp[N, T] >= 0, (container = Array)    # reactive power imported at time t
    p_exp[N, T] >= 0, (container = Array) # active power exported at time t
    q_exp[N, T] >= 0, (container = Array)   # reactive power exported at time t
    s_grid_max[N] >= 0, (container = Array) # grid imp/exp capacity
    p_pv[N, T] >= 0, (container = Array) # active power generated at time t
    #q_pv[N, T], (container = Array)   # reactive power generated at time t # DEACTIVATED in a first time to reimplement and test later.
    s_conv_pv[N] >= 0, (container = Array) # PV converter max capacity
    0 <= p_pv_max[N] <= MAX_PV_CAPACITY_PER_NODE, (container = Array) # PV capacity
    user_costs[N], (container = Array) # Total cost of each user, i.e. PV_costs + energy_costs - energy_revenues + grid_costs
    PV_costs[N] >= 0, (container = Array) # Costs related to PV investment for each user
    energy_costs[N] >= 0, (container = Array)  # Costs related to energy imported, for each user
    energy_revenues[N] >= 0, (container = Array)  # Costs related to energy imported, for each user
    grid_costs[N] >= 0, (container = Array) # Costs related to network for each user, both capacy and energy-based
end)

for i in Ns
    for t in T
        fix(p_imp[i, t], 0.0; force=true)
        fix(q_imp[i, t], 0.0; force=true)
        fix(p_exp[i, t], 0.0; force=true)
        fix(q_exp[i, t], 0.0; force=true)
        fix(s_grid_max[i], 0.0; force=true)
        fix(p_pv[i, t], 0.0; force=true)
        # fix(q_pv[i, t], 0.0; force=true) # DEACTIVATED in a first time to reimplement and test later.
    end
    fix(s_conv_pv[i], 0.0; force=true)
    fix(p_pv_max[i], 0.0; force=true)
    fix(user_costs[i], 0.0; force=true)
    fix(PV_costs[i], 0.0; force=true)
    fix(energy_costs[i], 0.0; force=true)
    fix(energy_revenues[i], 0.0; force=true)
    fix(grid_costs[i], 0.0; force=true)
end



@constraints(LOWER_MODEL, begin
    [i in Nu], user_costs[i] == energy_costs[i] - energy_revenues[i] + PV_costs[i] + grid_costs[i]
    [i in Nu], PV_costs[i] == BASE_POWER * (s_conv_pv[i] * CONVERTER_COST / AMORTIZATION_CONVERTER + p_pv_max[i] * PV_COST / AMORTIZATION_PV)
    [i in Nu], energy_costs[i] == BASE_POWER * sum(p_imp[i, t] * IMP_ELECTRICITY_ENRG_COST for t in T) * TIME_STEP * DAYS_A_YEAR / NB_PROFILES
    [i in Nu], energy_revenues[i] == BASE_POWER * sum(p_exp[i, t] * EXP_ELECTRICITY_ENRG_COST for t in T) * TIME_STEP * DAYS_A_YEAR / NB_PROFILES
    [i in Nu], grid_costs[i] == BASE_POWER * (s_grid_max[i] * GRID_CONNECTION_COST + sum((p_imp[i, t] * IMP_ELECTRICITY_DSO_COST + p_exp[i, t] * EXP_ELECTRICITY_DSO_COST) for t in T) * TIME_STEP * DAYS_A_YEAR / NB_PROFILES)
    [i in Nu, t in T], p_imp[i, t] - p_exp[i, t] == P_CONSUMPTION[i, t] - p_pv[i, t]
    [i in Nu, t in T], q_imp[i, t] - q_exp[i, t] == Q_CONSUMPTION[i, t] # - q_pv[i, t] DEACTIVATED in a first time to reimplement and test later.
    # [i in Nu, t in T], [s_grid_max[i], p_imp-p_exp[i,t], q_imp-q_exp[i,t]] in SecondOrderCone() # see next constraints
    [i in Nu, t in T], p_imp[i, t] <= s_grid_max[i]
    [i in Nu, t in T], q_imp[i, t] <= s_grid_max[i] # q_imp > 0
    [i in Nu, t in T], p_exp[i, t] <= s_grid_max[i]
    [i in Nu, t in T], q_exp[i, t] <= s_grid_max[i] # q_exp > 0

    # [i in Nu, t in T], q_pv[i, t] <= PV_MAX_Q * PV_PRODUCTION[i, t] * p_pv_max[i] # DEACTIVATED in a first time to reimplement and test later.
    # [i in Nu, t in T], -q_pv[i, t] <= PV_MAX_Q * PV_PRODUCTION[i, t] * p_pv_max[i] # DEACTIVATED in a first time to reimplement and test later.
    # [i in Nu, t in T], p_pv[i, t] <= PV_PRODUCTION[i, t] * p_pv_max[i] + PV_SLOPE * q_pv[i, t] # DEACTIVATED in a first time to reimplement and test later.
    # [i in Nu, t in T], p_pv[i, t] <= PV_PRODUCTION[i, t] * p_pv_max[i] - PV_SLOPE * q_pv[i, t] # DEACTIVATED in a first time to reimplement and test later.
    [i in Nu, t in T], p_pv[i, t] <= s_conv_pv[i] # The PV power output at time t is always bounded by the capacity of the converter
    [i in Nu, t in T], p_pv[i, t] <= PV_PRODUCTION[i, t] * p_pv_max[i] # The PV power output at time t is also bounded by the available sun power PV_PRODUCTION is in [W/wp]
    [i in Nu, t in T], p_exp[i, t] <= p_pv[i, t]  # NEW : because otherwise there are cases where p_exp > p_pv (Manon)
    #[i in Nu, t in T], q_exp[i, t] <= q_pv[i, t] # NEW : because otherwise there are cases where q_exp > q_pv (Manon)
end)