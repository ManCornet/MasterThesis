using JuMP, BilevelJuMP, Gurobi
using GraphRecipes, Graphs, Plots, Dates, DataFrames
if (!@isdefined LAUNCH_SENSITIVITY_ANALYSIS) || (LAUNCH_SENSITIVITY_ANALYSIS != true)
    include("parameters_BFM_1P.jl")
end
include("export_xlsx.jl")


SINGLE_LEVEL = true

if SINGLE_LEVEL
    model = Model(Gurobi.Optimizer)
    UPPER_MODEL = model
    LOWER_MODEL = model
else
    #model = BilevelModel(Gurobi.Optimizer, mode=BilevelJuMP.SOS1Mode())
    model = BilevelModel(Gurobi.Optimizer, mode=BilevelJuMP.StrongDualityMode())
    UPPER_MODEL = Upper(model)
    LOWER_MODEL = Lower(model)
end

# set_silent(model)
set_optimizer_attribute(model, "TimeLimit", 600)
set_optimizer_attribute(model, "MIPGap", 1e-2)
set_optimizer_attribute(model, "MIPFocus", 1)
set_optimizer_attribute(model, "NonConvex", 2) # set this to test conic constraints equality
#set_optimizer_attribute(model, "DualReductions", 0)

@variables(UPPER_MODEL, begin
    MIN_VOLTAGE^2 <= V_sqr[N, T] <= MAX_VOLTAGE^2, (container = Array)
    I_sqr_k[L, K, T] >= 0, (container = Array)
    #I_sqr[L, T] >= 0, (container = Array)
    P_sub[Ns, T], (container = Array)
    Q_sub[Ns, T], (container = Array)
    S_sub[Ns, T] >= 0, (container = Array) # Because apparent power
    P_ij_k[L, K, T], (container = Array)
    Q_ij_k[L, K, T], (container = Array)
    #P_ij[L, T], (container = Array)
    #Q_ij[L, T], (container = Array)
    S_sub_capa[Ns] >= 0, (container = Array)
    Y[L, T], (container=Array, binary=true)
    Alpha[L, K], (container=Array, binary=true) # rechange here: remove the T
    Beta[Ns], (container = Array, binary=true)
    Curr_limit[L, K, T], (container=Array, binary=true)
    Slack[L, K, T], (container = Array)
    DSO_fixed_costs
    DSO_loss_costs
    DSO_op_limits
end)
if ACTIVATE_RADIALITY_CONSTRAINTS
    # SINGLE COMMODITY FLOW VARIABLE
    @variable(UPPER_MODEL, k_ij[L, T], container = Array)
end
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

@objective(UPPER_MODEL, Min, DSO_fixed_costs / AMORTIZATION_DSO + DSO_loss_costs + DSO_op_limits + (SINGLE_LEVEL ? sum(user_costs) : 0))
#@objective(UPPER_MODEL, Min, DSO_fixed_costs + DSO_loss_costs * AMORTIZATION_DSO + (SINGLE_LEVEL ? sum(user_costs) : 0))
if !(SINGLE_LEVEL)
    @objective(LOWER_MODEL, Min, sum(user_costs))
end


@constraints(UPPER_MODEL, begin
    # DSO costs PER YEAR
    DSO_fixed_costs == (sum(Alpha[l, k] * line_cost[l, k] for l in L, k in K) +
                        sum(S_sub_capa[i] * SUB_INSTALL_COST for i in Ns)) # TOTAL
    DSO_loss_costs == DAYS_A_YEAR * LOSS_COST * TIME_STEP / NB_PROFILES *
                      BASE_POWER * sum(R[l, k] * I_sqr_k[l, k, t] for l in L, k in K, t in T) # PER YEAR
    DSO_fixed_costs * (1 + DSO_INTEREST_RATE)^AMORTIZATION_DSO + DSO_loss_costs * AMORTIZATION_DSO <= sum(grid_costs[i] for i in Nu) * AMORTIZATION_DSO

    # Operational limits violation term
    DSO_op_limits == WEIGHT_OP_LIMITS * sum(Curr_limit[l, k, t] for l in L, k in K, t in T) * DAYS_A_YEAR / NB_PROFILES

    # Power balances
    [i = Ns, t = T], -P_sub[i, t] == sum(P_ij_k[l, k, t] - R[l, k] * I_sqr_k[l, k, t] for l in Omega_receiving[i], k in K) -
                                     sum(P_ij_k[l, k, t] for l in Omega_sending[i], k in K)
    [i = Ns, t = T], -Q_sub[i, t] == sum(Q_ij_k[l, k, t] - X[l, k] * I_sqr_k[l, k, t] for l in Omega_receiving[i], k in K) -
                                     sum(Q_ij_k[l, k, t] for l in Omega_sending[i], k in K)
    [i = Nu, t = T], p_imp[i, t] - p_exp[i, t] == sum(P_ij_k[l, k, t] - R[l, k] * I_sqr_k[l, k, t] for l in Omega_receiving[i], k in K) -
                                                  sum(P_ij_k[l, k, t] for l in Omega_sending[i], k in K)
    [i = Nu, t = T], q_imp[i, t] - q_exp[i, t] == sum(Q_ij_k[l, k, t] - X[l, k] * I_sqr_k[l, k, t] for l in Omega_receiving[i], k in K) -
                                                  sum(Q_ij_k[l, k, t] for l in Omega_sending[i], k in K)

    # Distflow
    [l = L, t = T], V_sqr[line_ends[l][2], t] - V_sqr[line_ends[l][1], t] <=
                    -sum(2 * (R[l, k] * P_ij_k[l, k, t] + X[l, k] * Q_ij_k[l, k, t]) +
                         (R[l, k]^2 + X[l, k]^2) * I_sqr_k[l, k, t] for k in K) +
                    M * (1 - Y[l])
    [l = L, t = T], V_sqr[line_ends[l][2], t] - V_sqr[line_ends[l][1], t] >=
                    -sum(2 * (R[l, k] * P_ij_k[l, k, t] + X[l, k] * Q_ij_k[l, k, t]) +
                         (R[l, k]^2 + X[l, k]^2) * I_sqr_k[l, k, t] for k in K) -
                    M * (1 - Y[l])

    #[l = L, t = T], I_sqr[l, t] == sum(I_sqr_k[l, k, t] for k in K)
    #[l = L, t = T], P_ij[l, t] == sum(P_ij_k[l, k, t] for k in K)
    #[l = L, t = T], Q_ij[l, t] == sum(Q_ij_k[l, k, t] for k in K)

    #[l = L, t = T], [V_sqr[line_ends[l][1], t] / 2; sum(I_sqr_k[l, k, t] for k in K); sum(P_ij_k[l, k, t] for k in K); sum(Q_ij_k[l, k, t] for k in K)] in RotatedSecondOrderCone()
    [l = L, t = T], V_sqr[line_ends[l][1], t] * sum(I_sqr_k[l, k, t] for k in K) == sum(P_ij_k[l, k, t] for k in K)^2 + sum(Q_ij_k[l, k, t] for k in K)^2
    #[l = L, k = K, t = T], [V_sqr[line_ends[l][1], t] / 2; I_sqr_k[l, k, t]; P_ij_k[l, k, t]; Q_ij_k[l, k, t]] in RotatedSecondOrderCone()

    #[i = Nu], sum(Y[l] for l in Omega_sending[i]) + sum(Y[l] for l in Omega_receiving[i]) >= 1

    # Operational limits
    #lorsque slack est négative => intger variable must be one 
    #si slack positive => binary variable == 0
    #si slack négatife => binary variable == 1 (idée forcer la variable à valoir 1 si nega)
    #-slack <= M * bianry variable 
    [l = L, k = K, t = T], Slack[l, k, t] <= MAX_CURRENT[l, k]^2 * Curr_limit[l, k, t]
    [l = L, k = K, t = T], I_sqr_k[l, k, t] - Slack[l, k, t] <= MAX_CURRENT[l, k]^2 * Alpha[l, k] # indispensable mais relaxée
    #[l = L, k = K, t = T], I_sqr_k[l, k, t] <= MAX_CURRENT[l, k] * Alpha[l, k] # indispensable
    [l = L, k = K, t = T], P_ij_k[l, k, t] <= MAX_CURRENT[l, k] * MAX_VOLTAGE * Alpha[l, k] # indispensable
    [l = L, k = K, t = T], P_ij_k[l, k, t] >= -MAX_CURRENT[l, k] * MAX_VOLTAGE * Alpha[l, k] # indispensable
    [l = L, k = K, t = T], Q_ij_k[l, k, t] <= MAX_CURRENT[l, k] * MAX_VOLTAGE * Alpha[l, k] # indispensable
    [l = L, k = K, t = T], Q_ij_k[l, k, t] >= -MAX_CURRENT[l, k] * MAX_VOLTAGE * Alpha[l, k] # indispensable

    # Substation apparent power limits
    [i = Ns, t = T], [S_sub[i, t]; P_sub[i, t]; Q_sub[i, t]] in SecondOrderCone()
    #[i = Ns, t = T], S_sub[i, t]^2 ==  P_sub[i, t]^2 + Q_sub[i, t]^2
    [i = Ns, t = T], S_sub[i, t] <= S_rating_init[i] + S_sub_capa[i]
    [i = Ns], S_sub_capa[i] <= Beta[i] * S_rating_max[i]

    # Voltage references
    [i in Ns_init, t in T], V_sqr[i, t] == 1
    [i in Ns_notinit, t in T], V_sqr[i, t] - 1 <= (MAX_VOLTAGE^2 - 1) * (1 - Beta[i])
    [i in Ns_notinit, t in T], V_sqr[i, t] - 1 >= (MIN_VOLTAGE^2 - 1) * (1 - Beta[i])

    # Load over-satisfaction
    #[t in T], sum(P_sub[i, t] for i in Ns) + sum(p_exp[i, t] - p_imp[i, t] for i in Nu) >= 0
    #[t in T], sum(Q_sub[i, t] for i in Ns) + sum(q_exp[i, t] - q_imp[i, t] for i in Nu) >= 0

    # Radiality constraint
    [l = L, t = T], sum(Alpha[l, k] for k in K) == Y[l, t]
    [t in T], sum(Y[l, t] for l in L) == Nu_size

    # Bound on s_grid_max
    #[i=Nu], s_grid_max[i] <= sum(MAX_CURRENT[l, k] * MAX_VOLTAGE * Alpha[l, k] for l in Omega_receiving[i], k in K) + 
    #                        sum(MAX_CURRENT[l, k] * MAX_VOLTAGE * Alpha[l, k] for l in Omega_sending[i], k in K) # NEW : doesn't change anything to the solution except faster !
end)
if ACTIVATE_RADIALITY_CONSTRAINTS
    @constraint(UPPER_MODEL, single_comm_constr1[i=Nu, t=T],
        -sum(k_ij[l, t] for l in Omega_receiving[i])
        +
        sum(k_ij[l, t] for l in Omega_sending[i]) == -1)

    @constraint(UPPER_MODEL, single_comm_constr2[i=Ns, t=T],
        -sum(k_ij[l, t] for l in Omega_receiving[i])
        +
        sum(k_ij[l, t] for l in Omega_sending[i]) >= 0)

    @constraint(UPPER_MODEL, single_comm_constr3[i=Ns_init,t=T],
        -sum(k_ij[l, t] for l in Omega_receiving[i])
        +
        sum(k_ij[l, t] for l in Omega_sending[i]) <= Nu_size)


    @constraint(UPPER_MODEL, single_comm_constr5[i=Ns_notinit,t=T],
        -sum(k_ij[l,t] for l in Omega_receiving[i])
        +
        sum(k_ij[l,t] for l in Omega_sending[i]) <= Nu_size * Beta[i])

    @constraint(UPPER_MODEL,
        single_comm_constr6[l=L,t=T],
        k_ij[l,t] <= Nu_size * Y[l,t]
    )

    @constraint(UPPER_MODEL,
        single_comm_constr7[l=L,t=T],
        k_ij[l,t] >= -Nu_size * Y[l,t]
    )

end

# if MAX_CO2_BUDGET > 0
#     @constraint(UPPER_MODEL), CO2_budget <= MAX_CO2_BUDGET)
# end
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

# m2=BilevelJuMP._build_single_model(model)

optimize!(model)

include("results_tables.jl")
include("plot_results.jl")