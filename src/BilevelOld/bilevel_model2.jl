#-----------------------------------------------------------------------------
#
#                           - TFE : Bilevel DNEP - 
#                             University of Liege
#
#-----------------------------------------------------------------------------
# Created By  : Manon Cornet
# Created Date: Wednesday May 24 2023
#
# bilevel_model:
#   file containing the bilevel model
#
# =============================================================================
#                                   Imports
# =============================================================================
# Activating the julia environement
# Path: I must add in terminal julia --project -- src/main.jl --

using JuMP, BilevelJuMP, Gurobi

# =============================================================================
#                                   Functions
# =============================================================================

function bilevel!(network::Network, 
                  DSO_costs::DSOCosts,
                  User_costs::UserCosts; 
                  BILEVEL::Bool=false, RADIALITY::Bool=true)
    

    
    if BILEVEL
        #model = BilevelModel(Gurobi.Optimizer, mode=BilevelJuMP.SOS1Mode())
        model = BilevelModel(Gurobi.Optimizer, mode=BilevelJuMP.StrongDualityMode())
        UPPER_MODEL = Upper(model)
        LOWER_MODEL = Lower(model)
    else
        model = Model(Gurobi.Optimizer)
        UPPER_MODEL = model
        LOWER_MODEL = model
    end

    # set_silent(model)
    set_optimizer_attribute(model, "TimeLimit", 600)
    set_optimizer_attribute(model, "MIPGap", 1e-2)
    set_optimizer_attribute(model, "MIPFocus", 1)
    #set_optimizer_attribute(model, "DualReductions", 0)

    @variables(UPPER_MODEL, begin
        MIN_VOLTAGE^2 <= V_sqr[T, N] <= MAX_VOLTAGE^2, (container = Array)
        I_sqr_k[T, L, K] >= 0, (container = Array)
        I_sqr[T, L] >= 0, (container = Array)
        P_sub[T, Ns], (container = Array)
        Q_sub[T, Ns], (container = Array)
        S_sub[T, Ns] >= 0, (container = Array) # Because apparent power
        P_ij_k[T, L, K], (container = Array)
        Q_ij_k[T, L, K], (container = Array)
        P_ij[T, L], (container = Array)
        Q_ij[T, L], (container = Array)
        S_sub_capa[Ns] >= 0, (container = Array)
        Y[L], (container=Array, binary=true)
        Alpha[L, K], (container=Array, binary=true)
        Beta[Ns], (container = Array, binary=true)
        Curr_limit[T, L, K], (container=Array, binary=true)
        Slack[T, L, K], (container = Array)
        DSO_fixed_costs
        DSO_loss_costs
        DSO_op_limits
    end)
    if ACTIVATE_RADIALITY_CONSTRAINTS
        # SINGLE COMMODITY FLOW VARIABLE
        @variable(UPPER_MODEL, k_ij[L], container = Array)
    end
    @variables(LOWER_MODEL, begin
        p_imp[T, N] >= 0, (container = Array) # active power imported at time t
        q_imp[T, N] >= 0, (container = Array)    # reactive power imported at time t
        p_exp[T, N] >= 0, (container = Array) # active power exported at time t
        q_exp[T, N] >= 0, (container = Array)   # reactive power exported at time t
        s_grid_max[N] >= 0, (container = Array) # grid imp/exp capacity
        p_pv[T, N] >= 0, (container = Array) # active power generated at time t
        #q_pv[T, N], (container = Array)   # reactive power generated at time t # DEACTIVATED in a first time to reimplement and test later.
        s_conv_pv[N] >= 0, (container = Array) # PV converter max capacity
        0 <= p_pv_max[N] <= MAX_PV_CAPACITY_PER_NODE, (container = Array) # PV capacity
        user_costs[N], (container = Array) # Total cost of each user, i.e. PV_costs + energy_costs - energy_revenues + grid_costs
        PV_costs[N] >= 0, (container = Array) # Costs related to PV investment for each user
        energy_costs[N] >= 0, (container = Array)  # Costs related to energy imported, for each user
        energy_revenues[N] >= 0, (container = Array)  # Costs related to energy imported, for each user
        grid_costs[N] >= 0, (container = Array) # Costs related to network for each user, both capacy and energy-based
    end)

    for t in T
        for i in Ns
            fix(p_imp[t, i], 0.0; force=true)
            fix(q_imp[t, i], 0.0; force=true)
            fix(p_exp[t, i], 0.0; force=true)
            fix(q_exp[t, i], 0.0; force=true)
            fix(s_grid_max[i], 0.0; force=true)
            fix(p_pv[t, i], 0.0; force=true)
            # fix(q_pv[t, i], 0.0; force=true) # DEACTIVATED in a first time to reimplement and test later.
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
                        BASE_POWER * sum(R[l, k] * I_sqr_k[t, l, k] for t in T, l in L, k in K) # PER YEAR
        DSO_fixed_costs * (1 + DSO_INTEREST_RATE)^AMORTIZATION_DSO + DSO_loss_costs * AMORTIZATION_DSO <= sum(grid_costs[i] for i in Nu) * AMORTIZATION_DSO

        # Operational limits violation term
        DSO_op_limits == WEIGHT_OP_LIMITS * sum(Curr_limit[t, l, k] for t in T, l in L, k in K) * DAYS_A_YEAR / NB_PROFILES

        # Power balances
        [t = T, i = Ns], -P_sub[t, i] == sum(P_ij_k[t, l, k] - R[l, k] * I_sqr_k[t, l, k] for l in Omega_receiving[i], k in K) -
                                        sum(P_ij_k[t, l, k] for l in Omega_sending[i], k in K)
        [t = T, i = Ns], -Q_sub[t, i] == sum(Q_ij_k[t, l, k] - X[l, k] * I_sqr_k[t, l, k] for l in Omega_receiving[i], k in K) -
                                        sum(Q_ij_k[t, l, k] for l in Omega_sending[i], k in K)
        [t = T, i = Nu], p_imp[t, i] - p_exp[t, i] == sum(P_ij_k[t, l, k] - R[l, k] * I_sqr_k[t, l, k] for l in Omega_receiving[i], k in K) -
                                                    sum(P_ij_k[t, l, k] for l in Omega_sending[i], k in K)
        [t = T, i = Nu], q_imp[t, i] - q_exp[t, i] == sum(Q_ij_k[t, l, k] - X[l, k] * I_sqr_k[t, l, k] for l in Omega_receiving[i], k in K) -
                                                    sum(Q_ij_k[t, l, k] for l in Omega_sending[i], k in K)

        # Distflow
        [t = T, l = L], V_sqr[t, line_ends[l][2]] - V_sqr[t, line_ends[l][1]] <=
                        -sum(2 * (R[l, k] * P_ij_k[t, l, k] + X[l, k] * Q_ij_k[t, l, k]) +
                            (R[l, k]^2 + X[l, k]^2) * I_sqr_k[t, l, k] for k in K) +
                        M * (1 - Y[l])

        [t = T, l = L], V_sqr[t, line_ends[l][2]] - V_sqr[t, line_ends[l][1]] >=
                        -sum(2 * (R[l, k] * P_ij_k[t, l, k] + X[l, k] * Q_ij_k[t, l, k]) +
                            (R[l, k]^2 + X[l, k]^2) * I_sqr_k[t, l, k] for k in K) -
                        M * (1 - Y[l])

        [t = T, l = L], I_sqr[t, l] == sum(I_sqr_k[t, l, k] for k in K)
        [t = T, l = L], P_ij[t, l] == sum(P_ij_k[t, l, k] for k in K)
        [t = T, l = L], Q_ij[t, l] == sum(Q_ij_k[t, l, k] for k in K)

        [t = T, l = L], [V_sqr[t, line_ends[l][1]] / 2; I_sqr[t, l]; P_ij[t, l]; Q_ij[t, l]] in RotatedSecondOrderCone()
        #[l = L, k = K, t = T], [V_sqr[line_ends[l][1], t] / 2; I_sqr_k[t, l, k]; P_ij_k[t, l, k]; Q_ij_k[t, l, k]] in RotatedSecondOrderCone()

        #[i = Nu], sum(Y[l] for l in Omega_sending[i]) + sum(Y[l] for l in Omega_receiving[i]) >= 1

        # Operational limits
        #lorsque slack est négative => intger variable must be one 
        #si slack positive => binary variable == 0
        #si slack négatife => binary variable == 1 (idée forcer la variable à valoir 1 si nega)
        #-slack <= M * bianry variable 
        [t = T, l = L, k = K], Slack[t, l, k] <= MAX_CURRENT[l, k]^2 * Curr_limit[t, l, k]
        [t = T, l = L, k = K], I_sqr_k[t, l, k] - Slack[t, l, k] <= MAX_CURRENT[l, k]^2 * Alpha[l, k] # indispensable mais relaxée
        #[l = L, k = K, t = T], I_sqr_k[t, l, k] <= MAX_CURRENT[l, k] * Alpha[l, k] # indispensable
        [t = T, l = L, k = K], P_ij_k[t, l, k] <= MAX_CURRENT[l, k] * MAX_VOLTAGE * Alpha[l, k] # indispensable
        [t = T, l = L, k = K], P_ij_k[t, l, k] >= -MAX_CURRENT[l, k] * MAX_VOLTAGE * Alpha[l, k] # indispensable
        [t = T, l = L, k = K], Q_ij_k[t, l, k] <= MAX_CURRENT[l, k] * MAX_VOLTAGE * Alpha[l, k] # indispensable
        [t = T, l = L, k = K], Q_ij_k[t, l, k] >= -MAX_CURRENT[l, k] * MAX_VOLTAGE * Alpha[l, k] # indispensable

        # Substation apparent power limits
        [t = T, i = Ns], [S_sub[t, i]; P_sub[t, i]; Q_sub[t, i]] in SecondOrderCone()
        [t = T, i = Ns], S_sub[t, i] <= S_rating_init[i] + S_sub_capa[i]
        [i = Ns], S_sub_capa[i] <= Beta[i] * S_rating_max[i]

        # Voltage references
        [t = T, i = Ns_init], V_sqr[t, i] == 1
        [t = T, i = Ns_notinit], V_sqr[t, i] - 1 <= (MAX_VOLTAGE^2 - 1) * (1 - Beta[i])
        [t = T, i = Ns_notinit, ], V_sqr[t, i] - 1 >= (MIN_VOLTAGE^2 - 1) * (1 - Beta[i])

        # Load over-satisfaction
        #[t in T], sum(P_sub[t, i] for i in Ns) + sum(p_exp[t, i] - p_imp[t, i] for i in Nu) >= 0
        #[t in T], sum(Q_sub[t, i] for i in Ns) + sum(q_exp[t, i] - q_imp[t, i] for i in Nu) >= 0

        # Radiality constraint
        [l = L], sum(Alpha[l, k] for k in K) == Y[l]
        [t = T], sum(Y[l] for l in L) == Nu_size

        # Bound on s_grid_max
        #[i=Nu], s_grid_max[i] <= sum(MAX_CURRENT[l, k] * MAX_VOLTAGE * Alpha[l, k] for l in Omega_receiving[i], k in K) + 
        #                        sum(MAX_CURRENT[l, k] * MAX_VOLTAGE * Alpha[l, k] for l in Omega_sending[i], k in K) # NEW : doesn't change anything to the solution except faster !
    end)
    if ACTIVATE_RADIALITY_CONSTRAINTS
        @constraint(UPPER_MODEL, single_comm_constr1[i=Nu],
            -sum(k_ij[l] for l in Omega_receiving[i])
            +
            sum(k_ij[l] for l in Omega_sending[i]) == -1)

        @constraint(UPPER_MODEL, single_comm_constr2[i=Ns],
            -sum(k_ij[l] for l in Omega_receiving[i])
            +
            sum(k_ij[l] for l in Omega_sending[i]) >= 0)

        @constraint(UPPER_MODEL, single_comm_constr3[i=Ns_init],
            -sum(k_ij[l] for l in Omega_receiving[i])
            +
            sum(k_ij[l] for l in Omega_sending[i]) <= Nu_size)


        @constraint(UPPER_MODEL, single_comm_constr5[i=Ns_notinit],
            -sum(k_ij[l] for l in Omega_receiving[i])
            +
            sum(k_ij[l] for l in Omega_sending[i]) <= Nu_size * Beta[i])

        @constraint(UPPER_MODEL,
            single_comm_constr6[l=L],
            k_ij[l] <= Nu_size * Y[l]
        )

        @constraint(UPPER_MODEL,
            single_comm_constr7[l=L],
            k_ij[l] >= -Nu_size * Y[l]
        )

    end

    # if MAX_CO2_BUDGET > 0
    #     @constraint(UPPER_MODEL), CO2_budget <= MAX_CO2_BUDGET)
    # end
    @constraints(LOWER_MODEL, begin
        [i = Nu], user_costs[i] == energy_costs[i] - energy_revenues[i] + PV_costs[i] + grid_costs[i]
        [i = Nu], PV_costs[i] == BASE_POWER * (s_conv_pv[i] * CONVERTER_COST / AMORTIZATION_CONVERTER + p_pv_max[i] * PV_COST / AMORTIZATION_PV)
        [i = Nu], energy_costs[i] == BASE_POWER * sum(p_imp[t, i] * IMP_ELECTRICITY_ENRG_COST for t in T) * TIME_STEP * DAYS_A_YEAR / NB_PROFILES
        [i = Nu], energy_revenues[i] == BASE_POWER * sum(p_exp[t, i] * EXP_ELECTRICITY_ENRG_COST for t in T) * TIME_STEP * DAYS_A_YEAR / NB_PROFILES
        [i = Nu], grid_costs[i] == BASE_POWER * (s_grid_max[i] * GRID_CONNECTION_COST + sum((p_imp[t, i] * IMP_ELECTRICITY_DSO_COST + p_exp[t, i] * EXP_ELECTRICITY_DSO_COST) for t in T) * TIME_STEP * DAYS_A_YEAR / NB_PROFILES)
        [t = T, i = Nu], p_imp[t, i] - p_exp[t, i] == P_CONSUMPTION[t, i] - p_pv[t, i]
        [t = T, i = Nu], q_imp[t, i] - q_exp[t, i] == Q_CONSUMPTION[t, i] # - q_pv[t, i] DEACTIVATED in a first time to reimplement and test later.
        # [t = T, i = Nu], [s_grid_max[i], p_imp-p_exp[t, i], q_imp-q_exp[t, i]] in SecondOrderCone() # see next constraints
        [t = T, i = Nu], p_imp[t, i] <= s_grid_max[i]
        [t = T, i = Nu], q_imp[t, i] <= s_grid_max[i] # q_imp > 0
        [t = T, i = Nu], p_exp[t, i] <= s_grid_max[i]
        [t = T, i = Nu], q_exp[t, i] <= s_grid_max[i] # q_exp > 0

        # [t = T, i = Nu], q_pv[t, i] <= PV_MAX_Q * PV_PRODUCTION[t, i] * p_pv_max[i] # DEACTIVATED in a first time to reimplement and test later.
        # [t = T, i = Nu], -q_pv[t, i] <= PV_MAX_Q * PV_PRODUCTION[t, i] * p_pv_max[i] # DEACTIVATED in a first time to reimplement and test later.
        # [t = T, i = Nu], p_pv[t, i] <= PV_PRODUCTION[t, i] * p_pv_max[i] + PV_SLOPE * q_pv[t, i] # DEACTIVATED in a first time to reimplement and test later.
        # [t = T, i = Nu], p_pv[t, i] <= PV_PRODUCTION[t, i] * p_pv_max[i] - PV_SLOPE * q_pv[t, i] # DEACTIVATED in a first time to reimplement and test later.
        [t = T, i = Nu], p_pv[t, i] <= s_conv_pv[i] # The PV power output at time t is always bounded by the capacity of the converter
        [t = T, i = Nu], p_pv[t, i] <= PV_PRODUCTION[t, i] * p_pv_max[i] # The PV power output at time t is also bounded by the available sun power PV_PRODUCTION is in [W/wp]
        [t = T, i = Nu], p_exp[t, i] <= p_pv[t, i]  # NEW : because otherwise there are cases where p_exp > p_pv (Manon)
        #t = T, i = Nu], q_exp[t, i] <= q_pv[t, i] # NEW : because otherwise there are cases where q_exp > q_pv (Manon)
    end)

    # m2=BilevelJuMP._build_single_model(model)

    optimize!(model)


    #include("results_tables.jl")
    #include("plot_results.jl")

end