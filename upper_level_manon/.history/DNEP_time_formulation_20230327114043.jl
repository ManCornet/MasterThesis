#-----------------------------------------------------------------------------
#
#                   - TFE : Upper level problem formulation - 
#                             University of Liege
#
#-----------------------------------------------------------------------------
# Created By  : Manon Cornet
# Created Date: Tuesday March 14 2023
#
# upper_level_models:
#   File containing the upper level model functions
#
# =============================================================================
#                                   Imports
# =============================================================================
using JuMP, Gurobi
using DataFrames
using XLSX

include("utils.jl")

# =============================================================================
#                                   MODELS
# =============================================================================

# ======================= 1. Time-dependent formulation =======================

# Formulation that does not take into account the
function time_dependent_formulation(sets, costs, substation_param, conductor_param, 
                                    demand, delta_t)   

    # ATTENTION FOR THE FORMULATION MUST SPECIFY RANGE
    # ========================= Fetch the parameters ========================
    # Sets
    N, Ns, Ns_init, K, L, L_init, Y, T, Omega_sending, Omega_receiving = sets
    Ns_not_init = setdiff(Ns, Ns_init)
    Nu          = setdiff(N, Ns)
    L_not_init  = setdiff(L, L_init)

    N_size, Ns_size , K_size, L_size, Y_size, T_size = length(N), length(Ns), length(K), length(L), length(Y), length(T)
    # Costs
    sub_expan_cost, sub_install_cost, line_cost, losses_cost, DSO_INTEREST_RATE = costs

    # Substation parameters
    S_rating_init, S_rating_max = substation_param

    # Conductor parameters
    max_current, conductance, susceptance = conductor_param

    # Demand profiles
    P_D, tan_phi = demand
    Q_D = P_D .* tan_phi

    # Definition of weight terms of objective function
    DAYS_IN_A_YEAR = 365

    # ======================== Set up the Gurobi solver =======================
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "TimeLimit", 100)
    set_optimizer_attribute(model, "Presolve", 0)

    # ============================== Variables ================================ 
    # Other solution is to use this
    # Question: how to know the set of a variable to plot it ?
    #@variable(model, I_squared[time=T, line=L, conductor=K, container = Array) 
    @variable(model, I_squared[t=1:T_size, l=1:L_size, k=1:K_size])       # Squared current
    @variable(model, P_G[1:T_size, 1:N_size] >= 0)           # Active power generated
    @variable(model, Q_G[1:T_size, 1:N_size])                # Reactive power generated
    @variable(model, S_G[1:T_size, 1:N_size])                # Apparent power generated
    @variable(model, P_cond_forward[1:T_size, 1:L_size, 1:K_size])  # Forward direction active power flow
    @variable(model, P_cond_backward[1:T_size, 1:L_size, 1:K_size]) # Backward direction active power flow
    @variable(model, Q_cond_forward[1:T_size, 1:L_size, 1:K_size])  # Forward direction reactive power flow
    @variable(model, Q_cond_backward[1:T_size, 1:L_size, 1:K_size]) # Backward direction reactive power flow
    @variable(model, MIN_VOLTAGE^2 <= V_squared[1:T_size, 1:N_size] <= MAX_VOLTAGE^2)
    @variable(model, x[1:L_size], Bin)
    @variable(model, X_i_ij[1:T_size, 1:L_size, 1:K_size, 1:N_size] >= 0)
    @variable(model, X_ij_re[1:T_size, 1:L_size, 1:K_size] >= 0)
    @variable(model, X_ij_im[1:T_size, 1:L_size, 1:K_size])
    @variable(model, alpha[1:L_size, 1:K_size], Bin)
    @variable(model, beta[1:Ns_size], Bin)
    @variable(model, S_allocated[1:Ns_size] >= 0) 
    @variable(model, active_losses[1:T_size] >= 0) 

    for t in T, i in Nu
        fix(P_G[t, i], 0.0; force=true) 
        fix(Q_G[t, i], 0.0)
        fix(S_G[t, i], 0.0)
    end

    for t in T, l in L, k in K, i in N 
        if i in line_ends[l]
            continue
        end
        fix(X_i_ij[t, l, k, i], 0.0; force=true)
    end

    # ========================== Objective function =============================


    @objective(model, Min,  sum(alpha[l, k] * line_cost[k][l] * line_length[l] 
                                for k in K, l in L) 
                            + sum(S_allocated[i] * BASE_POWER * sub_install_cost for i in Ns_not_init)
                            + sum(S_allocated[i] * BASE_POWER * sub_expan_cost for i in Ns)
                            + DAYS_IN_A_YEAR * sum((1/(1 + DSO_INTEREST_RATE)^(t ÷ DAYS_IN_A_YEAR) * 
                              active_losses[t] * BASE_POWER * losses_cost * delta_t/60 for t in T))
    )

    # ============================== Constraints ================================ 

    # CONSTRAINT (2) -> In definition of alpha

    # CONSTRAINT (3)
    @constraint(model, line_constructed[l=L], x[l] == sum(alpha[l, k] for k in K))

    # CONSTRAINT (4)

    @constraint(model, 
    substation_apparent_power[t=T, i=Ns], 
        [S_G[t, i], P_G[t, i], Q_G[t, i]] in SecondOrderCone()
    )

    @constraint(model, 
                substation_capacity_limit[t=T, i=Ns], 
                S_G[t, i] <= S_rating_init[i] + S_allocated[i]
    )

    @constraint(model, 
                substation_capacity[i=Ns], 
                S_allocated[i] <= beta[i] * S_rating_max[i]
    ) 

    # CONSTRAINT (11)
    @constraint(model, number_of_lines, sum(x[l] for l in L) == length(N) - length(Ns))

    # CONSTRAINT (18)
    @constraint(model, 
                active_balance[t=T, i=N], 
                P_G[t, i] - P_D[t, i] == sum(P_cond_forward[t, l, k]  
                                        for l in Omega_sending[i], k in K)
                                        + sum(P_cond_backward[t, l, k] 
                                        for l in Omega_receiving[i], k in K)
    )

    # CONSTRAINT (19)
    @constraint(model, 
                reactive_balance[t=T, i=N], 
                Q_G[t, i] - Q_D[t, i] == sum(Q_cond_forward[t, l, k]
                                        for l in Omega_sending[i], k in K)
                                        + sum(Q_cond_backward[t, l, k]
                                        for l in Omega_receiving[i], k in K)
    )

    # CONSTRAINT (20)
    @constraint(model, 
                current_limit[t=T, l=L, k=K], 
                I_squared[t, l, k] <= alpha[l, k] * max_current[k][l]^2
    )

    for t in T, l in L, k in K 
        ifrom = line_ends[l][1]
        ito = line_ends[l][2]

        # CONSTRAINT (21)
        @constraint(model, P_cond_forward[t, l, k] 
                           == conductance[k][l] * (X_i_ij[t, l, k, ifrom] 
                                                   - X_ij_re[t, l, k])
                           - susceptance[k][l] * X_ij_im[t, l, k]
        )

        @constraint(model, P_cond_backward[t, l, k] 
                           == conductance[k][l] * (X_i_ij[t, l, k, ito] 
                                                   - X_ij_re[t, l, k])
                            - susceptance[k][l] * (- X_ij_im[t, l, k])) 

        # CONSTRAINT (22)
        @constraint(model, Q_cond_forward[t, l, k] 
                            == - susceptance[k][l] * (X_i_ij[t, l, k, ifrom] 
                                                     - X_ij_re[t, l, k])
                            - conductance[k][l] * X_ij_im[t, l, k])

        @constraint(model, Q_cond_backward[t, l, k] 
                            == - susceptance[k][l] * (X_i_ij[t, l, k, ito] 
                                                     - X_ij_re[t, l, k])
                            - conductance[k][l] * (- X_ij_im[t, l, k])) 


        # CONSTRAINT (23)
        @constraint(model, I_squared[t, l, k] 
                            == (conductance[k][l]^2 + susceptance[k][l]^2) 
                                * (X_i_ij[t, l, k, ifrom] + X_i_ij[t, l, k, ito] 
                                - 2 * X_ij_re[t, l, k]))


        # CONSTRAINT (24)
        @constraint(model, MIN_VOLTAGE^2 * alpha[l, k] <= X_i_ij[t, l, k, ifrom])
        @constraint(model, MAX_VOLTAGE^2 * alpha[l, k] >= X_i_ij[t, l, k, ifrom])
        @constraint(model, MIN_VOLTAGE^2 * alpha[l, k] <= X_i_ij[t, l, k, ito])
        @constraint(model, MAX_VOLTAGE^2 * alpha[l, k] >= X_i_ij[t, l, k, ito])

        # CONSTRAINT (25)
        @constraint(model, X_ij_re[t, l, k] <= MAX_VOLTAGE^2 * alpha[l, k])

        # CONSTRAINT (26)
        @constraint(model, X_ij_im[t, l, k] <= MAX_VOLTAGE^2 * alpha[l, k])
        @constraint(model, X_ij_im[t, l, k] >= - MAX_VOLTAGE^2 * alpha[l, k])

        # CONSTRAINT (27)
        @constraint(model, V_squared[t, ifrom] - X_i_ij[t, l, k, ifrom] 
                            >= MIN_VOLTAGE^2 * (1 - alpha[l, k]))

        @constraint(model, V_squared[t, ifrom] - X_i_ij[t, l, k, ifrom] 
                            <= MAX_VOLTAGE^2 * (1 - alpha[l, k]))

        @constraint(model, V_squared[t, ito] - X_i_ij[t, l, k, ito] 
                            >= MIN_VOLTAGE^2 * (1 - alpha[l, k]))

        @constraint(model, V_squared[t, ito] - X_i_ij[t, l, k, ito] 
                            <= MAX_VOLTAGE^2 * (1 - alpha[l, k]))
    end

    # CONSTRAINT (28)
    @constraint(model, 
                cone_28[t=T, l=L, k=K], 
                [X_i_ij[t, l, k, line_ends[l][1]] / 2, 
                X_i_ij[t, l, k, line_ends[l][2]],
                X_ij_re[t, l, k], X_ij_im[t, l, k]] 
                in RotatedSecondOrderCone())

    # CONSTRAINT (29)
    @constraint(model, 
                losses[t=T],
                active_losses[t] == sum(I_squared[t, l, k]/conductance[k][l] 
                                        for l in L, k in K))

    #print(model)
    optimize!(model)

    #solution_summary(model, verbose=false)

    if termination_status(model) == MOI.OPTIMAL
        return model

    elseif termination_status(model) == DUAL_INFEASIBLE
        println("problem unbounded")

    elseif termination_status(model) == MOI.INFEASIBLE
        println("problem infeasible")
    end
end