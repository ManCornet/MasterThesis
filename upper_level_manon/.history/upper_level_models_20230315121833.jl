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

# =============================================================================
#                                   Models
# =============================================================================


# =============================================================================
# ============================= MISOCP formulation ============================
# =============================================================================

function MISOCP_formulation(N, n_s, n_s_init, L, K ,
                            K_l, K_s, 
                            substation_utilization,
                            S_rating_init, S_rating_max,
                            conductance, susceptance, max_current, 
                            line_cost, line_ends, line_length, line_loss, 
                            interest_rate_losses, interest_rate_substation,
                            substation_op_cost, substation_fixed_cost,cost_unit_loss,  
                            Omega_sending, Omega_receiving, 
                            P_demand, Q_demand
                        )       

    # ======================== Set up the Gurobi solver =======================
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "TimeLimit", 100)
    set_optimizer_attribute(model, "Presolve", 0)

    # ============================== Variables ================================ 

    @variable(model, I_squared[1:K, 1:L] >= 0) # squared current
    @variable(model, P_s[1:N] >= 0)
    @variable(model, Q_s[1:N])
    @variable(model, S_s[1:n_s]) # Auxiliary variable
    @variable(model, P_cond_forward[1:K, 1:L]) # Forward direction
    @variable(model, P_cond_backward[1:K, 1:L]) # Backward direction
    @variable(model, Q_cond_forward[1:K, 1:L]) # Forward direction
    @variable(model, Q_cond_backward[1:K, 1:L]) # Backward direction
    @variable(model, MIN_VOLTAGE^2 <= V_squared[1:N] <= MAX_VOLTAGE^2)
    @variable(model, x[1:L], Bin)
    @variable(model, X_i_ij[1:K, 1:N, 1:L] >= 0) # TODO only for l connected to i, normally, because all others are zero.
    @variable(model, X_ij_re[1:K, 1:L] >= 0)
    @variable(model, X_ij_im[1:K, 1:L])
    @variable(model, alpha[1:K, 1:L], Bin)
    @variable(model, beta[1:N], Bin)
    @variable(model, active_losses[1:K, 1:L]) 
    @variable(model, reactive_losses[1:K, 1:L]) 
    @variable(model, losses1) 
    @variable(model, losses2) 
    @variable(model, epigraph)



    # It is easier to define P_s and Q_s for all nodes although they should be 
    # zero where it is not possible to put a substation.
    # force = true required to override the bound >= 0 given in P_s definition

    for i = n_s+1:N
        fix(P_s[i], 0.0; force=true) 
        fix(Q_s[i], 0.0)
    end

    for l in 1:L, k in 1:K, i in 1:N
        if i in line_ends[l]
            continue
        end
        fix(X_i_ij[k, i, l], 0.0; force=true)
    end

    # ============================= Constraints ===============================

    # CONSTRAINT (1) -> means constraint (1) in paper from Jabr (Polyhedral formulations ...)

   @objective(model, Min,  K_l               * sum(alpha[k, l] * line_cost[k][l] 
                                             * line_length[l] for k in 1:K, l in 1:L) 
                          + K_s              * sum(beta[i] * substation_fixed_cost[i] for i in 1:n_s)
                          + HOURS_IN_A_YEAR * (1 + interest_rate_losses) 
                                             * line_loss * cost_unit_loss 
                                             * sum((P_s[i] - P_demand[i]) * BASE_POWER for i in 1:N) 
                          +  HOURS_IN_A_YEAR * (1 + interest_rate_substation) 
                                             * substation_utilization * sum(substation_op_cost[i]
                                             * (P_s[i]^2 + Q_s[i]^2) * BASE_POWER^2 for i in 1:n_s)
                          + epigraph
                )

    @constraint(model, 
                    epigraph >= HOURS_IN_A_YEAR * (1 + interest_rate_substation) 
                    * substation_utilization * sum(substation_op_cost[i]
                    * (P_s[i]^2 + Q_s[i]^2) * BASE_POWER^2 for i in 1:n_s)
    )

    @constraint(model, 
                active_losses_constraints[k=1:K, l=1:L],
                active_losses[k, l] == I_squared[k, l]/conductance[k][l]
                )
    @constraint(model, 
                reactive_losses_constraints[k=1:K, l=1:L],
                reactive_losses[k, l] == I_squared[k, l]/susceptance[k][l]
            )

    @constraint(model, 
                losses_constraint1,
                losses1 == sum(I_squared[k, l]/conductance[k][l] for k in 1:K, l in 1:L)
                )
    @constraint(model, 
                losses_constraint2,
                losses2 == sum((P_s[i] - P_demand[i]) for i in 1:N)
                )

    # CONSTRAINT (2) -> In definition of alpha

    # CONSTRAINT (3)
    @constraint(model, line_constructed[l=1:L], x[l] == sum(alpha[k, l] for k in 1:K))

    # CONSTRAINT (4)
    @constraint(model, 
                substation_capacity[i=1:n_s], 
                S_s[i] == S_rating_init[i] + beta[i] * S_rating_max[i]
                ) 

    @constraint(model, 
                substation_capacity_limit[i=1:n_s], 
                [S_s[i], P_s[i], Q_s[i]] in SecondOrderCone()
                )

    # CONSTRAINT (11)
    @constraint(model, number_of_lines, sum(x[l] for l in 1:L) == N - n_s)

    # CONSTRAINT (18)
    @constraint(model, 
                active_balance[i=1:N], 
                P_s[i] - P_demand[i] ==   sum(P_cond_forward[k, l]
                                            for k in 1:K, l in Omega_sending[i])
                                        + sum(P_cond_backward[k, l] 
                                            for k in 1:K, l in Omega_receiving[i])
                )
    # CONSTRAINT (19)

    @constraint(model, 
                reactive_balance[i=1:N], 
                Q_s[i] - Q_demand[i] ==   sum(Q_cond_forward[k, l]
                                            for k in 1:K, l in Omega_sending[i])
                                        + sum(Q_cond_backward[k, l]
                                            for k in 1:K, l in Omega_receiving[i])
                )

    # CONSTRAINT (20)
    @constraint(model, 
                current_limit[k=1:K, l=1:L], 
                I_squared[k, l] <= alpha[k, l] * max_current[k][l]^2
                )

    for k in 1:K, l in 1:L
        ifrom = line_ends[l][1]
        ito = line_ends[l][2]

        # CONSTRAINT (21)
        @constraint(model, P_cond_forward[k, l] == conductance[k][l] * (X_i_ij[k, ifrom, l] - X_ij_re[k, l])
                                                    -
                                                    susceptance[k][l] * X_ij_im[k, l])

        # CONSTRAINT (21)
        @constraint(model, P_cond_backward[k, l] == conductance[k][l] * (X_i_ij[k, ito, l] - X_ij_re[k, l])
                                                     -
                                                    susceptance[k][l] * (-X_ij_im[k, l])) # X is hermitian => X_ij_im = - X_ji_im but we define only X_ij_im for simplicity.
    

        # CONSTRAINT (22)
        @constraint(model, 
                    Q_cond_forward[k,l] == - susceptance[k][l] 
                                           * (X_i_ij[k, ifrom, l] - X_ij_re[k, l]) 
                                           - conductance[k][l] * X_ij_im[k, l]
                    )

        
        # CONSTRAINT (22)
        @constraint(model, 
                    Q_cond_backward[k,l] == - susceptance[k][l] 
                                            * (X_i_ij[k, ito, l] - X_ij_re[k, l]) 
                                            - conductance[k][l] * (- X_ij_im[k, l])
                    )

        # CONSTRAINT (23)
        @constraint(model, 
        I_squared[k, l] ==    (conductance[k][l]^2 + susceptance[k][l]^2) 
                            * (X_i_ij[k, ifrom, l] + X_i_ij[k, ito, l] 
                            - 2 * X_ij_re[k, l])
        )
        # CONSTRAINT (24)
        @constraint(model, MIN_VOLTAGE^2 * alpha[k, l] <= X_i_ij[k, ifrom, l])
        @constraint(model, MAX_VOLTAGE^2 * alpha[k, l] >= X_i_ij[k, ifrom, l])
        @constraint(model, MIN_VOLTAGE^2 * alpha[k, l] <= X_i_ij[k, ito, l])
        @constraint(model, MAX_VOLTAGE^2 * alpha[k, l] >= X_i_ij[k, ito, l])

        # CONSTRAINT (25)
        @constraint(model, X_ij_re[k, l] <= MAX_VOLTAGE^2 * alpha[k, l])

        # CONSTRAINT (26)
        @constraint(model, X_ij_im[k, l] <= MAX_VOLTAGE^2 * alpha[k, l])
        @constraint(model, X_ij_im[k, l] >= - MAX_VOLTAGE^2 * alpha[k, l])

        # CONSTRAINT (27)
        @constraint(model,   V_squared[ifrom] - X_i_ij[k, ifrom, l] 
                          >= MIN_VOLTAGE^2 * (1 - alpha[k, l])
        )
        @constraint(model,   V_squared[ifrom] - X_i_ij[k, ifrom, l] 
                          <= MAX_VOLTAGE^2 * (1 - alpha[k, l])
        )
        @constraint(model,   V_squared[ito] - X_i_ij[k, ito, l] 
                          >= MIN_VOLTAGE^2 * (1 - alpha[k, l])
        )
        @constraint(model,   V_squared[ito] - X_i_ij[k, ito, l] 
                          <= MAX_VOLTAGE^2 * (1 - alpha[k, l])
        )
    end

    # CONSTRAINT (28)
    @constraint(model, cone_28[k=1:K, l=1:L], [X_i_ij[k, line_ends[l][1], l] / 2, X_i_ij[k, line_ends[l][2], l],
                X_ij_re[k, l], X_ij_im[k, l]] in RotatedSecondOrderCone())

    print(model)

    optimize!(model)

    solution_summary(model, verbose=true)

    if termination_status(model) == MOI.OPTIMAL
        I_squared   = transpose(value.(I_squared))
        V_squared   = value.(V_squared)
        x           = value.(x)
        alpha       = transpose(value.(alpha))
        beta        = value.(beta)
        X_ij_re     = transpose(value.(X_ij_re))
        X_ij_im     = transpose(value.(X_ij_im))
        X_i_ij      = value.(X_i_ij)
        P_s         = value.(P_s)
        Q_s         = value.(Q_s)
        P_cond_forward = transpose(value.(P_cond_forward))
        P_cond_backward = transpose(value.(P_cond_backward))
        Q_cond_forward = transpose(value.(Q_cond_forward))
        Q_cond_backward = transpose(value.(Q_cond_backward))
        losses1      = value.(losses1)
        losses2      = value.(losses2)
        obj         = objective_value(model)
        time        = solve_time(model)

        return I_squared, V_squared, x, alpha, beta, X_ij_re, 
        X_ij_im, X_i_ij, P_s, Q_s, P_cond_forward, Q_cond_forward,
        P_cond_backward, Q_cond_backward, obj, time, losses1, losses2

    elseif termination_status(model) == DUAL_INFEASIBLE
        println("problem unbounded")

    elseif termination_status(model) == MOI.INFEASIBLE
        println("problem infeasible")
    end

end