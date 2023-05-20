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
include("utils.jl")

# =============================================================================
#                                   MODELS
# =============================================================================

# ============================= 1. Jabr formulation ===========================

function Jabr_formulation(N, n_s, n_s_init, L, K ,
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

    @variable(model, I_squared[1:K, 1:L]) # squared current
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
                                             * sum((P_s[i] - P_demand[i]) for i in 1:N) 
                          +  HOURS_IN_A_YEAR * (1 + interest_rate_substation) 
                                             * substation_utilization * sum(substation_op_cost[i]
                                             * (P_s[i]^2 + Q_s[i]^2) for i in 1:n_s)
                )

    #@constraint(model, 
    #                epigraph >= HOURS_IN_A_YEAR * (1 + interest_rate_substation) 
    #                * substation_utilization * sum(substation_op_cost[i]
    #                * (P_s[i]^2 + Q_s[i]^2) * BASE_POWER^2 for i in 1:n_s)
    #)

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

    # Fix the substation voltage
    #=
    @constraint(model, 
                substation_voltage[i=1:n_s],
                V_squared[i] == 1
                )
    =#       
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
        @constraint(model, Q_cond_forward[k, l] == -susceptance[k][l] * (X_i_ij[k, ifrom, l] - X_ij_re[k, l])
                                              -
                                              conductance[k][l] * X_ij_im[k, l])
        @constraint(model, Q_cond_backward[k, l] == -susceptance[k][l] * (X_i_ij[k, ito, l] - X_ij_re[k, l])
                                              -
                                              conductance[k][l] * (-X_ij_im[k, l])) # X is hermitian => X_ij_im = - X_ji_im but we define only X_ij_im for simplicity.
    

        # CONSTRAINT (23)
        @constraint(model, I_squared[k, l] == (conductance[k][l]^2 + susceptance[k][l]^2) * (X_i_ij[k, ifrom, l] + X_i_ij[k, ito, l] - 2 * X_ij_re[k, l]))

        
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
        include("export_xlsx.jl")
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

        check_rotated_cones(cone_28)
        check_cones(substation_capacity_limit) 

        return I_squared, V_squared, x, alpha, beta, X_ij_re, 
        X_ij_im, X_i_ij, P_s, Q_s, P_cond_forward, Q_cond_forward,
        P_cond_backward, Q_cond_backward, obj, time, losses1, losses2

    elseif termination_status(model) == DUAL_INFEASIBLE
        println("problem unbounded")

    elseif termination_status(model) == MOI.INFEASIBLE
        println("problem infeasible")
    end

end



# Formulation that does not take into account the
function time_dependent_formulation(sets, costs, substation_param, conductor_param, 
                                    demand, delta_t)   

    # ========================= Fetch the parameters ========================
    # Sets
    N, Ns, Ns_init, K, L, L_init, Y, T, Omega_sending, Omega_receiving = sets
    Ns_not_init = setdiff(Ns, Ns_init)
    Nu          = setdiff(N, Ns)
    L_not_init  = setdiff(L, L_init)

    # Costs
    sub_expan_cost, sub_install_cost, line_cost, losses_cost = costs

    # Substation parameters
    S_rating_init, S_rating_max = substation_param

    # Conductor parameters
    max_current, conductance, susceptance = conductor_param

    # Demand profiles
    demand_profiles, tan_phi = demand
    P_D = [zeros(length(Ns)); demand_profiles .* 1e3 ./ BASE_POWER]
    Q_D = P_D .* tan_phi
 
    # Definition of weight terms of objective function
    DAYS_IN_A_YEAR = 365

    # ======================== Set up the Gurobi solver =======================
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "TimeLimit", 100)
    set_optimizer_attribute(model, "Presolve", 0)

    # ============================== Variables ================================ 

    @variable(model, I_squared[Y, T, L, K])       # Squared current
    @variable(model, P_G[Y, T, N] >= 0)           # Active power generated
    @variable(model, Q_G[Y, T, N])                # Reactive power generated
    @variable(model, S_G[Y, T, N])                # Apparent power generated
    @variable(model, P_cond_forward[Y, T, L, K])  # Forward direction active power flow
    @variable(model, P_cond_backward[Y, T, L, K]) # Backward direction active power flow
    @variable(model, Q_cond_forward[Y, T, L, K])  # Forward direction reactive power flow
    @variable(model, Q_cond_backward[Y, T, L, K]) # Backward direction reactive power flow
    @variable(model, MIN_VOLTAGE^2 <= V_squared[Y, T, N] <= MAX_VOLTAGE^2)
    @variable(model, x[L], Bin)
    @variable(model, X_i_ij[Y, T, L, K, N] >= 0)
    @variable(model, X_ij_re[Y, T, L, K] >= 0)
    @variable(model, X_ij_im[Y, T, L, K])
    @variable(model, alpha[L, K], Bin)
    @variable(model, beta[Ns], Bin)
    @variable(model, S_allocated[Ns] >= 0) 
    @variable(model, active_losses[Y, T] >= 0) 

    for i in Nu
        fix(P_G[i], 0.0; force=true) 
        fix(Q_G[i], 0.0)
        fix(S_G[i], 0.0)
    end

    for y in Y, t in T, l in L, k in K, i in N 
        if i in line_ends[l]
            continue
        end
        fix(X_i_ij[y, t, l, k, i], 0.0; force=true)
    end

    # ========================== Objective function =============================
    

    @objective(model, Min,  sum(alpha[l, k] * line_cost[k][l] * line_length[l] 
                                for k in K, l in L) 
                          + sum(S_allocated[i] * BASE_POWER * sub_install_cost for i in Ns_not_init)
                          + sum(S_allocated[i] * BASE_POWER * sub_expan_cost for i in Ns)
                          + DAYS_IN_A_YEAR * sum((1/(1 + DSO_INTEREST_RATE)^(y-1) 
                          * sum(active_losses[y, t] * BASE_POWER * losses_cost * delta_t/60 for t in T) 
                                for y in Y))
                )

    # ============================== Constraints ================================ 

    # CONSTRAINT (2) -> In definition of alpha

    # CONSTRAINT (3)
    @constraint(model, line_constructed[l=L], x[l] == sum(alpha[l, k] for k in K))

    # CONSTRAINT (4)
   
    @constraint(model, 
                substation_capacity_limit[y=Y, t=T, i=Ns], 
                [S_G[y, t, i], P_G[y, t, i], Q_G[y, t, i]] in SecondOrderCone()
                )

    @constraint(model, 
                substation_capacity_limit[y=Y, t=T, i=Ns], 
                S_G[y, t, i] <= S_rating_init[i] + S_allocated[i]
                )

    @constraint(model, 
                substation_capacity[i=Ns], 
                S_allocated[i] <= beta[i] * S_rating_max[i]
                ) 

    # CONSTRAINT (11)
    @constraint(model, number_of_lines, sum(x[l] for l in L) == length(N) - length(Ns))

    # CONSTRAINT (18)
    @constraint(model, 
                active_balance[y=Y, t=T, i=N], 
                P_G[y, t, i] - P_D[y, t, i] == sum(P_cond_forward[y, t, l, k]  
                                                 for k in K, l in Omega_sending[i])
                                                 + sum(P_cond_backward[y, t, l, k] 
                                                 for k in K, l in Omega_receiving[i])
                )

    # CONSTRAINT (19)
    @constraint(model, 
                reactive_balance[y=Y, t=T, i=N], 
                Q_G[y, t, i] - Q_D[y, t, i] == sum(Q_cond_forward[y, t, l, k]
                                                for k in K, l in Omega_sending[i])
                                                + sum(Q_cond_backward[y, t, l, k]
                                                for k in K, l in Omega_receiving[i])
                )

    # CONSTRAINT (20)
    @constraint(model, 
                current_limit[y=Y, t=T, l=L, k=K], 
                I_squared[y, t, l, k] <= alpha[l, k] * max_current[k][l]^2
                )

    for k in K, l in L, y in Y, t in T
        ifrom = line_ends[l][1]
        ito = line_ends[l][2]

        # CONSTRAINT (21)
        @constraint(model, P_cond_forward[k, l, y, t] 
                           == 
                           conductance[k][l] * (X_i_ij[k, ifrom, l, y, t] - X_ij_re[k, l, y, t])
                           - susceptance[k][l] * X_ij_im[k, l, y, t]
                    )

        @constraint(model, P_cond_backward[k, l, y, t] 
                           ==
                           conductance[k][l] * (X_i_ij[k, ito, l, y, t] - X_ij_re[k, l, y, t])
                           - susceptance[k][l] * (- X_ij_im[k, l, y, t])) 
    
        # CONSTRAINT (22)
        @constraint(model, Q_cond_forward[k, l, y, t] 
                           == 
                           - susceptance[k][l] * (X_i_ij[k, ifrom, l, y, t] - X_ij_re[k, l, y, t])
                           - conductance[k][l] * X_ij_im[k, l, y, t])

        @constraint(model, Q_cond_backward[k, l, y, t] 
                           == 
                           - susceptance[k][l] * (X_i_ij[k, ito, l, y, t] - X_ij_re[k, l, y, t])
                           - conductance[k][l] * (- X_ij_im[k, l, y, t])) 
    

        # CONSTRAINT (23)
        @constraint(model, I_squared[k, l, y, t] 
                           == 
                           (conductance[k][l]^2 + susceptance[k][l]^2) 
                           * (X_i_ij[k, ifrom, l, y, t] + X_i_ij[k, ito, l, y, t] 
                           - 2 * X_ij_re[k, l, y, t]))

        
        # CONSTRAINT (24)
        @constraint(model, MIN_VOLTAGE^2 * alpha[k, l] <= X_i_ij[k, ifrom, l, y, t])
        @constraint(model, MAX_VOLTAGE^2 * alpha[k, l] >= X_i_ij[k, ifrom, l, y, t])
        @constraint(model, MIN_VOLTAGE^2 * alpha[k, l] <= X_i_ij[k, ito, l, y, t])
        @constraint(model, MAX_VOLTAGE^2 * alpha[k, l] >= X_i_ij[k, ito, l, y, t])

        # CONSTRAINT (25)
        @constraint(model, X_ij_re[k, l, y, t] <= MAX_VOLTAGE^2 * alpha[k, l])

        # CONSTRAINT (26)
        @constraint(model, X_ij_im[k, l, y, t] <= MAX_VOLTAGE^2 * alpha[k, l])
        @constraint(model, X_ij_im[k, l, y, t] >= - MAX_VOLTAGE^2 * alpha[k, l])

        # CONSTRAINT (27)
        @constraint(model,   V_squared[ifrom, y, t] - X_i_ij[k, ifrom, l, y, t] 
                          >= MIN_VOLTAGE^2 * (1 - alpha[k, l])
        )
        @constraint(model,   V_squared[ifrom, y, t] - X_i_ij[k, ifrom, l, y, t] 
                          <= MAX_VOLTAGE^2 * (1 - alpha[k, l])
        )
        @constraint(model,   V_squared[ito, y, t] - X_i_ij[k, ito, l, y, t] 
                          >= MIN_VOLTAGE^2 * (1 - alpha[k, l])
        )
        @constraint(model,   V_squared[ito, y, t] - X_i_ij[k, ito, l, y, t] 
                          <= MAX_VOLTAGE^2 * (1 - alpha[k, l])
        )
    end

    # CONSTRAINT (28)
    @constraint(model, 
                cone_28[k=K, l=L, y=Y, t=T], 
                [X_i_ij[k, line_ends[l][1], l, y, t] / 2, 
                 X_i_ij[k, line_ends[l][2], l, y, y],
                 X_ij_re[k, l, y, t], X_ij_im[k, l, y, t]] 
                 in RotatedSecondOrderCone())

    # CONSTRAINT (29)
    @constraint(model, 
               losses[y=Y, t=T],
               active_losses[y, t] 
               == sum(I_squared[k, l, y, t]/conductance[k][l] for k in K, l in L)
    )
    print(model)

    optimize!(model)

    solution_summary(model, verbose=true)

    if termination_status(model) == MOI.OPTIMAL
        I_squared       = value.(I_squared)
        V_squared       = value.(V_squared)
        x               = value.(x)
        alpha           = value.(alpha)
        beta            = value.(beta)
        X_ij_re         = value.(X_ij_re)
        X_ij_im         = value.(X_ij_im)
        X_i_ij          = value.(X_i_ij)
        P_G             = value.(P_G)
        Q_G             = value.(Q_G)
        S_G             = value.(S_G)
        S_allocated     = value.(S_allocated)
        P_cond_forward  = value.(P_cond_forward)
        P_cond_backward = value.(P_cond_backward)
        Q_cond_forward  = value.(Q_cond_forward)
        Q_cond_backward = value.(Q_cond_backward)
        active_losses   = value.(active_losses)
        obj             = objective_value(model)
        time            = solve_time(model)

        check_rotated_cones(cone_28)
        check_cones(substation_capacity_limit) 

        return I_squared, V_squared, x, alpha, beta, X_ij_re, 
        X_ij_im, X_i_ij, P_G, Q_G, S_G, S_allocated, P_cond_forward, 
        Q_cond_forward, P_cond_backward, Q_cond_backward, active_losses,
        obj, time

    elseif termination_status(model) == DUAL_INFEASIBLE
        println("problem unbounded")

    elseif termination_status(model) == MOI.INFEASIBLE
        println("problem infeasible")
    end

  
end
