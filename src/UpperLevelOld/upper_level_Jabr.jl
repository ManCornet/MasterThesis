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

function check_rotated_cones(constraint)
    println("Rotated cone constraints that are not tight:")
    for k in 1:K, l in 1:K
        x = value(constraint[k, l])
        slack = 2 * x[1] * x[2] - (x[3]^2 + x[4]^2)
        if abs(slack) > 1e-3
            println(k, l)
        end
    end
end

function check_cones(constraint)
    println("Cone constraints that are not tight:")
    for n in 1:n_s
        x = value(constraint[n])
        slack = x[1]^2 - (x[2]^2 + x[3]^2)
        if abs(slack) > 1e-3
            println(n)
        end
    end
end
# =============================================================================
# ============================= MISOCP formulation ============================
# =============================================================================

function UL_Jabr(sets, costs, substation_param, conductor_param, losses, others, demand)       

    # ========================= Fetch the parameters ========================
    # Sets
    N, Ns, K, L, Omega_sending, Omega_receiving = sets
    Nu          = setdiff(N, Ns)

    N_size, Ns_size , K_size, L_size = length(N), length(Ns), length(K), length(L)
    # Costs
    substation_fixed_cost, substation_op_cost, line_cost, cost_unit_loss = costs

    # Losses 
    line_loss, interest_rate_losses = losses

    # Substation parameters
    S_rating_init, S_rating_max, substation_utilization, interest_rate_substation = substation_param

    # Conductor parameters
    max_current, conductance, susceptance = conductor_param

    # Demand profiles
    P_D, Q_D = demand
    println(P_D)
    println(Q_D)
    K_s, K_l = others

    # ======================== Set up the Gurobi solver =======================
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "TimeLimit", 100)
    set_optimizer_attribute(model, "Presolve", 0)

    # ============================== Variables ================================ 

    @variables( model,   
                begin 
                MIN_VOLTAGE^2 <= V_squared[1:N_size] <= MAX_VOLTAGE^2
                I_squared[1:L_size, 1:K_size]
                P_ij[1:L_size, 1:K_size] # forward direction
                P_ji[1:L_size, 1:K_size] # backward direction
                Q_ij[1:L_size, 1:K_size]  # forward direction
                Q_ji[1:L_size, 1:K_size] # backward direction
                X_i_ij[1:L_size, 1:K_size, 1:N_size]
                X_ij_re[1:L_size, 1:K_size]             >= 0
                X_ij_im[1:L_size, 1:K_size]
                P_G[1:N_size]                           >= 0
                Q_G[1:N_size]
                S_G[1:N_size]                          
                alpha[1:L_size, 1:K_size]               , Bin
                beta[1:Ns_size]                         , Bin 
                x[1:L_size]                             , Bin
                end
            )

    for i in Nu
        fix(S_G[i], 0.0; force=true)
        fix(P_G[i], 0.0; force=true) 
        fix(Q_G[i], 0.0)
    end

    for l in L, k in K, i in N 
        if i in line_ends[l]
            continue
        else
            fix(X_i_ij[l, k, i], 0.0; force=true)
        end
    end
    


    # ============================= Constraints ===============================

    # CONSTRAINT (1) -> means constraint (1) in paper from Jabr (Polyhedral formulations ...)

   @objective(model, Min,   K_l             * sum(alpha[l, k] * line_cost[l, k] for l in L, k in K) 
                          + K_s             * sum(beta[i] * substation_fixed_cost[i] for i in Ns)
                          + HOURS_IN_A_YEAR * (1 + interest_rate_losses) 
                                            * line_loss * cost_unit_loss * BASE_POWER
                                            * sum((P_G[i] - P_D[i]) for i in N) 
                          + HOURS_IN_A_YEAR * (1 + interest_rate_substation) 
                                            * substation_utilization * sum(substation_op_cost[i]
                                            * (P_G[i]^2 + Q_G[i]^2) * BASE_POWER^2 for i in Ns)
                )

    # CONSTRAINT (2) -> In definition of alpha

    # CONSTRAINT (3)
    @constraint(model, line_constructed[l=L], x[l] == sum(alpha[l, k] for k in K))

    # CONSTRAINT (4)
    @constraint(model, 
                substation_capacity[i=Ns], 
                S_G[i] == S_rating_init[i] + beta[i] * S_rating_max[i]
                ) 

    @constraint(model, 
                substation_capacity_limit[i=Ns], 
                [S_G[i]; P_G[i]; Q_G[i]] in SecondOrderCone()
                )

    # CONSTRAINT (11)
    @constraint(model, number_of_lines, sum(x[l] for l in L) ==  N_size - Ns_size)

    # CONSTRAINT (18)
    @constraint(model, 
                active_balance[i=N], 
                P_G[i] - P_D[i] ==   sum(P_ij[l, k]  
                                    for l in Omega_sending[i], k in K)
                                    + sum(P_ji[l, k] 
                                    for l in Omega_receiving[i], k in K)
                )
    # CONSTRAINT (19)

    @constraint(model, 
                reactive_balance[i=N], 
                Q_G[i] - Q_D[i] == sum(Q_ij[l, k]
                                        for l in Omega_sending[i], k in K)
                                        + sum(Q_ji[l, k]
                                        for l in Omega_receiving[i], k in K)
                )

     # CONSTRAINT (20)
     @constraint(model, 
                current_limit[l=L, k=K], 
                I_squared[l, k] <= alpha[l, k] * max_current[l, k]^2
                )

    for l in L, k in K
        ifrom = line_ends[l][1]
        ito = line_ends[l][2]

        # CONSTRAINT (21)
        @constraint(model, P_ij[l, k] 
                           == conductance[l, k] * (X_i_ij[l, k, ifrom] 
                                                   - X_ij_re[l, k])
                           + susceptance[l, k] * X_ij_im[l, k]
        )

        
        @constraint(model, P_ji[l, k] 
                           == conductance[l, k] * (X_i_ij[l, k, ito] 
                                                   - X_ij_re[l, k])
                            - susceptance[l, k] * X_ij_im[l, k]) 
                            
        # This constraint does not work
        #=
        @constraint(model, -P_ji[l, k] 
                        == P_ij[l, k] + I_squared[l,k] * resistance[l, k]) 
        =#
        # CONSTRAINT (22)
        @constraint(model, Q_ij[l, k] 
                            == susceptance[l, k] * (X_i_ij[l, k, ifrom] 
                                                     - X_ij_re[l, k])
                            - conductance[l, k] * X_ij_im[l, k])

        
        @constraint(model, Q_ji[l, k] 
                            ==  susceptance[l, k] * (X_i_ij[l, k, ito] 
                                                     - X_ij_re[l, k])
                             + conductance[l, k] * (X_ij_im[l, k])) 
    
        # this constraint does not work because since we minimize I 
        # the problem wants to 
        #=
        @constraint(model, -Q_ji[l, k] 
                            == Q_ij[l, k] + I_squared[l,k]/susceptance[l,k]) 
        =#
        # CONSTRAINT (23)
        @constraint(model, I_squared[l, k] 
                            == (conductance[l, k]^2 + susceptance[l, k]^2) 
                                * (X_i_ij[l, k, ifrom] + X_i_ij[l, k, ito] 
                                - 2 * X_ij_re[l, k]))

        
        # CONSTRAINT (24)
        @constraint(model, MIN_VOLTAGE^2 * alpha[l, k] <= X_i_ij[l, k, ifrom])
        @constraint(model, MAX_VOLTAGE^2 * alpha[l, k] >= X_i_ij[l, k, ifrom])
        @constraint(model, MIN_VOLTAGE^2 * alpha[l, k] <= X_i_ij[l, k, ito])
        @constraint(model, MAX_VOLTAGE^2 * alpha[l, k] >= X_i_ij[l, k, ito])

        # CONSTRAINT (25)
        @constraint(model, X_ij_re[l, k] <= MAX_VOLTAGE^2 * alpha[l, k])

        # CONSTRAINT (26)
        @constraint(model, X_ij_im[l, k] <= MAX_VOLTAGE^2 * alpha[l, k])
        @constraint(model, X_ij_im[l, k] >= - MAX_VOLTAGE^2 * alpha[l, k])

        # CONSTRAINT (27)
        @constraint(model, V_squared[ifrom] - X_i_ij[l, k, ifrom] 
                            >= MIN_VOLTAGE^2 * (1 - alpha[l, k]))

        @constraint(model, V_squared[ifrom] - X_i_ij[l, k, ifrom] 
                            <= MAX_VOLTAGE^2 * (1 - alpha[l, k]))

        @constraint(model, V_squared[ito] - X_i_ij[l, k, ito] 
                            >= MIN_VOLTAGE^2 * (1 - alpha[l, k]))

        @constraint(model, V_squared[ito] - X_i_ij[l, k, ito] 
                            <= MAX_VOLTAGE^2 * (1 - alpha[l, k]))
    end


    # CONSTRAINT (28)
    @constraint(model, 
                cone_28[l=L, k=K], 
                [X_i_ij[l, k, line_ends[l][1]] / 2; 
                X_i_ij[l, k, line_ends[l][2]];
                X_ij_re[l, k]; 
                X_ij_im[l, k]] 
                in RotatedSecondOrderCone())
    
    # Additional constraints

    # CONSTRAINT (29)
    @constraint(model, 
                power_balance1, 
                sum(P_G[i] - P_D[i] for i in N) >= 0
    )

    @constraint(model, 
                power_balance2, 
                sum(Q_G[i] - Q_D[i] for i in N) >= 0
    )

    # Reference voltage of slack node
    @constraint(model, 
                ref_voltage1[i=Ns], 
                V_squared[i] - 1 <= (MAX_VOLTAGE^2-1)*(1-beta[i])
    )

    @constraint(model, 
                ref_voltage2[i=Ns], 
                V_squared[i] - 1 >= (MIN_VOLTAGE^2 - 1)*(1-beta[i])
    )
    
    # One possibility: force the powers to be equal to this
    
    @constraint(model, 
                losses1[l=L], 
                sum(I_squared[l, k] * resistance[l,k] for k in K) == sum(P_ij[l, k] + P_ji[l, k] for k in K)
    )

    @constraint(model, 
                losses2[l=L], 
                sum(I_squared[l, k] * resistance[l,k] for k in K) == sum(Q_ij[l, k] + Q_ji[l, k] for k in K)
    )
    
    


    # CONSTRAINT (30)
    @constraint(model, 
                cone_30[l=L, k=K], 
                [V_squared[line_ends[l][1]] / 2;
                I_squared[l, k];  
                P_ij[l, k];
                Q_ij[l, k]] 
                in RotatedSecondOrderCone())

    # CONSTRAINT (31)
    @constraint(model, 
                cone_31[l=L, k=K], 
                [V_squared[line_ends[l][2]] / 2; 
                I_squared[l, k]; 
                P_ji[l, k];
                Q_ji[l, k]] 
                in RotatedSecondOrderCone())
   
    #print(model)

    optimize!(model)

    solution_summary(model, verbose=false)

    function check_rotated_cones()
        println("Rotated cone constraints that are not tight:")
        for l in L, k in K
            x = value(cone_28[l, k])
            println(2 * x[1] * x[2])
            println(x[3]^2 + x[4]^2)
            slack = 2 * x[1] * x[2] - (x[3]^2 + x[4]^2)
            if abs(slack) > 1e-6
                println(k, l)
            end
        end
    end
    
    function check_cones()
        println("Cone constraints that are not tight:")
        for n in Ns
            x = value(substation_capacity_limit[n])
            slack = x[1]^2 - (x[2]^2 + x[3]^2)
            if abs(slack) > 1e-6
                println(n)
            end
        end
    end

    function check_losses()
        I_2 = value.(I_squared)
        p_ij = value.(P_ij)
        p_ji = value.(P_ji)
        losses1 = Vector{Float64}(undef, L_size)
        losses2 = Vector{Float64}(undef, L_size)
        for l in L
            println(line_ends[l])
            losses1[l] = sum([I_2[l, k]/conductance[l, k] for k in K])
            losses2[l] = sum([abs(p_ij[l, k] + p_ji[l, k]) for k in K])
        end
       

        x_ij_re = sum(value.(X_ij_re)[3, :])
        x_i_ij = sum(value.(X_i_ij)[3, :, line_ends[3][1]])
        x_j_ij = sum(value.(X_i_ij)[3, :, line_ends[3][2]])

        println(conductance)
        println(susceptance)
        println(2*x_ij_re)
        println(x_i_ij^2 + x_j_ij^2)
        println(x_i_ij^2 + x_j_ij^2 - 2*x_ij_re)


        println("losses 1: $losses1")
        println("losses 2: $losses2")
       
    end
    
    check_rotated_cones()
    check_cones()
    check_losses()

    if termination_status(model) == MOI.OPTIMAL
        var_values = Dict(  string(k) => value.(v) for (k, v) 
                            in object_dictionary(model) 
                            if (v isa AbstractArray{VariableRef} || 
                                v isa VariableRef)
                        )

        
        var_sets = Dict("X_i_ij" => ["node"], 
                        "X_ij_re" => ["line", "conductor"],
                        "X_ij_im" => ["line", "conductor"],
                        "P_ij" => ["line", "conductor"],
                        "P_ji" => ["line", "conductor"],
                        "Q_ij" => ["line", "conductor"],
                        "Q_ji" => ["line", "conductor"],
                        "I_squared" => ["line", "conductor"],
                        "V_squared" => ["node"],
                        "P_G" => ["node"],
                        "Q_G" => ["node"],
                        "S_G" => ["node"],
                        "alpha" => ["line", "conductor"],
                        "beta" => ["node"],
                        "x" => ["line"]
                        )

        #check_rotated_cones(cone_28, L, K)
        #check_cones(substation_apparent_power, Ns)

        return var_values, var_sets 

    elseif termination_status(model) == DUAL_INFEASIBLE
        println("problem unbounded")

    elseif termination_status(model) == MOI.INFEASIBLE
        println("problem infeasible")
    end
end