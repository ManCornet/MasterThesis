#-----------------------------------------------------------------------------
#
#                   - TFE : Upper level problem formulation - 
#                             University of Liege
#
#-----------------------------------------------------------------------------
# Created By  : Manon Cornet
# Created Date: Friday April 10th 2023
#
# upper_level_models:
#   File containing the distflow upper level model function
#
# =============================================================================
#                                   Imports
# =============================================================================

using JuMP, Gurobi

# =============================================================================
#                                  Functions
# =============================================================================

function check_rotated_cones(constraint)
    println("Rotated cone constraints that are not tight:")
    for k in 1:L, l in 1:K
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

function f(tau, lambda)
    return (1 - 1/(1 + tau)^lambda)/tau  
end

# =============================================================================
#                                   Models
# =============================================================================

# =============================================================================
# ============================= MISOCP formulation ============================
# =============================================================================

# Implementation of the paper:  
# A mixed-integer quadratically-constrained programming model for 
# the distribution system expansion planning
# They use the Nahman Peric network to test it 

function UL_BFM(sets, costs, substation_param, conductor_param, losses, others, demand)       

    # ========================= Fetch the parameters ========================
    #= The set of parameters is found in the paper:
        "A Constructive Heuristic Algorithm for Distribution System Planning"
    =#

    # Sets
    N, Ns, K, L, Omega_sending, Omega_receiving = sets
    Nu          = setdiff(N, Ns)

    N_size, Ns_size , K_size, L_size = length(N), length(Ns), length(K), length(L)

    lambda = 20 # nb of years in the planning horizon
    K_l    = 0.1
    K_s    = 0.1
    
    loss_factor_circuits = 0.35
    loss_factor_subs     = 0.35
    HOURS_PER_YEAR          = 8760
    C_fs =  1e6 #  [US$] substation fixed cost ($)
    C_vs = 0.1 #[$/(kVAh^2)] substation operation (variable) cost ($/(kWh))
    # c_l: energy cost 
    C_l = 0.05 # [$/kWh]
    S_rating_init = [4 0]
    S_rating_add_max = []
    tau_l = 0.1
    tau_s = 0.1
    cos_phi = 0.9 # power factor of lods 
    b_max = 10


    lambda = 1 # year
    tau_s  = 0.01 # interest rate of a substation
    tau_l  = 0.01 # interest rate of lines

    # ======================== Set up the Gurobi solver =======================
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "TimeLimit", 100)
    set_optimizer_attribute(model, "Presolve", 0)

    # ============================== Variables ================================ 

    @variables( model,   
                begin 
                MIN_VOLTAGE^2 <= V_sqr[1:N_size]  <= MAX_VOLTAGE^2
                I_sqr_k[1:L_size, 1:K_size]       >= 0
                I_sqr[1:L_size]
                P_G[1:N_size]                           >= 0
                Q_G[1:N_size]
                S_G[1:N_size]  
                P_ij_k[1:L_size, 1:K_size] # forward direction
                Q_ij_k[1:L_size, 1:K_size]  # forward direction
                P_ij[1:L_size] # forward direction
                Q_ij[1:L_size]  # forward direction
                b[1:L_size]                     , Bin
                y_send[1:L_size]                , Bin
                y_rec[1:L_size]                 , Bin
                alpha[1:L_size, 1:K_size]           , Bin
                beta[1:Ns_size]                 , Bin      
                end
            )

    for i in Nu
        fix(S_G[i], 0.0; force=true)
        fix(P_G[i], 0.0) 
        fix(Q_G[i], 0.0)
    end

    # ============================= Constraints ===============================

    # CONSTRAINT (6) -> means constraint (6) in the paper
    # The cost of the line already embeds the line length
    @objective(model, Min,  K_l * sum(alpha[l, k] * line_cost[l, k] for l in L, k in K)
                            + HOURS_PER_YEAR * f(tau_l, lambda) * loss_factor_circuits * C_l 
                                * sum(R[l, k] * I_sqr_k[l, k] for l in L, k in K)
                            + K_s * sum(beta[i] * C_fs for i in Ns) 
                            + HOURS_PER_YEAR * f(tau_s, lambda) * loss_factor_subs * C_vs 
                                * sum(S_G[i]^2 for i in Ns)
                )

    
    # P_ij is defined as being the sending-end power
    # this formulation is from the paper active dnep
    # CONSTRAINT (7)
    @constraint(model, 
                active_balance[i=N], 
                P_D[i] - P_G[i]
                == sum(P_ij_k[l, k] - R[l, k] * I_sqr_k[l, k] 
                        for l in Omega_receiving[i], k in K)
                -  sum(P_ij_k[l, k]  
                       for l in Omega_sending[i], k in K)
            )

    # CONSTRAINT (8)
    @constraint(model, 
                reactive_balance[i=N], 
                Q_D[i] - Q_G[i]
                == sum(Q_ij_k[l, k] - X[l, k] * I_sqr_k[l, k] 
                        for l in Omega_receiving[i], k in K)
                -  sum(Q_ij_k[l, k] 
                        for l in Omega_sending[i], k in K)
                )
    
     # CONSTRAINT (9)
     @constraint(model, 
                voltage[l=L],
                V_sqr[line_ends[2]] - V_sqr[line_ends[1]] 
                == sum(-2 * (R[l, k] * P_ij_k[l, k] + X[l, k] * Q_ij_k[l, k]) 
                + (R[l, k]^2 * X[l, k]^2) * I_2_k[l, k] for k in K) 
                + b[l]
            )
    
    # CONSTRAINT (10)
    @constraint(model, 
                rotated_cone[l=L], 
                [V_sqr[line_ends[l][1]] / 2;
                I_sqr[l];  
                P_ij[l];
                Q_ij[l]] 
                in RotatedSecondOrderCone()
            )

    # CONSTRAINT (11)
    @constraint(model,
                line_current[l=L],
                I_sqr[l] == sum(I_sqr_k[l, k] for k in K) 
            )

    # CONSTRAINT (12)
    @constraint(model,
                line_active_power[l=L],
                P_ij[l] == sum(P_ij_k[l, k] for k in K) 
            )
    
    # CONSTRAINT (13)
    @constraint(model,
                line_reactive_power[l=L],
                Q_ij[l] == sum(Q_ij_k[l, k] for k in K) 
            )
    
    # CONSTRAINT (14) is in the definition of the variables


    # CONSTRAINT (15)
    @constraint(model,
                current_limit1[l=L, k=K],
                I_sqr_k[l, k] <= max_current[l, k]^2 * (y_send[l] + y_rec[l])
    )

    # CONSTRAINT (16)
    @constraint(model,
                conductor_choice[l=L],
                sum(alpha[l, k] for k in K) == y_send[l] + y_rec[l]
    )

    # CONSTRAINT (17)
    @constraint(model,
                one_powflow_per_line[l=L],
                y_send[l] + y_rec[l] <= 1
    )
    
    # CONSTRAINT (18)
    @constraint(model,
                radiality,
                sum(y_send[l] + y_rec[l] for l in L) == N_size - Ns_size
                )

    # CONSTRAINT (19)
    @constraint(model,
                loads_are_connected[i=Nu],
                sum(y_rec[l] for l in Omega_sending[i])
                + sum(y_send[l] for l in Omega_receiving[i]) >= 1
    )

    # CONSTRAINT (20)
    @constraint(model,
                b_limit1[l=L],
                b[l] <= b_max * (1- y_send[l] - y_rec[l])
    )

    @constraint(model,
                b_limit2[l=L],
                b[l] >= - b_max * (1 - y_send[l] - y_rec[l])
    )

    # CONSTRAINT (21)
    @constraint(model,
                P_limit1[l=L, k=K],
                P_ij_k[l, k] <= max_current[l, k] * MAX_VOLTAGE * M * y_send[l] 
                )

    # CONSTRAINT (22)
    @constraint(model,
                P_limit2[l=L, k=K],
                P_ij_k[l, k] >= - max_current[l, k] * MAX_VOLTAGE * y_rec[l] 
                )

    # CONSTRAINT (23)
    @constraint(model,
                Q_limit1[l=L, k=K],
                Q_ij_k[l, k] <= max_current[l, k] * MAX_VOLTAGE * (y_send[l] + y_rec[l]) 
                )

    @constraint(model,
                Q_limit2[l=L, k=K],
                Q_ij_k[l, k] >= - max_current[l, k] * MAX_VOLTAGE * (y_send[l] + y_rec[l]) 
                )

    # CONSTRAINT (24)
    @constraint(model,
                current_limit2[l=L, k=K],
                I_sqr_k[l, k] <= max_current[l, k]^2 * alpha[l, k] 
                )
    
    # CONSTRAINT (25)
    @constraint(model,
                P_limit3[l=L, k=K],
                P_ij_k[l, k] <= max_current[l, k] * MAX_VOLTAGE * alpha[l, k] 
                )

    @constraint(model,
                P_limit4[l=L, k=K],
                P_ij_k[l, k] >= - max_current[l, k] * MAX_VOLTAGE * alpha[l, k] 
                )
    
    # CONSTRAINT (26)
    @constraint(model,
                Q_limit3[l=L, k=K],
                Q_ij_k[l, k] <= max_current[l, k] * MAX_VOLTAGE * alpha[l, k] 
                )

    @constraint(model,
                Q_limit4[l=L, k=K],
                Q_ij_k[l, k] >= - max_current[l, k] * MAX_VOLTAGE * alpha[l, k] 
                )

    # CONSTRAINT (27)
    @constraint(model, 
                substation_power[i=Ns], 
                [S_G[i]; P_G[i]; Q_G[i]] in SecondOrderCone()
                )

    # CONSTRAINT (28)
    @constraint(model, 
                substation_capa_limit[i=Ns], 
                S_G[i] <= S_rating_init[i] + beta[i] * S_rating_max[i]
                )


    # Reference voltage of slack node
    @constraint(model, 
                ref_voltage1[i=Ns], 
                V_sqr[i] - 1 <= (MAX_VOLTAGE^2-1)*(1-beta[i])
    )

    @constraint(model, 
                ref_voltage2[i=Ns], 
                V_sqr[i] - 1 >= (MIN_VOLTAGE^2 - 1)*(1-beta[i])
    )
    
    
    

    
   
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