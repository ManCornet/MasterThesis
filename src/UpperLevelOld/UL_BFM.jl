#-----------------------------------------------------------------------------
#
#                   - TFE : Upper level problem formulation - 
#                             University of Liege
#
#-----------------------------------------------------------------------------
# Created By  : Manon Cornet
# Created Date: Friday April 10th 2023
#
# UL_BFM:
#   File containing the BFM DNEP model
#   Paper ref:
#   "Franco, J. F., Rider, M. J., & Romero, R. (2014). A mixed-integer 
#   quadratically-constrained programming model for the distribution system 
#   expansion planning. International Journal of Electrical Power & Energy 
#   Systems, 62, 265-272."
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

function UL_BFM_1P(network_dict::Dict, obj_dict::Dict; radiality::Integer=0, generation::Bool=false)       

    # ========================= Network parameters ========================

    N, Omega_sending, Omega_receiving, MIN_VOLTAGE, MAX_VOLTAGE = network_dict[:bus]
    Ns, S_rating_init, S_rating_max = network_dict[:sub_bus]
    Nu, P_D, Q_D, delta_t = network_dict[:load_bus]
    L, line_ends, max_i, R, X, G, B = network_dict[:line]
    K = network_dict[:conductor]
    M = MAX_VOLTAGE^2
    Ns_init = [1] 
    Ns_notinit = setdiff(Ns, Ns_init)

    # ===================== Obj. Function parameters =======================

    loss_cost, sub_install_cost, sub_op_cost, line_cost = obj_dict[:costs]
    line_loss_factor, sub_loss_factor = obj_dict[:LF]
    K_l, K_s = obj_dict[:CRF]
    f_l, f_s = obj_dict[:NPV_coeff]

    HOURS_PER_YEAR = 8760

    # ======================== Set up the Gurobi solver =======================
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "TimeLimit", 200)
    set_optimizer_attribute(model, "Presolve", 0)

    # ============================== Variables ================================ 

    @variables( model,   
                begin 
                MIN_VOLTAGE^2 <= V_sqr[N] <= MAX_VOLTAGE^2 , (container=Array)
                I_sqr_k[L, K] >= 0                         , (container=Array)
                I_sqr[L]      >= 0                         , (container=Array)
                P_G[N]                                     , (container=Array)
                Q_G[N]                                     , (container=Array)
                S_G[N] >= 0                                   , (container=Array)
                P_ij_k[L, K]                               , (container=Array) 
                Q_ij_k[L, K]                               , (container=Array)
                P_ij[L]                                    , (container=Array)
                Q_ij[L]                                    , (container=Array)
                y_send[L]                                  , (container=Array, binary=true)
                y_rec[L]                                   , (container=Array, binary=true)
                alpha[L, K]                                , (container=Array, binary=true)
                beta[Ns]                                   , (container=Array, binary=true)    
                end
            )

    # Test the 3 radiality constraints from Jabr Paper
    if radiality == 1 || radiality == 2
        # SINGLE COMMODITY FLOW CONSTRAINT
        @variable(model, k_ij[L], container=Array)
        radiality == 1 && @variable(model, K_i[N], container=Array)

    elseif radiality == 3
        @variable(model, k_ij[L, N], integer=true, container=Array)

    elseif radiality == 4
        @variable(model, z_ij[L, N], binary=true, container=Array)
        @variable(model, z_ji[L, N], binary=true, container=Array)
    elseif radiality == 5
        @variable(model, k_ij[L, N], integer= true, container=Array)
        @variable(model, x_ij[L, N], binary = true, container=Array)
        @variable(model, x_span_ij[L, N], binary = true, container=Array)
    end



    for i in Nu
        fix(P_G[i], 0.0) 
        fix(Q_G[i], 0.0)
        fix(S_G[i], 0.0; force=true)
    end

    # ============================= Constraints ===============================

    # CONSTRAINT (6) -> means constraint (6) in the paper
    @objective(model, Min,    K_l            * sum(alpha[l, k] * line_cost[l, k] for l in L, k in K)
                            + K_s            * sum(beta[i] * sub_install_cost for i in Ns) 
                            + HOURS_PER_YEAR * f_l * line_loss_factor * loss_cost * delta_t
                                             * sum(R[l, k] * I_sqr_k[l, k] * BASE_POWER for l in L, k in K)
                            + HOURS_PER_YEAR * f_s * sub_loss_factor * sub_op_cost 
                                             * sum((P_G[i]^2 + Q_G[i]^2) * (delta_t * BASE_POWER)^2 for i in Ns)
                )

    
    # P_ij is defined as being the sending-end power
    # this formulation is from the paper active dnep
    # CONSTRAINT (7)
    @constraint(model, 
                active_balance[i=N], 
                P_D[i] - P_G[i]
                == sum(P_ij_k[l, k] -  R[l, k] * I_sqr_k[l, k] 
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
                voltage_value1[l=L],
                V_sqr[line_ends[l][2]] - V_sqr[line_ends[l][1]] 
                <= sum(-2 * (R[l, k] * P_ij_k[l, k] + X[l, k] * Q_ij_k[l, k]) 
                + (R[l, k]^2 + X[l, k]^2) * I_sqr_k[l, k] for k in K) 
                + M * (1 - (y_send[l] + y_rec[l]))
            )
    @constraint(model, 
                voltage_value2[l=L],
                V_sqr[line_ends[l][2]] - V_sqr[line_ends[l][1]] 
                >= sum(-2 * (R[l, k] * P_ij_k[l, k] + X[l, k] * Q_ij_k[l, k]) 
                + (R[l, k]^2 + X[l, k]^2) * I_sqr_k[l, k] for k in K) 
                - M * (1 - (y_send[l] + y_rec[l]))
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

    # Here also look if we can gain some time by remowing those variables since two times more variable
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
                I_sqr_k[l, k] <= max_i[l, k]^2 * (y_send[l] + y_rec[l])
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
    
    



    # CONSTRAINT (19)
    @constraint(model,
                loads_are_connected[i=Nu],
                sum(y_rec[l] for l in Omega_sending[i])
                + sum(y_send[l] for l in Omega_receiving[i]) >= 1
    )
    
    #=
    # CONSTRAINT (20)
    @constraint(model,
                b_limit1[l=L],
                b[l] <= B_MAX * (1 - y_send[l] - y_rec[l])
    )

    @constraint(model,
                b_limit2[l=L],
                b[l] >= - B_MAX * (1 - y_send[l] - y_rec[l])
    )=#

    # CONSTRAINT (21)
    @constraint(model,
                P_limit1[l=L, k=K],
                P_ij_k[l, k] <= max_i[l, k] * MAX_VOLTAGE * y_send[l] 
                )

    # CONSTRAINT (22)
    @constraint(model,
                P_limit2[l=L, k=K],
                P_ij_k[l, k] >= (- max_i[l, k] * MAX_VOLTAGE) * y_rec[l] 
                )

    # Why not a direction for the reactive power ?
    # CONSTRAINT (23)
    @constraint(model,
                Q_limit1[l=L, k=K],
                Q_ij_k[l, k] <= max_i[l, k] * MAX_VOLTAGE * (y_send[l] + y_rec[l]) 
                )

    @constraint(model,
                Q_limit2[l=L, k=K],
                Q_ij_k[l, k] >= (- max_i[l, k] * MAX_VOLTAGE) * (y_send[l] + y_rec[l]) 
                )

    # CONSTRAINT (24)
    @constraint(model,
                current_limit2[l=L, k=K],
                I_sqr_k[l, k] <= max_i[l, k]^2 * alpha[l, k] 
                )
    
    # CONSTRAINT (25)
    @constraint(model,
                P_limit3[l=L, k=K],
                P_ij_k[l, k] <= max_i[l, k] * MAX_VOLTAGE * alpha[l, k] 
                )

    @constraint(model,
                P_limit4[l=L, k=K],
                P_ij_k[l, k] >= - max_i[l, k] * MAX_VOLTAGE * alpha[l, k] 
                )
    
    # CONSTRAINT (26)
    @constraint(model,
                Q_limit3[l=L, k=K],
                Q_ij_k[l, k] <= max_i[l, k] * MAX_VOLTAGE * alpha[l, k] 
                )

    @constraint(model,
                Q_limit4[l=L, k=K],
                Q_ij_k[l, k] >= - max_i[l, k] * MAX_VOLTAGE * alpha[l, k] 
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

    

    # ========================= Additional Constraints ========================

    if generation
        # CONSTRAINT : ADD GENERATION EVERYWHERE

        @constraint(model, 
                    generation[i=5:5], 
                    P_G[i] == 1.5
                )
        @constraint(model, 
                    generation2[i=Nu], 
                    P_G[i] <= 10
                )
    end
   
    # LOAD OVER-SATISFACTION
    @constraint(model, 
                power_balance1, 
                sum(P_G[i] - P_D[i] for i in N) >= 0
    )

    @constraint(model, 
                power_balance2, 
                sum(Q_G[i] - Q_D[i] for i in N) >= 0
    )


    # REFERENCE VOLTAGE AT SLACK BUS
    @constraint(model, 
                ref_voltage_sub_init[i=Ns_init], 
                V_sqr[i] == 1
    )

 
    # REFERENCE VOLTAGE AT SLACK BUS
    @constraint(model, 
                ref_voltage_sub_notinit1[i=Ns_notinit], 
                V_sqr[i] - 1 <= (MAX_VOLTAGE^2-1)*(1-beta[i])
    )

    @constraint(model, 
                ref_voltage_sub_notinit2[i=Ns_notinit], 
                V_sqr[i] - 1 >= (MIN_VOLTAGE^2 - 1)*(1-beta[i])
    ) 

    # CONSTRAINT (18)
    @constraint(model,
            radiality_constraint,
            sum(y_send[l] + y_rec[l] for l in L) == length(N) - length(Ns)
    )

    

    if radiality == 1
        # ATTENTION IL FAUT FORCER le flux à valoir 0 lorsque la sous-station n'est pas construite
        # ------------------------------------------
        # - SINGLE-COMMODITY FLOW CONSTRAINTS (V1) -
        # ------------------------------------------ 

        # CONSTRAINT (18)
        # @constraint(model,
        #             radiality_constraint,
        #             sum(y_send[l] + y_rec[l] for l in L) == length(N) - length(Ns)
        # )


        @constraint(model, single_comm_constr1[i=Ns_notinit], K_i[i] <= beta[i]*length(Nu))
        @constraint(model, single_comm_constr2[i=Ns_notinit], K_i[i] >= 0)

        @constraint(model, single_comm_constr1bis[i=Ns_init], K_i[i] <= length(Nu))
        @constraint(model, single_comm_constr2bis[i=Ns_init], K_i[i] >= 0)


        @constraint(model, single_comm_constr3[i=Nu], K_i[i] == -1)

        @constraint(model,
                    single_comm_constr4[i=N],
                    K_i[i] == - sum(k_ij[l] for l in Omega_receiving[i])
                            + sum(k_ij[l] for l in Omega_sending[i])
                )

        @constraint(model,
                    single_comm_constr6[l=L],
                    k_ij[l] <= (length(Nu)) * (y_send[l] + y_rec[l])
        )

        @constraint(model,
                    single_comm_constr7[l=L],
                    k_ij[l] >= - (length(Nu)) * (y_send[l] + y_rec[l])
        )

    elseif radiality == 2
        # ------------------------------------------
        # - SINGLE-COMMODITY FLOW CONSTRAINTS (V2) -
        # ------------------------------------------ 
        # CONSTRAINT (18)
        @constraint(model,
                    radiality_constraint,
                    sum(y_send[l] + y_rec[l] for l in L) == length(N) - length(Ns))
        
        @constraint(model, single_comm_constr1[i=Nu], 
                         -sum(k_ij[l] for l in Omega_receiving[i])
                        + sum(k_ij[l] for l in Omega_sending[i]) == -1)

        @constraint(model, single_comm_constr2[i=Ns], 
                        - sum(k_ij[l] for l in Omega_receiving[i])
                       + sum(k_ij[l] for l in Omega_sending[i]) >= 0)

        @constraint(model, single_comm_constr3[i=Ns_init], 
                       - sum(k_ij[l] for l in Omega_receiving[i])
                      + sum(k_ij[l] for l in Omega_sending[i]) <= length(Nu))

        
        @constraint(model, single_comm_constr5[i=Ns_notinit], 
                    - sum(k_ij[l] for l in Omega_receiving[i])
                    + sum(k_ij[l] for l in Omega_sending[i]) <= length(Nu) * beta[i])

        @constraint(model,
                    single_comm_constr6[l=L],
                    k_ij[l] <= length(Nu) * (y_send[l] + y_rec[l])
        )

        @constraint(model,
                    single_comm_constr7[l=L],
                    k_ij[l] >= - length(Nu) * (y_send[l] + y_rec[l])
        )

    elseif radiality == 3
        # ------------------------------------------
        # ---- MULTI-COMMODITY FLOW CONSTRAINTS ----
        # ------------------------------------------
         # CONSTRAINT (18)
         @constraint(model,
            radiality_constraint,
            sum(y_send[l] + y_rec[l] for l in L) == length(N) - length(Ns))
        
        @constraint(model, multi_comm_constr1[i=Ns, w=Nu], 
                    - sum(k_ij[l, w] for l in Omega_receiving[i])
                    + sum(k_ij[l, w] for l in Omega_sending[i]) >= 0)
        
        @constraint(model, multi_comm_constr2[i=Ns_init, w=Nu], 
                    - sum(k_ij[l, w] for l in Omega_receiving[i])
                    + sum(k_ij[l, w] for l in Omega_sending[i]) <= 1)
   
        @constraint(model, 
                    multi_comm_constr3[i=Ns_notinit, w=Nu],
                    - sum(k_ij[l, w] for l in Omega_receiving[i])
                    + sum(k_ij[l, w] for l in Omega_sending[i]) <= beta[i])

        @constraint(model, multi_comm_constr4[i=Nu], 
                    - sum(k_ij[l, i] for l in Omega_receiving[i])
                    + sum(k_ij[l, i] for l in Omega_sending[i]) == -1)

        for i in Nu, w in Nu
            if i !== w 
                @constraint(model, - sum(k_ij[l, w] for l in Omega_receiving[i])
                                    + sum(k_ij[l, w] for l in Omega_sending[i]) == 0)
            
            end
        end


        @constraint(model,
            multi_comm_constr5[l=L, w=Nu],
            k_ij[l, w] <= (y_send[l] + y_rec[l])
        )

        @constraint(model,
            multi_comm_constr6[l=L, w=Nu],
            k_ij[l, w] >= - (y_send[l] + y_rec[l])
        )
    elseif radiality == 4
        # ------------------------------------------
        # ------- SPANNING TREE  CONSTRAINTS -------
        # ------------------------------------------
        # CONSTRAINT (18)
        @constraint(model,
                    radiality_constraint,
                    sum(y_send[l] + y_rec[l] for l in L) == length(N) - length(Ns))
        
        @constraint(model,
                    span_tree1[l=L, w=N],
                    (y_send[l] + y_rec[l]) == z_ij[l, w] + z_ji[l, w] 
        )

        for i in N, w in N
            if i <= length(Ns_init)
                if i !== w 
                    @constraint(model, + sum(z_ji[l, w] for l in Omega_receiving[i])
                                       + sum(z_ij[l, w] for l in Omega_sending[i]) >= 0
                    )
                    @constraint(model, + sum(z_ji[l, w] for l in Omega_receiving[i])
                                       + sum(z_ij[l, w] for l in Omega_sending[i]) <= 1
                    )
                end
            elseif i > length(Ns_init) && i <= length(Ns)
                if i !== w 
                    @constraint(model,  + sum(z_ji[l, w] for l in Omega_receiving[i])
                                        + sum(z_ij[l, w] for l in Omega_sending[i]) >= 0
                    )
                    @constraint(model,  + sum(z_ji[l, w] for l in Omega_receiving[i])
                                        + sum(z_ij[l, w] for l in Omega_sending[i]) <= beta[i]
                    )
                end
            elseif i !== w 
                @constraint(model, + sum(z_ji[l, w] for l in Omega_receiving[i])
                                    + sum(z_ij[l, w] for l in Omega_sending[i]) == 1
                )
            end    
        end

        @constraint(model, 
                    span_tree2[w=N],
                    + sum(z_ji[l, w] for l in Omega_receiving[w])
                    + sum(z_ij[l, w] for l in Omega_sending[w]) == 0
    )
    elseif radiality == 5

        @constraint(model,
                    radiality_constraint,
                    sum(y_send[l] + y_rec[l] for l in L) == length(N) - length(Ns))
        
        #@constraint(model, 
        #            mcf_v2_1[l=L],
        #            y_send[l] + y_rec[l] == x_span[l]
        #)

        #@constraint(model, 
        #            mcf_v2_2[l=L],
        #            x_span[l] <= x[l]
        #)

        @constraint(model, mcf_v2_3[i=Ns, w=Nu], 
                    - sum(k_ij[l, w] for l in Omega_receiving[i])
                    + sum(k_ij[l, w] for l in Omega_sending[i]) >= 0)
        
        @constraint(model, mcf_v2_4[i=Ns_init, w=Nu], 
                    - sum(k_ij[l, w] for l in Omega_receiving[i])
                    + sum(k_ij[l, w] for l in Omega_sending[i]) <= 1)
   
        @constraint(model, 
                    mcf_v2_5[i=Ns_notinit, w=Nu],
                    - sum(k_ij[l, w] for l in Omega_receiving[i])
                    + sum(k_ij[l, w] for l in Omega_sending[i]) <= beta[i])

        @constraint(model, mcf_v2_6[i=Nu], 
                    - sum(k_ij[l, i] for l in Omega_receiving[i])
                    + sum(k_ij[l, i] for l in Omega_sending[i]) == -1)

        for i in Nu, w in Nu
            if i !== w 
                @constraint(model, - sum(k_ij[l, w] for l in Omega_receiving[i])
                                    + sum(k_ij[l, w] for l in Omega_sending[i]) == 0)
            
            end
        end

        @constraint(model, mcf_v2_7[l=L, w=Nu], 
                    k_ij[l, w] <= y_send[l])

        @constraint(model, mcf_v2_8[l=L, w=Nu], 
                    k_ij[l, w] >= - y_rec[l])
   
        

    end
    # REFERENCE VOLTAGE AT SLACK BUS
   
    #print(model)

    optimize!(model)

    solution_summary(model, verbose=false)

    if termination_status(model) == MOI.OPTIMAL
        var_values = Dict(  string(k) => value.(v) for (k, v) 
                            in object_dictionary(model) 
                            if (v isa AbstractArray{VariableRef} || 
                                v isa VariableRef)
                        )
        
        var_sets = Dict("V_sqr"     => ["node"], 
                        "I_sqr_k"   => ["line", "conductor"],
                        "I_sqr"     => ["line"],
                        "P_G"       => ["node"],
                        "Q_G"       => ["node"],
                        "S_G"       => ["node"],
                        "P_ij_k"    => ["line", "conductor"],
                        "Q_ij_k"    => ["line", "conductor"],
                        "P_ij"      => ["line"],
                        "Q_ij"      => ["line"],
                        "y_send"    => ["line"],
                        "y_rec"     => ["line"],
                        "alpha"     => ["line", "conductor"],
                        "beta"      => ["node"]
                        )

        if radiality == 1 || radiality == 2
            var_sets["k_ij"] = ["line"]
            radiality == 1 && (var_sets["K_i"] = ["node"])

        elseif radiality == 3
            var_sets["k_ij"] = ["line", "node"]
        elseif radiality == 4
            var_sets["z_ij"] = ["line", "node"]
            var_sets["z_ji"] = ["line", "node"]
        end

        return objective_value(model), var_values, var_sets, solve_time(model) 

    elseif termination_status(model) == DUAL_INFEASIBLE
        println("problem unbounded")

    elseif termination_status(model) == MOI.INFEASIBLE
        println("problem infeasible")
    end
end



function UL_BFM_2P(network_dict::Dict, obj_dict::Dict)       

    # ========================= Network parameters ========================

    N, Omega_sending, Omega_receiving, MIN_VOLTAGE, MAX_VOLTAGE = network_dict[:bus]
    Ns, S_rating_init, S_rating_max = network_dict[:sub_bus]
    Nu, P_D, Q_D, delta_t = network_dict[:load_bus]
    L, line_ends, max_i, R, X, G, B = network_dict[:line]
    K = network_dict[:conductor]
    M = MAX_VOLTAGE^2
    Ns_init = [1] 
    Ns_notinit = setdiff(Ns, Ns_init)

    # ===================== Obj. Function parameters =======================

    loss_cost, sub_install_cost, sub_op_cost, line_cost = obj_dict[:costs]
    line_loss_factor, sub_loss_factor = obj_dict[:LF]
    K_l, K_s = obj_dict[:CRF]
    f_l, f_s = obj_dict[:NPV_coeff]

    HOURS_PER_YEAR = 8760

    # ======================== Set up the Gurobi solver =======================
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "TimeLimit", 200)
    set_optimizer_attribute(model, "Presolve", 0)

    # ============================== Variables ================================ 

    @variables( model,   
                begin 
                MIN_VOLTAGE^2 <= V_sqr[N]  <= MAX_VOLTAGE^2, (container=Array)
                I_sqr_k[L, K] >= 0                         , (container=Array)
                I_sqr[L]      >= 0                         , (container=Array)
                P_G[N] >= 0                                , (container=Array)
                Q_G[N]                                     , (container=Array)
                S_G[N]                                     , (container=Array)
                P_ij_k[L, K]                               , (container=Array) 
                Q_ij_k[L, K]                               , (container=Array)
                P_ji_k[L, K]                               , (container=Array) 
                Q_ji_k[L, K]                               , (container=Array)
                P_ij[L]                                    , (container=Array)
                Q_ij[L]                                    , (container=Array)
                P_ji[L]                                    , (container=Array)
                Q_ji[L]                                    , (container=Array)
                b[L]                                       , (container=Array, binary=true)
                x[L]                                       , (container=Array, binary=true)
                alpha[L, K]                                , (container=Array, binary=true)
                beta[Ns]                                   , (container=Array, binary=true)    
                end
            )

    #for i in Nu
        #fix(P_G[i], 0.0; force=true) 
        #fix(Q_G[i], 0.0)
        #fix(S_G[i], 0.0)
    #end

    # ============================= Constraints ===============================

    # CONSTRAINT (6) -> means constraint (6) in the paper
    @objective(model, Min,    K_l            * sum(alpha[l, k] * line_cost[l, k] for l in L, k in K)
                            + K_s            * sum(beta[i] * sub_install_cost for i in Ns) 
                            + HOURS_PER_YEAR * f_l * line_loss_factor * loss_cost * delta_t
                                             * sum(R[l, k] * I_sqr_k[l, k] * BASE_POWER for l in L, k in K)
                            + HOURS_PER_YEAR * f_s * sub_loss_factor * sub_op_cost 
                                             * sum((P_G[i]^2 + Q_G[i]^2) * (delta_t * BASE_POWER)^2 for i in Ns)
                )

    
    # P_ij is defined as being the sending-end power
    # this formulation is from the paper active dnep
    # CONSTRAINT (7)
    @constraint(model, 
                active_balance[i=N], 
                P_D[i] - P_G[i]
                == sum(P_ji_k[l, k] 
                        for l in Omega_receiving[i], k in K)
                -  sum(P_ij_k[l, k] 
                       for l in Omega_sending[i], k in K)
            )

    @constraint(model,
                receiving_active_power[l=L, k=K],
                P_ji_k[l, k] == P_ij_k[l, k] -  R[l, k] * I_sqr_k[l, k]
    )
    # CONSTRAINT (8)
    @constraint(model, 
                reactive_balance[i=N], 
                Q_D[i] - Q_G[i]
                == sum(Q_ji_k[l, k]
                        for l in Omega_receiving[i], k in K)
                -  sum(Q_ij_k[l, k]
                        for l in Omega_sending[i], k in K)
                )

    @constraint(model,
                receiving_reactive_power[l=L, k=K],
                Q_ji_k[l, k] == Q_ij_k[l, k] -  X[l, k] * I_sqr_k[l, k]
    )
    
     # CONSTRAINT (9)
     @constraint(model, 
                voltage_value1[l=L],
                V_sqr[line_ends[l][2]] - V_sqr[line_ends[l][1]] 
                <= - sum(2 * (R[l, k] * P_ij_k[l, k] + X[l, k] * Q_ij_k[l, k]) 
                + (R[l, k]^2 + X[l, k]^2) * I_sqr_k[l, k] for k in K) 
                + M * (1 - x[l])
            )
    @constraint(model, 
                voltage_value2[l=L],
                V_sqr[line_ends[l][2]] - V_sqr[line_ends[l][1]] 
                >= - sum(2 * (R[l, k] * P_ij_k[l, k] + X[l, k] * Q_ij_k[l, k]) 
                + (R[l, k]^2 + X[l, k]^2) * I_sqr_k[l, k] for k in K) 
                - M * (1 - x[l])
        )
    
    # CONSTRAINT (10)
    @constraint(model, 
                rotated_cone_s[l=L], 
                [V_sqr[line_ends[l][1]] / 2;
                I_sqr[l];  
                P_ij[l];
                Q_ij[l]] 
                in RotatedSecondOrderCone()
            )
    @constraint(model, 
                rotated_cone_r[l=L], 
                [V_sqr[line_ends[l][2]] / 2;
                I_sqr[l];  
                P_ji[l];
                Q_ji[l]] 
                in RotatedSecondOrderCone()
        )

    # CONSTRAINT (11)
    @constraint(model,
                line_current[l=L],
                I_sqr[l] == sum(I_sqr_k[l, k] for k in K) 
            )

    # CONSTRAINT (12)
    @constraint(model,
                line_active_power_s[l=L],
                P_ij[l] == sum(P_ij_k[l, k] for k in K) 
            )
    @constraint(model,
                line_active_power_r[l=L],
                P_ji[l] == sum(P_ji_k[l, k] for k in K) 
                )
    
    # CONSTRAINT (13)
    @constraint(model,
                line_reactive_power_s[l=L],
                Q_ij[l] == sum(Q_ij_k[l, k] for k in K) 
            )

    @constraint(model,
                line_reactive_power_r[l=L],
                Q_ji[l] == sum(Q_ji_k[l, k] for k in K) 
    )
    
    # CONSTRAINT (14) is in the definition of the variables


    # CONSTRAINT (15)
    @constraint(model,
                current_limit1[l=L, k=K],
                I_sqr_k[l, k] <= max_i[l, k]^2 * x[l]
    )

    # CONSTRAINT (16)
    @constraint(model,
                conductor_choice[l=L],
                sum(alpha[l, k] for k in K) == x[l]
    )

    
    # CONSTRAINT (18)
    @constraint(model,
                radiality_constraint,
                sum(x[l] for l in L) == length(N) - length(Ns)
                )

    
    # CONSTRAINT (19)
    @constraint(model,
                loads_are_connected[i=Nu],
                sum(x[l] for l in Omega_sending[i])
                + sum(x[l] for l in Omega_receiving[i]) >= 1
    )
    
    
    #=
    # CONSTRAINT (20)
    @constraint(model,
                b_limit1[l=L],
                b[l] <= B_MAX * (1 - y_send[l] - y_rec[l])
    )

    @constraint(model,
                b_limit2[l=L],
                b[l] >= - B_MAX * (1 - y_send[l] - y_rec[l])
    )=#

    # CONSTRAINT (21)
    @constraint(model,
                P_limit1_s[l=L, k=K],
                P_ij_k[l, k] <= max_i[l, k] * MAX_VOLTAGE * x[l] 
                )

    @constraint(model,
                P_limit1_r[l=L, k=K],
                P_ji_k[l, k] <= max_i[l, k] * MAX_VOLTAGE * x[l] 
                )

    # CONSTRAINT (22)
    @constraint(model,
                P_limit2_s[l=L, k=K],
                P_ij_k[l, k] >= - max_i[l, k] * MAX_VOLTAGE * x[l] 
                )

    @constraint(model,
                P_limit2_r[l=L, k=K],
                P_ji_k[l, k] >= - max_i[l, k] * MAX_VOLTAGE * x[l] 
                )


    # CONSTRAINT (23)
    @constraint(model,
                Q_limit1_s[l=L, k=K],
                Q_ij_k[l, k] <= max_i[l, k] * MAX_VOLTAGE * x[l]  
                )

    @constraint(model,
                Q_limit2_s[l=L, k=K],
                Q_ij_k[l, k] >= - max_i[l, k] * MAX_VOLTAGE * x[l]  
                )

    @constraint(model,
                Q_limit1_r[l=L, k=K],
                Q_ji_k[l, k] <= max_i[l, k] * MAX_VOLTAGE * x[l]  
                )

    @constraint(model,
                Q_limit2_r[l=L, k=K],
                Q_ji_k[l, k] >= - max_i[l, k] * MAX_VOLTAGE * x[l]  
                )

    # CONSTRAINT (24)
    @constraint(model,
                current_limit2[l=L, k=K],
                I_sqr_k[l, k] <= max_i[l, k]^2 * alpha[l, k] 
                )
    
    # CONSTRAINT (25)
    @constraint(model,
                P_limit3_s[l=L, k=K],
                P_ij_k[l, k] <= max_i[l, k] * MAX_VOLTAGE * alpha[l, k] 
                )

    @constraint(model,
                P_limit4_s[l=L, k=K],
                P_ij_k[l, k] >= - max_i[l, k] * MAX_VOLTAGE * alpha[l, k] 
                )

    @constraint(model,
                P_limit3_r[l=L, k=K],
                P_ji_k[l, k] <= max_i[l, k] * MAX_VOLTAGE * alpha[l, k] 
                )

    @constraint(model,
                P_limit4_r[l=L, k=K],
                P_ji_k[l, k] >= - max_i[l, k] * MAX_VOLTAGE * alpha[l, k] 
                )
    
    # CONSTRAINT (26)
    @constraint(model,
                Q_limit3_s[l=L, k=K],
                Q_ij_k[l, k] <= max_i[l, k] * MAX_VOLTAGE * alpha[l, k] 
                )

    @constraint(model,
                Q_limit4_s[l=L, k=K],
                Q_ij_k[l, k] >= - max_i[l, k] * MAX_VOLTAGE * alpha[l, k] 
                )

    @constraint(model,
                Q_limit3_r[l=L, k=K],
                Q_ji_k[l, k] <= max_i[l, k] * MAX_VOLTAGE * alpha[l, k] 
                )

    @constraint(model,
                Q_limit4_r[l=L, k=K],
                Q_ji_k[l, k] >= - max_i[l, k] * MAX_VOLTAGE * alpha[l, k] 
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


    # ========================= Additional Constraints ========================

    # CONSTRAINT : ADD GENERATION EVERYWHERE
    @constraint(model, 
                generation[i=Nu], 
                P_G[i] <= 0.1
    )
   
    # LOAD OVER-SATISFACTION
    @constraint(model, 
                power_balance1, 
                sum(P_G[i] - P_D[i] for i in N) >= 0
    )

    @constraint(model, 
                power_balance2, 
                sum(Q_G[i] - Q_D[i] for i in N) >= 0
    )

    # REFERENCE VOLTAGE AT SLACK BUS
    @constraint(model, 
                ref_voltage_sub_init[i=Ns_init], 
                V_sqr[i] == 1
    )

 
    # REFERENCE VOLTAGE AT SLACK BUS
    @constraint(model, 
                ref_voltage_sub_notinit1[i=Ns_notinit], 
                V_sqr[i] - 1 <= (MAX_VOLTAGE^2 - 1)*(1 - beta[i])
    )

    @constraint(model, 
                ref_voltage_sub_notinit2[i=Ns_notinit], 
                V_sqr[i] - 1 >= (MIN_VOLTAGE^2 - 1)*(1 - beta[i])
    ) 
   
    #print(model)

    optimize!(model)

    solution_summary(model, verbose=false)

    if termination_status(model) == MOI.OPTIMAL
        var_values = Dict(  string(k) => value.(v) for (k, v) 
                            in object_dictionary(model) 
                            if (v isa AbstractArray{VariableRef} || 
                                v isa VariableRef)
                        )
        
        var_sets = Dict("V_sqr"     => ["node"], 
                        "I_sqr_k"   => ["line", "conductor"],
                        "I_sqr"     => ["line"],
                        "P_G"       => ["node"],
                        "Q_G"       => ["node"],
                        "S_G"       => ["node"],
                        "P_ij_k"    => ["line", "conductor"],
                        "Q_ij_k"    => ["line", "conductor"],
                        "P_ji_k"    => ["line", "conductor"],
                        "Q_ji_k"    => ["line", "conductor"],
                        "P_ij"      => ["line"],
                        "Q_ij"      => ["line"],
                        "P_ji"      => ["line"],
                        "Q_ji"      => ["line"],
                        "b"         => ["line"],
                        "x"         => ["line"],
                        "alpha"     => ["line", "conductor"],
                        "beta"      => ["node"]
                        )

        return objective_value(model), var_values, var_sets, solve_time(model) 

    elseif termination_status(model) == DUAL_INFEASIBLE
        println("problem unbounded")

    elseif termination_status(model) == MOI.INFEASIBLE
        println("problem infeasible")
    end
end