#-----------------------------------------------------------------------------
#
#                   - TFE : Upper level problem formulation - 
#                             University of Liege
#
#-----------------------------------------------------------------------------
# Created By  : Manon Cornet
# Created Date: Tuesday March 14 2023
#
# UL_Jabr:
#   File containing the Jabr DNEP model 
#   Paper ref: 
#   "Jabr, R. A. (2012). Polyhedral formulations and loop elimination constraints 
#   for distribution network expansion planning. IEEE Transactions on Power 
#   Systems, 28(2), 1888-1897."
# =============================================================================
#                                   Imports
# =============================================================================

using JuMP, Gurobi


# =============================================================================
#                                   Model
# =============================================================================

function UL_Jabr(network_dict::Dict, obj_dict::Dict)       

    # ========================= Network parameters ========================

    N, Omega_sending, Omega_receiving, MIN_VOLTAGE, MAX_VOLTAGE, _ , _ = network_dict[:bus]
    Ns, S_rating_init, S_rating_max = network_dict[:sub_bus]
    Nu, P_D, Q_D, delta_t = network_dict[:load_bus]
    L, line_ends, max_i, R, X, G, B = network_dict[:line]
    K = network_dict[:conductor]
    Ns_init = [1] 
    Ns_notinit = setdiff(Ns, Ns_init)
    
    # ===================== Obj. Function parameters =======================

    loss_cost, sub_install_cost, sub_op_cost, line_cost = obj_dict[:costs]
    line_loss_factor, sub_loss_factor = obj_dict[:LF]
    K_l, K_s = obj_dict[:CRF]
    tau_l, tau_s = obj_dict[:tau]

    HOURS_PER_YEAR = 8760

    # ======================== Set up the Gurobi solver =======================

    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "TimeLimit", 200)
    set_optimizer_attribute(model, "MIPGap", 1e-6)
    set_optimizer_attribute(model, "Presolve", 0)

    # ============================== Variables ================================ 

    @variables( model,   
                begin 
                MIN_VOLTAGE^2 <= V_sqr[N] <= MAX_VOLTAGE^2, (container=Array)
                I_sqr[L, K]                               , (container=Array)
                P_ij[L, K]                                , (container=Array)
                P_ji[L, K]                                , (container=Array)
                Q_ij[L, K]                                , (container=Array)
                Q_ji[L, K]                                , (container=Array)
                X_i_ij[L, K, N]                           , (container=Array)
                X_ij_re[L, K] >= 0                        , (container=Array)
                X_ij_im[L, K]                             , (container=Array)
                P_G[N]                                     , (container=Array)
                Q_G[N]                                    , (container=Array)
                S_G[N] >= 0                               , (container=Array)                          
                alpha[L, K]                               , (container=Array, binary=true)
                beta[Ns]                                  , (container=Array, binary=true)
                x[L]                                      , (container=Array, binary=true)
                end
            )

    for i in Nu
        fix(P_G[i], 0.0) 
        fix(Q_G[i], 0.0)
        fix(S_G[i], 0.0; force=true)
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

    @objective(model, Min,    K_l             * sum(alpha[l, k] * line_cost[l, k] for l in L, k in K) 
                            + K_s             * sum(beta[i] * sub_install_cost for i in Ns)
                            + HOURS_PER_YEAR * (1 + tau_l) 
                                              * line_loss_factor * loss_cost * delta_t
                                              * sum((P_G[i] - P_D[i]) * BASE_POWER for i in N) 
                            + HOURS_PER_YEAR * (1 + tau_s) 
                                              * sub_loss_factor * sum(sub_op_cost
                                              * (S_G[i]^2) * (delta_t * BASE_POWER)^2 for i in Ns)
                    )

    # CONSTRAINT (2) -> In definition of alpha

    # CONSTRAINT (3)
    @constraint(model, line_constructed[l=L], x[l] == sum(alpha[l, k] for k in K))

    # CONSTRAINT (4)
    @constraint(model, 
                substation_capacity[i=Ns], 
                S_G[i] <= S_rating_init[i] + beta[i] * S_rating_max[i]
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
                I_sqr[l, k] <= alpha[l, k] * max_i[l, k]^2
                )

    for l in L, k in K
        ifrom = line_ends[l][1]
        ito = line_ends[l][2]

        # CONSTRAINT (21)
        @constraint(model, P_ij[l, k] 
                           == G[l, k] * (X_i_ij[l, k, ifrom] 
                                                   - X_ij_re[l, k])
                           + B[l, k] * X_ij_im[l, k]
        )

        
        @constraint(model, P_ji[l, k] 
                           == G[l, k] * (X_i_ij[l, k, ito] 
                                                   - X_ij_re[l, k])
                            - B[l, k] * X_ij_im[l, k]) 
                            
      
        # CONSTRAINT (22)
        @constraint(model, Q_ij[l, k] 
                            == B[l, k] * (X_i_ij[l, k, ifrom] 
                                                     - X_ij_re[l, k])
                            - G[l, k] * X_ij_im[l, k])

        
        @constraint(model, Q_ji[l, k] 
                            ==  B[l, k] * (X_i_ij[l, k, ito] 
                                                     - X_ij_re[l, k])
                             + G[l, k] * (X_ij_im[l, k])) 
    
        
        # CONSTRAINT (23)
        @constraint(model, I_sqr[l, k] 
                            == (G[l, k]^2 + B[l, k]^2) 
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
        @constraint(model, V_sqr[ifrom] - X_i_ij[l, k, ifrom] 
                            >= MIN_VOLTAGE^2 * (1 - alpha[l, k]))

        @constraint(model, V_sqr[ifrom] - X_i_ij[l, k, ifrom] 
                            <= MAX_VOLTAGE^2 * (1 - alpha[l, k]))

        @constraint(model, V_sqr[ito] - X_i_ij[l, k, ito] 
                            >= MIN_VOLTAGE^2 * (1 - alpha[l, k]))

        @constraint(model, V_sqr[ito] - X_i_ij[l, k, ito] 
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
    
    # ========================= Additional Constraints ========================

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

    #print(model)

    optimize!(model)

    solution_summary(model, verbose=false)

    if termination_status(model) == MOI.OPTIMAL
        var_values = Dict(  string(k) => value.(v) for (k, v) 
                            in object_dictionary(model) 
                            if (v isa AbstractArray{VariableRef} || 
                                v isa VariableRef)
                        )

        
        var_sets = Dict("X_i_ij"    => ["node"], 
                        "X_ij_re"   => ["line", "conductor"],
                        "X_ij_im"   => ["line", "conductor"],
                        "P_ij"      => ["line", "conductor"],
                        "P_ji"      => ["line", "conductor"],
                        "Q_ij"      => ["line", "conductor"],
                        "Q_ji"      => ["line", "conductor"],
                        "I_sqr"     => ["line", "conductor"],
                        "V_sqr"     => ["node"],
                        "P_G"       => ["node"],
                        "Q_G"       => ["node"],
                        "S_G"       => ["node"],
                        "alpha"     => ["line", "conductor"],
                        "beta"      => ["node"],
                        "x"         => ["line"]
                        )

        return model 

    elseif termination_status(model) == DUAL_INFEASIBLE
        println("problem unbounded")

    elseif termination_status(model) == MOI.INFEASIBLE
        println("problem infeasible")
    end
return
end