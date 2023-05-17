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

# ======================= 1. Time-dependent formulation =======================

# Formulation that does not take into account the

function time_dependent_formulation(sets, costs, substation_param, conductor_param, 
                                    demand, delta_t)   

    # ATTENTION FOR THE FORMULATION MUST SPECIFY RANGE
    # ========================= Fetch the parameters ========================
    # Sets
    N, Ns, Ns_init, K, L, L_init, T, Omega_sending, Omega_receiving = sets
    Ns_not_init = setdiff(Ns, Ns_init)
    Nu          = setdiff(N, Ns)
    L_not_init  = setdiff(L, L_init)

    N_size, Ns_size , K_size, L_size, T_size = length(N), length(Ns), length(K), length(L), length(T)
    # Costs
    sub_expan_cost, sub_install_cost, line_cost, losses_cost, DSO_INTEREST_RATE = costs

    # Substation parameters
    S_rating_init, S_rating_max = substation_param

    # Conductor parameters
    max_current, conductance, susceptance = conductor_param

    # Demand profiles
    P_D, tan_phi = demand
    println(P_D[:,1])
    println(P_D[:,2])
    println(P_D[:,3])
    println(P_D[:,4])
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
    
    # small prob here 
    # i want arrays so table indexed with set beginning at 1
    @variables( model,   
                begin 
                V_squared[1:T_size, 1:N_size]
                I_squared[1:T_size, 1:L_size, 1:K_size]
                P_sending[1:T_size, 1:L_size, 1:K_size] 
                P_receiving[1:T_size, 1:L_size, 1:K_size] 
                Q_sending[1:T_size, 1:L_size, 1:K_size]
                Q_receiving[1:T_size, 1:L_size, 1:K_size]
                X_i_ij[1:T_size, 1:L_size, 1:K_size, 1:N_size]  >= 0
                X_ij_re[1:T_size, 1:L_size, 1:K_size]           >= 0
                X_ij_im[1:T_size, 1:L_size, 1:K_size]
                P_G[1:T_size, 1:N_size]                         >= 0
                Q_G[1:T_size, 1:N_size]
                S_G[1:T_size, 1:N_size]                         >= 0
                S_allocated[1:Ns_size]                          >= 0
                P_losses[1:T_size, 1:L_size]                    >= 0          
                alpha[1:L_size, 1:K_size]                       , Bin
                #beta[1:Ns_size]                                , Bin 
                x[1:L_size]                                     , Bin
                end
            )

    # need to post_process beta and alpha
    for t in T, i in Nu
        fix(S_G[t, i], 0.0; force=true)
        fix(P_G[t, i], 0.0; force=true) 
        fix(Q_G[t, i], 0.0)
    end

    # Maybe there is another way to write this variables ?
    for t in T, l in L, k in K, i in N 
        if i in line_ends[l]
            continue
        else
            fix(X_i_ij[t, l, k, i], 0.0; force=true)
        end
    end

    # ========================== Objective function =============================
    
    #steps_per_day  = size(P_D)[1]
    #steps_per_year = steps_per_day * DAYS_IN_A_YEAR
    # everything is done per day then needs to adapt for a year

    @objective(model, Min,  sum(alpha[l, k] * line_cost[l, k] for l in L, k in K) 
                            + sum(S_allocated[i] * sub_install_cost for i in Ns_not_init)
                            + sum(S_allocated[i] * sub_expan_cost for i in Ns_init)
                            + DAYS_IN_A_YEAR * (1/(1 + DSO_INTEREST_RATE)) 
                            * sum(P_losses[t, l] * losses_cost * delta_t/60 for t in T, l in L)
                )
    # ============================== Constraints ================================ 

    # CONSTRAINT (2) -> In definition of alpha

    # CONSTRAINT (3)
    @constraint(model, line_constructed[l=L], x[l] == sum(alpha[l, k] for k in K))

    # CONSTRAINT (4)

    @constraint(model, 
                substation_apparent_power[t=T, i=Ns], 
                [S_G[t, i]; P_G[t, i]; Q_G[t, i]] in SecondOrderCone()
    )

    @constraint(model, 
                substation_capacity_limit[t=T, i=Ns], 
                S_G[t, i] <= S_rating_init[i] + S_allocated[i]
    )

    # Problem here because S_allocated is not bounded below 
    # So we can have beta = 1 and S_allocated =0 because we want to minimize it 
    # Solution ? beta should be 0 when S_allocated is 0
    # S_rating max represents the max total apparent power 
    @constraint(model, 
                substation_capacity[i=Ns], 
                S_allocated[i] <= S_rating_max[i] - S_rating_init[i]
    ) 
    
    #=
    @constraint(model, 
                substation_voltage[t=T, i=Ns], 
                V_squared[t, i] == 1
    )=#

    # due to the voltages => issue

    # CONSTRAINT (11)
    @constraint(model, number_of_lines, sum(x[l] for l in L) == N_size - Ns_size)
    
    # CONSTRAINT POWER BALANCE => put this since not in objective 
    # this constraint would probably be in the lower level
  
    @constraint(model, 
                power_balance1[t=T], 
                sum(P_G[t, i] - P_D[t, i] for i in N) >= 0
    )

    @constraint(model, 
                power_balance2[t=T], 
                sum(Q_G[t, i] - Q_D[t, i] for i in N) >= 0
    )

    # CONSTRAINT (18)
    @constraint(model, 
                active_balance[t=T, i=N], 
                P_G[t, i] - P_D[t, i] == sum(P_sending[t, l, k]  
                                        for l in Omega_sending[i], k in K)
                                        - sum(P_receiving[t, l, k] 
                                        for l in Omega_receiving[i], k in K)
    )

    # CONSTRAINT (19)
    @constraint(model, 
                reactive_balance[t=T, i=N], 
                Q_G[t, i] - Q_D[t, i] == sum(Q_sending[t, l, k]
                                        for l in Omega_sending[i], k in K)
                                        - sum(Q_receiving[t, l, k]
                                        for l in Omega_receiving[i], k in K)
    )

    # CONSTRAINT (20)
    @constraint(model, 
                current_limit[t=T, l=L, k=K], 
                I_squared[t, l, k] <= alpha[l, k] * max_current[l, k]^2
    )


    for t in T, l in L, k in K 
        ifrom = line_ends[l][1]
        ito = line_ends[l][2]

        # CONSTRAINT (21)
        @constraint(model, P_sending[t, l, k] 
                           == conductance[l, k] * (X_i_ij[t, l, k, ifrom] 
                                                   - X_ij_re[t, l, k])
                           - susceptance[l, k] * X_ij_im[t, l, k]
                    )

        @constraint(model, P_receiving[t, l, k] 
                           == - conductance[l, k] * (X_i_ij[t, l, k, ito] 
                                                   - X_ij_re[t, l, k])
                            - susceptance[l, k] * X_ij_im[t, l, k]
                    ) 
        

        # CONSTRAINT (22)
        @constraint(model, Q_sending[t, l, k] 
                            == - susceptance[l, k] * (X_i_ij[t, l, k, ifrom] 
                                                     - X_ij_re[t, l, k])
                            - conductance[l, k] * X_ij_im[t, l, k]
                    )

        @constraint(model, Q_receiving[t, l, k] 
                            == susceptance[l, k] * (X_i_ij[t, l, k, ito] 
                                                     - X_ij_re[t, l, k])
                             - conductance[l, k] * X_ij_im[t, l, k]
                    ) 


        # CONSTRAINT (23)
        @constraint(model, I_squared[t, l, k] 
                            == (conductance[l, k]^2 + susceptance[l, k]^2) 
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
                [X_i_ij[t, l, k, line_ends[l][1]] / 2; 
                X_i_ij[t, l, k, line_ends[l][2]];
                X_ij_re[t, l, k]; 
                X_ij_im[t, l, k]] 
                in RotatedSecondOrderCone())
    
    # CONSTRAINT (28)bis
    @constraint(model, 
                cone_28bis[t=T, l=L, k=K], 
                [V_squared[t, line_ends[l][2]]/2;
                I_squared[t, l, k];  
                P_sending[t, l, k];
                Q_sending[t, l, k]] 
                in RotatedSecondOrderCone())

    # CONSTRAINT (28)bisbis
    @constraint(model, 
                cone_28bisbis[t=T, l=L, k=K], 
                [V_squared[t, line_ends[l][1]] /2; 
                I_squared[t, l,k]; 
                P_receiving[t, l, k];
                Q_receiving[t, l, k]] 
                in RotatedSecondOrderCone())

    # CONSTRAINT (29)
    @constraint(model, 
                active_losses1[t=T, l=L],
                P_losses[t, l] >= sum(P_sending[t, l, k] - P_receiving[t, l, k] for k in K))

    @constraint(model, 
                active_losses2[t=T, l=L],
                P_losses[t, l] >= sum(P_receiving[t, l, k] - P_sending[t, l, k] for k in K))

    #print(model)
    optimize!(model)

    #solution_summary(model, verbose=false)

    if termination_status(model) == MOI.OPTIMAL
        var_values = Dict(string(k) => value.(v) for (k, v) in object_dictionary(model) 
                        if (v isa AbstractArray{VariableRef} || v isa VariableRef))

        var_sets = Dict("X_i_ij" => ["time","node"],
                        "X_ij_re" => ["time", "line", "conductor"],
                        "X_ij_im" => ["time", "line", "conductor"],
                        "P_sending" => ["time", "line", "conductor"],
                        "P_receiving" => ["time", "line", "conductor"],
                        "Q_sending" => ["time", "line", "conductor"],
                        "Q_receiving" => ["time", "line", "conductor"],
                        "I_squared" => ["time", "line", "conductor"],
                        "V_squared" => ["time", "node"],
                        "P_G" => ["time", "node"],
                        "Q_G" => ["time", "node"],
                        "S_G" => ["time", "node"],
                        "alpha" => ["line", "conductor"],
                        "S_allocated" => ["node"],
                        "P_losses" => ["time", "line"],
                        #"beta" => ["node"],
                        "x" => ["line"]
                        )
        check_rotated_cones(cone_28, T, L, K)
        check_cones(substation_apparent_power, T, Ns)
        return var_values, var_sets

    elseif termination_status(model) == DUAL_INFEASIBLE
        println("problem unbounded")

    elseif termination_status(model) == MOI.INFEASIBLE
        println("problem infeasible")
    end
end