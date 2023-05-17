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

function time_dependent_formulation(network::Network, cost_functions::Cost, planning_horizon)  
#function time_dependent_formulation(sets, costs, substation_param, conductor_param, 
#                                    demand, delta_t)   

    # ======================== Definition of the sets ========================

    # Assumptions: the initial lines and initial substations are indexed first
    
    # ---- Line set ----
     L      = 1:network.nb_lines
     L_init = [line.id for line in network.lines if line.built]

    # ---- Node sets ----
    # All nodes
    N  = 1:network.nb_nodes

    # Substation nodes
    Ns          = 1:network.nb_substations
    Ns_init     = [sub.node.id for sub in network.sub if sub.built]
    Ns_not_init = setdiff(Ns, Ns_init)

    # User nodes
    Nu = setdiff(N, Ns)

    # ---- Substation param ----
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
    
    # small prob here 
    # i want arrays so table indexed with set beginning at 1
    @variable(model, I_squared[time=T, line=L, conductor=K], container = Array) 

    @variables( model,   
                begin 
                V_squared[T, N],             (container = Array)
                I_squared[T, L, K],          (container = Array)
                P_cond_forward[T, L, K],     (container = Array)
                P_cond_backward[T, L, K],    (container = Array)
                Q_cond_forward[T, L, K],     (container = Array)
                Q_cond_backward[T, L, K], (container=Array)
                X_i_ij[T, L, K, N]  >= 0, (container=Array)
                X_ij_re[T, L, K]    >= 0,
                X_ij_im[T, L, K]
                P_G[T, L, K]                         >= 0
                Q_G[T, L, K]
                S_G[T, L, K]
                S_alloc[N]                          >= 0
                P_loss[T]                         >= 0
                alpha[L, K]                       , Bin
                beta[N]                           , Bin 
                x[L]                              , Bin
                end
            )

    # S_allocated 

    #
    for i in Nu
        fix(S_alloc[i], 0.0; force=true)
        fix(beta[i], 0.0) 
    end
    
    for t in T, i in Nu
        fix(P_G[t, i], 0.0; force=true) 
        fix(Q_G[t, i], 0.0)
        fix(S_G[t, i], 0.0)
    end

    # Maybe there is another way to write this variables ?
    for t in T, l in L, k in K, i in N 
        if i in line_ends[l]
            continue
        end
        fix(X_i_ij[t, l, k, i], 0.0; force=true)
    end

    # ========================== Objective function =============================
    

    @objective(model, Min,  sum(alpha[l, k] * line_cost[l, k] * line_length[l] 
                                for k in K, l in L) 
                            + sum(S_alloc[i] * sub_install_cost  for i in Ns_not_init)
                            + sum(S_alloc[i] * sub_expan_cost for i in Ns)
                            + DAYS_IN_A_YEAR * sum((1/(1 + DSO_INTEREST_RATE)^(t ÷ DAYS_IN_A_YEAR) * 
                              active_losses[t] * losses_cost * delta_t/60 for t in T))
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
                S_G[t, i] <= S_rating_init[i] + S_alloc[i]
    )

    @constraint(model, 
                substation_capacity[i=Ns], 
                S_alloc[i] <= beta[i] * S_rating_max[i]
    ) 

    # CONSTRAINT (11)
    @constraint(model, number_of_lines, sum(x[l] for l in L) == length(N) - length(Ns))

    # CONSTRAINT (18)
    # node with id i
    @constraint(model, 
                active_balance[t=T, i=N], 
                P_G[t, i] - P_D[t, i] == sum(P_cond_forward[t, l, k]  
                                        for l in network.nodes[i].send_lines, k in K)
                                        + sum(P_cond_backward[t, l, k] 
                                        for l in network.nodes[i].rec_lines, k in K)
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
                I_squared[t, l, k] <= alpha[l, k] * max_current[l, k]^2
    )

    for t in T, l in L, k in K 
        ifrom = line_ends[l][1]
        ito = line_ends[l][2]

        # CONSTRAINT (21)
        @constraint(model, P_cond_forward[t, l, k] 
                           == conductance[l, k] * (X_i_ij[t, l, k, ifrom] 
                                                   - X_ij_re[t, l, k])
                           - susceptance[l, k] * X_ij_im[t, l, k]
        )

        @constraint(model, P_cond_backward[t, l, k] 
                           == conductance[l, k] * (X_i_ij[t, l, k, ito] 
                                                   - X_ij_re[t, l, k])
                            - susceptance[l, k] * (- X_ij_im[t, l, k])) 

        # CONSTRAINT (22)
        @constraint(model, Q_cond_forward[t, l, k] 
                            == - susceptance[l, k] * (X_i_ij[t, l, k, ifrom] 
                                                     - X_ij_re[t, l, k])
                            - conductance[l, k] * X_ij_im[t, l, k])

        @constraint(model, Q_cond_backward[t, l, k] 
                            == - susceptance[l, k] * (X_i_ij[t, l, k, ito] 
                                                     - X_ij_re[t, l, k])
                            - conductance[l, k] * (- X_ij_im[t, l, k])) 


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
                [X_i_ij[t, l, k, line_ends[l][1]] / 2, 
                X_i_ij[t, l, k, line_ends[l][2]],
                X_ij_re[t, l, k], X_ij_im[t, l, k]] 
                in RotatedSecondOrderCone())

    # CONSTRAINT (29)
    @constraint(model, 
                active_losses[t=T],
                P_loss[t] == sum(I_squared[t, l, k]/conductance[l, k] 
                                        for l in L, k in K))

    #print(model)
    optimize!(model)

    #solution_summary(model, verbose=false)

    if termination_status(model) == MOI.OPTIMAL

        var_values = Dict(string(k) => value.(v) for (k, v) in object_dictionary(model) 
                        if v isa AbstractArray{VariableRef} || v isa VariableRef)

        var_sets = Dict("X_i_ij" => ["time", "line", "conductor", "node"],
                        "X_ij_re" => ["time", "line", "conductor"],
                        "X_ij_im" => ["time", "line", "conductor"],
                        "P_cond_forward" => ["time", "line", "conductor"],
                        "P_cond_backward" => ["time", "line", "conductor"],
                        "Q_cond_forward" => ["time", "line", "conductor"],
                        "Q_cond_backward" => ["time", "line", "conductor"],
                        "I_squared" => ["time", "line", "conductor"],
                        "V_squared" => ["time", "node"],
                        "P_G" => ["time", "node"],
                        "Q_G" => ["time", "node"],
                        "S_G" => ["time", "node"],
                        "alpha" => ["line", "conductor"],
                        "S_alloc" => ["node"],
                        "P_loss" => ["time"],
                        "beta" => ["node"],
                        "x" => ["line"]
                        )

        return var_values, var_sets

    elseif termination_status(model) == DUAL_INFEASIBLE
        println("problem unbounded")

    elseif termination_status(model) == MOI.INFEASIBLE
        println("problem infeasible")
    end
end