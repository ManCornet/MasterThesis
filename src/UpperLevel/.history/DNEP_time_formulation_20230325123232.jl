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

function time_dependent_formulation(sets, costs, substation_param, conductor_param, 
    demand, delta_t)   

# ========================= Fetch the parameters ========================
# Sets
N, Ns, Ns_init, K, L, L_init, Y, T, Omega_sending, Omega_receiving = sets
Ns_not_init = setdiff(Ns, Ns_init)
Nu          = setdiff(N, Ns)
L_not_init  = setdiff(L, L_init)

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

for y in Y, t in T, i in Nu
fix(P_G[y, t, i], 0.0; force=true) 
fix(Q_G[y, t, i], 0.0)
fix(S_G[y, t, i], 0.0)
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
substation_apparent_power[y=Y, t=T, i=Ns], 
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
                 for l in Omega_sending[i], k in K)
                 + sum(P_cond_backward[y, t, l, k] 
                 for l in Omega_receiving[i], k in K)
)

# CONSTRAINT (19)
@constraint(model, 
reactive_balance[y=Y, t=T, i=N], 
Q_G[y, t, i] - Q_D[y, t, i] == sum(Q_cond_forward[y, t, l, k]
                for l in Omega_sending[i], k in K)
                + sum(Q_cond_backward[y, t, l, k]
                for l in Omega_receiving[i], k in K)
)

# CONSTRAINT (20)
@constraint(model, 
current_limit[y=Y, t=T, l=L, k=K], 
I_squared[y, t, l, k] <= alpha[l, k] * max_current[k][l]^2
)

for y in Y, t in T, l in L, k in K 
ifrom = line_ends[l][1]
ito = line_ends[l][2]

# CONSTRAINT (21)
@constraint(model, P_cond_forward[y, t, l, k] 
== 
conductance[k][l] * (X_i_ij[y, t, l, k, ifrom] - X_ij_re[y, t, l, k])
- susceptance[k][l] * X_ij_im[y, t, l, k]
)

@constraint(model, P_cond_backward[y, t, l, k] 
==
conductance[k][l] * (X_i_ij[y, t, l, k, ito] - X_ij_re[y, t, l, k])
- susceptance[k][l] * (- X_ij_im[y, t, l, k])) 

# CONSTRAINT (22)
@constraint(model, Q_cond_forward[y, t, l, k] 
== 
- susceptance[k][l] * (X_i_ij[y, t, l, k, ifrom] - X_ij_re[y, t, l, k])
- conductance[k][l] * X_ij_im[y, t, l, k])

@constraint(model, Q_cond_backward[y, t, l, k] 
== 
- susceptance[k][l] * (X_i_ij[y, t, l, k, ito] - X_ij_re[y, t, l, k])
- conductance[k][l] * (- X_ij_im[y, t, l, k])) 


# CONSTRAINT (23)
@constraint(model, I_squared[y, t, l, k] 
== 
(conductance[k][l]^2 + susceptance[k][l]^2) 
* (X_i_ij[y, t, l, k, ifrom] + X_i_ij[y, t, l, k, ito] 
- 2 * X_ij_re[y, t, l, k]))


# CONSTRAINT (24)
@constraint(model, MIN_VOLTAGE^2 * alpha[l, k] <= X_i_ij[y, t, l, k, ifrom])
@constraint(model, MAX_VOLTAGE^2 * alpha[l, k] >= X_i_ij[y, t, l, k, ifrom])
@constraint(model, MIN_VOLTAGE^2 * alpha[l, k] <= X_i_ij[y, t, l, k, ito])
@constraint(model, MAX_VOLTAGE^2 * alpha[l, k] >= X_i_ij[y, t, l, k, ito])

# CONSTRAINT (25)
@constraint(model, X_ij_re[y, t, l, k] <= MAX_VOLTAGE^2 * alpha[l, k])

# CONSTRAINT (26)
@constraint(model, X_ij_im[y, t, l, k] <= MAX_VOLTAGE^2 * alpha[l, k])
@constraint(model, X_ij_im[y, t, l, k] >= - MAX_VOLTAGE^2 * alpha[l, k])

# CONSTRAINT (27)
@constraint(model,   V_squared[y, t, ifrom] - X_i_ij[y, t, l, k, ifrom] 
>= MIN_VOLTAGE^2 * (1 - alpha[l, k])
)
@constraint(model,   V_squared[y, t, ifrom] - X_i_ij[y, t, l, k, ifrom] 
<= MAX_VOLTAGE^2 * (1 - alpha[l, k])
)
@constraint(model,   V_squared[y, t, ito] - X_i_ij[y, t, l, k, ito] 
>= MIN_VOLTAGE^2 * (1 - alpha[l, k])
)
@constraint(model,   V_squared[y, t, ito] - X_i_ij[y, t, l, k, ito] 
<= MAX_VOLTAGE^2 * (1 - alpha[l, k])
)
end

# CONSTRAINT (28)
@constraint(model, 
cone_28[y=Y, t=T, l=L, k=K], 
[X_i_ij[y, t, l, k, line_ends[l][1]] / 2, 
X_i_ij[y, t, l, k, line_ends[l][2]],
X_ij_re[y, t, l, k], X_ij_im[y, t, l, k]] 
in RotatedSecondOrderCone())

# CONSTRAINT (29)
@constraint(model, 
losses[y=Y, t=T],
active_losses[y, t] 
== sum(I_squared[y, t, l, k]/conductance[k][l] for l in L, k in K)
)
#print(model)

optimize!(model)

solution_summary(model, verbose=true)

if termination_status(model) == MOI.OPTIMAL
#include("export_xlsx.jl")
dict_output = Dict(k => typeof(v) for (k, v) in object_dictionary(model) if v is a )
println(dict_output)
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

#check_rotated_cones(cone_28)
#check_cones(substation_capacity_limit) 

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