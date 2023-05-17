using JuMP
using Printf

#using Juniper

# Parameters 
#=
HOURS_IN_A_YEAR = 8760 # DO NOT CHANGE
BASE_VOLTAGE = 15 # kV
BASE_POWER = 1e3 # kVA
BASE_CURRENT = BASE_POWER / BASE_VOLTAGE # A
BASE_ADMITANCE = 1e-3 * BASE_CURRENT / BASE_VOLTAGE # Siemens
MIN_VOLTAGE = 0.95 # pu
MAX_VOLTAGE = 1.05 # pu
=#
#include("toy_data.jl")
include("read_xlsx_net.jl")

# Objective function related data
#K_l = 1
#K_s = 1
#line_loss = 1 #??
#substation_utilization = 1
#interest_rate_losses = 0.02
#interest_rate_substation = 0.02

# TODO understand substation related costs and capacity.

# setup solvers for Pajarito
# using HiGHS # Does not handle cone constraints
# using Hypatia # Handles cone constraints
# using Pajarito # Provides the MI functionality.
# oa_solver = optimizer_with_attributes(HiGHS.Optimizer,
#      MOI.Silent() => true,
#     "mip_feasibility_tolerance" => 1e-8,
#     "mip_rel_gap" => 1e-6,
# )
# conic_solver = optimizer_with_attributes(Hypatia.Optimizer, 
#     MOI.Silent() => true,
# )
# opt = optimizer_with_attributes(Pajarito.Optimizer,
#     "time_limit" => 60, 
#     "oa_solver" => oa_solver, 
#     "conic_solver" => conic_solver,
# )
# model = Model(opt)

# setup solvers for Gurobi
using Gurobi
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "TimeLimit", 100)
set_optimizer_attribute(model, "Presolve", 0)

# Build the model 

@variable(model, I_squared[1:K, 1:L] >= 0) # squared current
@variable(model, P_s[1:N] >= 0)
@variable(model, Q_s[1:N])
@variable(model, S_s[1:n_s]) # Auxiliary variable
@variable(model, P_conductor_f[1:K, 1:L]) # Forward direction
@variable(model, P_conductor_b[1:K, 1:L]) # Backward direction
@variable(model, Q_conductor_f[1:K, 1:L]) # Forward direction
@variable(model, Q_conductor_b[1:K, 1:L]) # Backward direction
@variable(model, MIN_VOLTAGE^2 <= v_squared[1:N] <= MAX_VOLTAGE^2)
@variable(model, x[1:L], Bin)
@variable(model, X_i_ij[1:K, 1:N, 1:L] >= 0) # TODO only for l connected to i, normally, because all others are zero.
@variable(model, X_ij_re[1:K, 1:L] >= 0)
@variable(model, X_ij_im[1:K, 1:L])
@variable(model, alpha[1:K, 1:L], Bin)
@variable(model, beta[1:N], Bin)

# It is easier to define P_s and Q_s for all nodes although they should be zero where it is not possible to put a substation.
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

#1 -> means constraint (1) in paper from Jabr (Polyhedral formulations ...)
@objective(model, Min, K_l * sum(alpha[k, l] * line_cost[k][l] * line_length[l] for k in 1:K, l in 1:L)
                       +
                       K_s * sum(beta[i] * substation_cost[i] for i in 1:n_s)
                       + (1 + interest_rate_losses) * line_loss * sum((P_s[i] - P_demand[i]) for i in 1:N) * BASE_POWER
                       + (1 + interest_rate_substation) * substation_utilization * sum((P_s[i]^2 + Q_s[i]^2)  for i in 1:n_s)* BASE_POWER^2) 

#2 is in the definition of alpha
#3
@constraint(model, line_constructed[l=1:L], x[l] == sum(alpha[k, l] for k in 1:K)) # Test alternative sum(alpha, dims=1)
#4
@constraint(model, substation_capacity[i=1:n_s], S_s[i] == S_rating_init[i] + beta[i] * S_rating_max[i]) # Define S_s as an auxiliary variable
@constraint(model, substation_capacity_limit[i=1:n_s], [S_s[i], P_s[i], Q_s[i]] in SecondOrderCone()) # To state the SOC constraint in a simple way.
#11
@constraint(model, number_of_lines, sum(x[l] for l in 1:L) == N - n_s_init - sum(beta[i] for i in n_s_init+1:n_s))
#18
@constraint(model, active_balance[i=1:N], P_s[i] - P_demand[i] == sum(P_conductor_f[k, l] for k in 1:K, l in Omega_sending[i])
                                                                  +
                                                                  sum(P_conductor_b[k, l] for k in 1:K, l in Omega_receiving[i]))
#19 
@constraint(model, reacti_balance[i=1:N], Q_s[i] - Q_demand[i] == sum(Q_conductor_f[k, l] for k in 1:K, l in Omega_sending[i])
                                                                  +
                                                                  sum(Q_conductor_b[k, l] for k in 1:K, l in Omega_receiving[i]))
#20
@constraint(model, current_limit[k=1:K, l=1:L], I_squared[k, l] <= alpha[k, l] * max_current[k][l]^2)
for k = 1:K, l = 1:L
    ifrom = line_ends[l][1]
    ito = line_ends[l][2]
    #21
    @constraint(model, P_conductor_f[k, l] == conductance[k][l] * (X_i_ij[k, ifrom, l] - X_ij_re[k, l])
                                              -
                                              susceptance[k][l] * X_ij_im[k, l])
    @constraint(model, P_conductor_b[k, l] == conductance[k][l] * (X_i_ij[k, ito, l] - X_ij_re[k, l])
                                              -
                                              susceptance[k][l] * (-X_ij_im[k, l])) # X is hermitian => X_ij_im = - X_ji_im but we define only X_ij_im for simplicity.
    #22
    @constraint(model, Q_conductor_f[k, l] == -susceptance[k][l] * (X_i_ij[k, ifrom, l] - X_ij_re[k, l])
                                              -
                                              conductance[k][l] * X_ij_im[k, l])
    @constraint(model, Q_conductor_b[k, l] == -susceptance[k][l] * (X_i_ij[k, ito, l] - X_ij_re[k, l])
                                              -
                                              conductance[k][l] * (-X_ij_im[k, l])) # X is hermitian => X_ij_im = - X_ji_im but we define only X_ij_im for simplicity.
    #23
    @constraint(model, I_squared[k, l] == (conductance[k][l]^2 + susceptance[k][l]^2) * (X_i_ij[k, ifrom, l] + X_i_ij[k, ito, l] - 2 * X_ij_re[k, l]))
    #24
    @constraint(model, MIN_VOLTAGE^2 * alpha[k, l] <= X_i_ij[k, ifrom, l])
    @constraint(model, X_i_ij[k, ifrom, l] <= MAX_VOLTAGE^2 * alpha[k, l])
    @constraint(model, MIN_VOLTAGE^2 * alpha[k, l] <= X_i_ij[k, ito, l])
    @constraint(model, X_i_ij[k, ito, l] <= MAX_VOLTAGE^2 * alpha[k, l])
    #25
    @constraint(model, X_ij_re[k, l] <= MAX_VOLTAGE^2 * alpha[k, l]) # Lower bound is in constraint definition.s
    #26
    @constraint(model, -MAX_VOLTAGE^2 * alpha[k, l] <= X_ij_im[k, l])
    @constraint(model, X_ij_im[k, l] <= MAX_VOLTAGE^2 * alpha[k, l])
    #27
    @constraint(model, MIN_VOLTAGE^2 * (1 - alpha[k, l]) <= v_squared[ifrom] - X_i_ij[k, ifrom, l])
    @constraint(model, v_squared[ifrom] - X_i_ij[k, ifrom, l] <= MAX_VOLTAGE^2 * (1 - alpha[k, l])) # TODO add name to constraint
    @constraint(model, MIN_VOLTAGE^2 * (1 - alpha[k, l]) <= v_squared[ito] - X_i_ij[k, ito, l])
    @constraint(model, v_squared[ito] - X_i_ij[k, ito, l] <= MAX_VOLTAGE^2 * (1 - alpha[k, l]))  # TODO add name to constraint
end

#28
@constraint(model, cone_28[k=1:K, l=1:L], [X_i_ij[k, line_ends[l][1], l] / 2, X_i_ij[k, line_ends[l][2], l],
    X_ij_re[k, l], X_ij_im[k, l]] in RotatedSecondOrderCone())

print(model)

optimize!(model)

solution_summary(model, verbose=true)

function check_rotated_cones()
    println("Rotated cone constraints that are not tight:")
    for k in 1:K, l in 1:K
        x = value(cone_28[k, l])
        slack = 2 * x[1] * x[2] - (x[3]^2 + x[4]^2)
        if abs(slack) > 1e-3
            println(k, l)
        end
    end
end

function check_cones()
    println("Cone constraints that are not tight:")
    for n in 1:n_s
        x = value(substation_capacity_limit[n])
        slack = x[1]^2 - (x[2]^2 + x[3]^2)
        if abs(slack) > 1e-3
            println(n)
        end
    end
end


I_squared   = transpose(value.(I_squared))
V_squared   = value.(v_squared)
x           = value.(x)
alpha       = transpose(value.(alpha))
beta        = value.(beta)
X_ij_re     = transpose(value.(X_ij_re))
X_ij_im     = transpose(value.(X_ij_im))
X_i_ij      = value.(X_i_ij)
P_s         = value.(P_s)
Q_s         = value.(Q_s)
P_cond_forward = transpose(value.(P_conductor_f))
P_cond_backward = transpose(value.(P_conductor_b))
Q_cond_forward = transpose(value.(Q_conductor_f))
Q_cond_backward = transpose(value.(Q_conductor_b))
obj         = objective_value(model)
time        = solve_time(model)

function print_table(table, table_name, unity)

    end_string = 12 - length(table_name)
    @printf("%s", table_name)
    for i in 1:end_string
        @printf(" ")
    end 
    @printf(": [ ")
    for (index, value) in enumerate(table)
        if typeof(value) == ComplexF64
            @printf("%.4f + i%.4f ", real(value), imag(value)) 

        elseif typeof(value) == Float64
            if value - floor(value) == 0
                @printf("%d ", Int(round(value))) 
            else
                @printf("%.4f ", value) 
            end

        else typeof(value) == Int64
            @printf("%d", value) 
        end
        if index < length(table)
            @printf(" ")
        end
    end
    @printf("]")
    @printf(" [%s]\n", unity) 
end

function print_2Dtable(table, table_name, unity)
    
    end_string = 12 - length(table_name)
    @printf("%s", table_name)
    for i in 1:end_string
        @printf(" ")
    end 
    @printf(": [ ")

    s = size(table)
    for i in 1:s[1]
        for j in 1:s[2]
            value = table[i, j]
            if typeof(value) == ComplexF64
                @printf("%.4f + i%.4f ", real(value), imag(value)) 
            elseif typeof(value) == Float64
                if value - floor(value) < 1e-6
                    @printf("%d ", Int(round(value)))
                else
                    @printf("%.4f ", value) 
                end
            else typeof(value) == Int64
                @printf("%d ", value) 
            end
        end
        if i < s[1]
            @printf("; ")
        end
    end
    @printf("]")
    @printf(" [%s]\n", unity)  
end

@printf("-> Number of buses: %d \n\n", N)
@printf("-> Number of lines: %d \n\n", L)
@printf("-> Number of conductor types: %d \n\n", K)
@printf("-> Binary variables: \n\n")

print_2Dtable(alpha, "Alpha" , "-") # 2D dimensions
print_table(x, "x", "-")
print_table(beta, "Beta", "-")
@printf("\n-> Continuous variables: \n\n")
print_2Dtable(I_squared, "I_squared" , "pu^2")
print_table(V_squared, "V_squared", "pu^2")
print_table(abs.(P_s .* beta), "P_s", "kW")
print_table(abs.(Q_s .* beta), "Q_s", "kVar")
print_2Dtable(abs.(P_cond_forward .* alpha), "P_cond_forward", "kW")
print_2Dtable(abs.(P_cond_backward .* alpha), "P_cond_backward", "kW")
print_2Dtable(abs.(Q_cond_forward .* alpha), "Q_cond_forward", "kVar")
print_2Dtable(abs.(Q_cond_backward .* alpha), "Q_cond_backward", "kVar")
#print_2Dtable(X, "X", "pu^2")
println(obj)

check_rotated_cones()
check_cones()