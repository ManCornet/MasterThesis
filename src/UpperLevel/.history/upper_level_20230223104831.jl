import XLSX
import DataFrames
using JuMP, Gurobi
#include("functions.jl")

using Formatting
function printx(args...)
    for arg in args
        if typeof(arg) == Float64
            _print_float(arg)
        elseif typeof(arg) == Array{Float64} || typeof(arg) == Vector{Float64}
            print("[")
            for ar in eachindex(arg)
                _print_float(arg[ar])
                if ar != length(arg) print("  ") end
            end
            print("]")
        elseif typeof(arg) == Matrix{Float64}
            print("[")
            for ar in eachindex(arg)
                _print_float(arg[ar])
                if ar != length(arg) print("  ") end
            end
            print("]")
        else print(arg)
        end
    end
    println()
end

function _print_float(arg)
    if abs(arg)<1e-5 print("0")
    elseif iszero(arg-floor.(arg)) print(Int.(round.(arg)))
    elseif abs(arg) > 1e3 print(Int.(round.(arg)))
    else
        for i in 2:-1:-8
            if abs(arg) > 10.0^i
                print(sprintf1("%.$(3-i)f",arg))
                break
            end
        end
    end
end

# PARAMETERS 
HOURS_IN_A_YEAR = 8760 # DO NOT CHANGE
BASE_VOLTAGE = 15 # kV
BASE_POWER = 1e3 # kVA
BASE_CURRENT = BASE_POWER / BASE_VOLTAGE # A
BASE_ADMITANCE = 1e-3 * BASE_CURRENT / BASE_VOLTAGE # Siemens
MIN_VOLTAGE = 0.95 # pu
MAX_VOLTAGE = 1.05 # pu

# PREPROCESSING

XLSX_FILE_PATH = "model_2S2H.xlsx"

CONDUCTOR_COST_PER_MM2_PER_KM = 200

# extracts the lines and the ending nodes
df_line = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "line"))
L = size(df_line)[1]
println("L: $L")
line_ends = Dict(l => (df_line.from_bus[l], df_line.to_bus[l]) for l in 1:L) # check bus/lines indexes from the .xlsx file
line_length = df_line.length_km # km
println("line_ends: $line_ends")
println("line_length: $line_length")

# extracts the buses (nodes)
df_bus = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "bus"))
N = size(df_bus)[1]
println("N: $N")

# links the lines to the nodes
Omega_sending = Dict(n => [] for n in 1:N)
Omega_receiving = Dict(n => [] for n in 1:N)
for l in 1:L
    push!(Omega_sending[line_ends[l][1]], l)
    push!(Omega_receiving[line_ends[l][2]], l)
end

println("Omega_sending: $Omega_sending")
println("Omega_receiving: $Omega_receiving")

# extracts the lines properties (conductance & susceptance)
K = 3 # TODO select conductors to consider ?
df_conductor = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "line_std_types"))
max_current = Dict(k => [0.0 for l in 1:L] for k in 1:K) # absolute, pu
conductance = Dict(k => [0.0 for l in 1:L] for k in 1:K) # absolute, pu
susceptance = Dict(k => [0.0 for l in 1:L] for k in 1:K) # absolute, pu
line_cost = Dict(k => [0.0 for l in 1:L] for k in 1:K)  # EUR/km
function process_conductors(conductor_idx, line_idx)
    max_i = df_conductor.max_i_ka[conductor_idx]*1e3
    max_current[conductor_idx][line_idx] = max_i / BASE_CURRENT
    
    r = line_length[line_idx] * df_conductor.r_ohm_per_km[conductor_idx]
    x = line_length[line_idx] * df_conductor.x_ohm_per_km[conductor_idx]
    y = 1/(r+im*x) / BASE_ADMITANCE
    
    conductance[conductor_idx][line_idx] = real(y) 
    susceptance[conductor_idx][line_idx] = imag(y)

    section = df_conductor.q_mm2[conductor_idx]
    line_cost[conductor_idx][line_idx] = section * line_length[line_idx] * CONDUCTOR_COST_PER_MM2_PER_KM
end

for k in 1:K
    for l in 1:L
        process_conductors(k,l)
    end
end


println("max_current: $max_current")
println("conductance: $conductance")
println("susceptance: $susceptance")
println("line_cost: $line_cost")

# TODO adapt xlsx. file to extract this data
n_s_init = 0 # Should be <= n_s
n_s = 2 # susbstation nodes are numbered to correspond to the first n_s nodes in the network
S_rating_init = 0*ones(n_s) / BASE_POWER
S_rating_max = 200*ones(n_s) / BASE_POWER
# substation_cost = [1, 1, 1]*1e0 # EUR small value => should build substations and no line: OK.
substation_cost = 1e5*[10,1]#ones(n_s)  # EUR


df_load = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "load"))
P_demand = [zeros(n_s); df_load.p_mw * 1e3 / BASE_POWER]
Q_demand = [zeros(n_s); df_load.q_mvar * 1e3 / BASE_POWER]


# OPTIMIZATION MODEL

# Objective function related data
K_l = 1
K_s = 1
line_loss = 0.01 # phi_l
loss_cost = 1 # c_l
substation_utilization = 0.01 # phi_S
substation_operation_cost = ones(n_s) # C_i^v
interest_rate_losses = 0.02
interest_rate_substation = 0.02
# TODO understand substation related costs and capacity.

# setup solvers for Gurobi
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "TimeLimit", 100)
set_optimizer_attribute(model, "Presolve", 0)

# Build the model 

@variable(model, I_squared[1:K, 1:L]) # squared current
@variable(model, P_s[1:N] >= 0)
@variable(model, Q_s[1:N])
@variable(model, S_s[1:N]) # Auxiliary variable
@variable(model, P_conductor[1:K, 1:L]) # CHANGEMENT ICI
@variable(model, Q_conductor[1:K, 1:L]) # CHANGEMENT ICI
@variable(model, MIN_VOLTAGE^2 <= v_squared[1:N] <= MAX_VOLTAGE^2)
@variable(model, x[1:L], Bin)
@variable(model, X_i_ij[1:K, 1:N, 1:L] >= 0)
@variable(model, X_ij_re[1:K, 1:L] >= 0)
@variable(model, X_ij_im[1:K, 1:L])
@variable(model, alpha[1:K, 1:L], Bin)
@variable(model, beta[1:N], Bin)
@variable(model, objective_terms[1:4])

# It is easier to define P_s and Q_s for all nodes although they should be zero where it is not possible to put a substation.
for i = n_s+1:N
    fix(P_s[i], 0.0; force=true)
    fix(Q_s[i], 0.0)
    fix(S_s[i], 0.0)
end

for l in 1:L, k in 1:K, i in 1:N
    if i in line_ends[l]
        continue
    end
    fix(X_i_ij[k, i, l], 0.0; force=true)
end

#1 -> means constraint (1) in paper from Jabr (Polyhedral formulations ...)
# @objective(model, Min, K_l * sum(alpha[k, l] * line_cost[k][l] * line_length[l] for k in 1:K, l in 1:L)
#                        + K_s * sum(beta[i] * substation_cost[i] for i in 1:n_s)
#                        + HOURS_IN_A_YEAR * (1 + interest_rate_losses) * line_loss * sum((P_s[i] - P_demand[i]) for i in 1:N) * BASE_POWER
#                        + HOURS_IN_A_YEAR * (1 + interest_rate_substation) * substation_utilization * sum((P_s[i]^2 + Q_s[i]^2) for i in 1:n_s)* BASE_POWER^2)
@objective(model, Min, sum(objective_terms) )
@constraint(model, objective_terms[1] >= K_l * sum(alpha[k, l] * line_cost[k][l] * line_length[l] for k in 1:K, l in 1:L))
@constraint(model, objective_terms[2] >= K_s * sum(beta[i] * substation_cost[i] for i in 1:n_s))
@constraint(model, objective_terms[3] >= HOURS_IN_A_YEAR * (1 + interest_rate_losses) * loss_cost * line_loss * sum((P_s[i] - P_demand[i]) for i in 1:N) * BASE_POWER)
@constraint(model, objective_terms[4] >= HOURS_IN_A_YEAR * (1 + interest_rate_substation) * substation_utilization * sum(substation_operation_cost[i]*(P_s[i]^2 + Q_s[i]^2) for i in 1:n_s)* BASE_POWER^2)

#2 is in the definition of alpha
#3
@constraint(model, line_constructed[l=1:L], x[l] == sum(alpha[k, l] for k in 1:K))
#4
@constraint(model, substation_capacity[i=1:n_s], S_s[i] == S_rating_init[i] + beta[i] * S_rating_max[i]) # Define S_s as an auxiliary variable
@constraint(model, substation_capacity_limit[i=1:n_s], [S_s[i], P_s[i], Q_s[i]] in SecondOrderCone()) # To state the SOC constraint in a simple way.
#11
@constraint(model, number_of_lines, sum(x[l] for l in 1:L) == N - n_s)#_init - sum(beta[i] for i in n_s_init+1:n_s))
#18
@constraint(model, active_balance[i=1:N], P_s[i] - P_demand[i] == sum(P_conductor[k, l] for k in 1:K, l in Omega_sending[i])
                                                                  -
                                                                  sum(P_conductor[k, l] for k in 1:K, l in Omega_receiving[i]))
#19 
@constraint(model, reacti_balance[i=1:N], Q_s[i] - Q_demand[i] == sum(Q_conductor[k, l] for k in 1:K, l in Omega_sending[i])
                                                                  -
                                                                  sum(Q_conductor[k, l] for k in 1:K, l in Omega_receiving[i]))
#20
@constraint(model, current_limit[k=1:K, l=1:L], I_squared[k, l] <= alpha[k, l] * max_current[k][l]^2)
for k = 1:K, l = 1:L
    ifrom = line_ends[l][1]
    ito = line_ends[l][2]
    #21
    @constraint(model, P_conductor[k, l] == conductance[k][l] * (X_i_ij[k, ifrom, l] - X_ij_re[k, l])
                                              -
                                              susceptance[k][l] * X_ij_im[k, l]) # X is hermitian => X_ij_im = - X_ji_im but we define only X_ij_im for simplicity.
    #22
    @constraint(model, Q_conductor[k, l] == -susceptance[k][l] * (X_i_ij[k, ifrom, l] - X_ij_re[k, l])
                                              -
                                              conductance[k][l] * X_ij_im[k, l]) # X is hermitian => X_ij_im = - X_ji_im but we define only X_ij_im for simplicity.
    #23
    @constraint(model, I_squared[k, l] == (conductance[k][l]^2 + susceptance[k][l]^2) * (X_i_ij[k, ifrom, l] + X_i_ij[k, ito, l] - 2 * X_ij_re[k, l]))
    #24
    @constraint(model, MIN_VOLTAGE * MAX_VOLTAGE * alpha[k, l] <= X_i_ij[k, ifrom, l]) # CHANGEMENT ICI
    @constraint(model, X_i_ij[k, ifrom, l] <= MAX_VOLTAGE^2 * alpha[k, l])
    @constraint(model, MIN_VOLTAGE * MAX_VOLTAGE * alpha[k, l] <= X_i_ij[k, ito, l]) # CHANGEMENT ICI
    @constraint(model, X_i_ij[k, ito, l] <= MAX_VOLTAGE^2 * alpha[k, l])
    #25
    @constraint(model, X_ij_re[k, l] <= MAX_VOLTAGE^2 * alpha[k, l]) # Lower bound is in variable definition.
    #26
    @constraint(model, - MAX_VOLTAGE^2 * alpha[k, l] <= X_ij_im[k, l])
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

# print(model)

#print(model)
optimize!(model)

# solution_summary(model, verbose=false)

printx("")
printx("objective  ", value.(objective_terms)," = ",Int(round(objective_value(model))))
printx("")
printx("x            ",value.(x))
printx("I_squared    ",value.(I_squared).*value.(alpha))
printx("P_conductor  ",value.(P_conductor).*value.(alpha))
printx("Q_conductor  ",value.(Q_conductor).*value.(alpha))
printx("")
printx("beta      ",value.(beta))
printx("P_s       ",value.(P_s).*value.(beta))
printx("Q_s       ",value.(Q_s).*value.(beta))
printx("V_squared ",value.(v_squared))