import XLSX
import DataFrames
using JuMP
using Gurobi
include("functions.jl")

# Parameters 
HOURS_IN_A_YEAR = 8760 # DO NOT CHANGE
BASE_VOLTAGE = 15 # kV
BASE_POWER = 1e3 # kVA
BASE_CURRENT = BASE_POWER / BASE_VOLTAGE # A
BASE_ADMITTANCE = 1e-3 * BASE_CURRENT / BASE_VOLTAGE # Siemens
MIN_VOLTAGE = 0.95 # pu
MAX_VOLTAGE = 1.05 # pu


# Definition of the topology. To understand, not to change.
XLSX_FILE_PATH = "model_2S2H.xlsx"

CONDUCTOR_COST_PER_MM2_PER_KM = 200

# extracts the lines and the ending nodes
df_line = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "line")...)
L = size(df_line)[1]
line_ends = Dict(l => (df_line.from_bus[l], df_line.to_bus[l]) for l in 1:L)
line_length = df_line.length_km # km

# extracts the buses (nodes) 
df_bus = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "bus")...)
N = size(df_bus)[1]

# links the lines to the nodes
Omega_sending = Dict(n => [] for n in 1:N)
Omega_receiving = Dict(n => [] for n in 1:N)
for l in 1:L
    push!(Omega_sending[line_ends[l][1]], l)
    push!(Omega_receiving[line_ends[l][2]], l)
end

# extracts the lines properties (conductance & susceptance)
df_conductor = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "line_std_types")...)
K = 1  # You may play with this parameter, or even modify the code to chose specific conductors in the XLSX file
max_current = Dict(k => [0.0 for l in 1:L] for k in 1:K) # absolute, pu
conductance = Dict(k => [0.0 for l in 1:L] for k in 1:K) # absolute, pu
susceptance = Dict(k => [0.0 for l in 1:L] for k in 1:K) # absolute, pu
line_cost = Dict(k => [0.0 for l in 1:L] for k in 1:K)  # EUR/km
function process_conductors(conductor_idx, line_idx)
    max_i = df_conductor.max_i_ka[conductor_idx]*1e3
    max_current[conductor_idx][line_idx] = max_i / BASE_CURRENT
    
    r = line_length[line_idx] * df_conductor.r_ohm_per_km[conductor_idx]
    x = line_length[line_idx] * df_conductor.x_ohm_per_km[conductor_idx]
    y = 1/(r+im*x) / BASE_ADMITTANCE
    
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


n_s_init = 0 # Should be <= n_s
n_s = 2 # susbstation nodes are numbered to correspond to the first n_s nodes in the network
S_rating_init = 0*ones(n_s) / BASE_POWER
S_rating_max = 200*ones(n_s) / BASE_POWER
substation_cost = 1e5*[1,1]#1e5*ones(n_s)  # EUR


df_load = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "load")...)
P_demand = [zeros(n_s); df_load.p_mw * 1e3 / BASE_POWER]
Q_demand = [zeros(n_s); df_load.q_mvar * 1e3 / BASE_POWER]


# Optimization model

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
set_optimizer_attribute(model, "NonConvex", 2)

# Build the model 

@variable(model, I_squared[1:K, 1:L]) # squared current
@variable(model, P_s[1:N] >= 0)
@variable(model, Q_s[1:N])
@variable(model, S_s[1:n_s]) # Auxiliary variable
@variable(model, P_conductor[1:K, 1:L])
@variable(model, Q_conductor[1:K, 1:L])
@variable(model, x[1:L], Bin)
@variable(model, alpha[1:K, 1:L], Bin)
@variable(model, beta[1:N], Bin)
@variable(model, V_re[1:N])
@variable(model, V_im[1:N])
@variable(model, expression_1[1:K, 1:L])
@variable(model, expression_2[1:K, 1:L])
@variable(model, expression_3[1:K, 1:L])
@variable(model, objective_terms[1:4])

# It is easier to define P_s and Q_s for all nodes although they should be zero where it is not possible to put a substation.
for i = n_s+1:N
    fix(P_s[i], 0.0; force=true) #force=true is required to override the bound >=0 given in P_s definition
    fix(Q_s[i], 0.0)
end


#1 -> means constraint (1) in paper from Jabr (Polyhedral formulations ...)
# @objective(model, Min, K_l * sum(alpha[k, l] * line_cost[k][l] * line_length[l] for k in 1:K, l in 1:L)
                    #    + K_s * sum(beta[i] * substation_cost[i] for i in 1:n_s)
                    #    + HOURS_IN_A_YEAR * (1 + interest_rate_losses) * line_loss * sum((P_s[i] - P_demand[i]) for i in 1:N) * BASE_POWER
                    #    + HOURS_IN_A_YEAR * (1 + interest_rate_substation) * substation_utilization * sum((P_s[i]^2 + Q_s[i]^2) for i in 1:n_s)* BASE_POWER^2)
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
#5
@constraint(model, active_balance[i=1:N], P_s[i] - P_demand[i] == sum(alpha[k, l]*P_conductor[k, l] for k in 1:K, l in Omega_sending[i])
                                                                  -
                                                                  sum(alpha[k, l]*P_conductor[k, l] for k in 1:K, l in Omega_receiving[i]))
#6 
@constraint(model, reacti_balance[i=1:N], Q_s[i] - Q_demand[i] == sum(alpha[k, l]*Q_conductor[k, l] for k in 1:K, l in Omega_sending[i])
                                                                  -
                                                                  sum(alpha[k, l]*Q_conductor[k, l] for k in 1:K, l in Omega_receiving[i]))
for k = 1:K, l = 1:L
    ifrom = line_ends[l][1]
    ito = line_ends[l][2]
    # tensions & courants sont sous-contraints, fortes variations dans les résultats !!
    @constraint(model, expression_1[k,l] == V_re[ifrom]^2 + V_im[ifrom]^2 - V_re[ifrom]*V_re[ito] - V_im[ifrom]*V_im[ito])
    @constraint(model, expression_2[k,l] == V_re[ifrom]*V_im[ito] - V_im[ifrom]*V_re[ito])
    # @constraint(model, expression_3[k,l] == V_re[ito]^2   + V_im[ito]^2   - V_re[ifrom]*V_re[ito] - V_im[ifrom]*V_im[ito])
    #7 real
    @constraint(model, P_conductor[k,l] == conductance[k][l] * expression_1[k,l] + susceptance[k][l] * expression_2[k,l])
    #7 imag
    @constraint(model, Q_conductor[k,l] == conductance[k][l] * expression_2[k,l] - susceptance[k][l] * expression_1[k,l])
    #9
    @constraint(model, I_squared[k, l] == (conductance[k][l]^2 + susceptance[k][l]^2) * (expression_1[k,l] + V_re[ito]^2   + V_im[ito]^2   - V_re[ifrom]*V_re[ito] - V_im[ifrom]*V_im[ito]))
    # @constraint(model, max_current[k][l]^2 >= (conductance[k][l]^2 + susceptance[k][l]^2) * (V_re[ifrom]^2+V_im[ifrom]^2+V_re[ito]^2+V_im[ito]^2-2*(V_re[ifrom]*V_re[ito]+V_im[ifrom]*V_re[ito])))

end

#8
@constraint(model, current_limit[k=1:K, l=1:L], alpha[k, l] * I_squared[k, l] <= max_current[k][l]^2)

#10
@constraint(model, [i=1:N], V_re[i]^2+V_im[i]^2 <= MAX_VOLTAGE^2)
@constraint(model, [i=1:N], V_re[i]^2+V_im[i]^2 >= MIN_VOLTAGE^2)

# print(model)

optimize!(model)

# solution_summary(model, verbose=false)

V_squared=value.(V_re).^2+value.(V_im).^2

printx("")
printx("objective  ", value.(objective_terms)," = ",Int(round(objective_value(model))))
printx("")
printx("x            ",value.(x))
printx("I_squared    ",value.(I_squared).*value.(alpha))
printx("P_conductor  ",value.(P_conductor).*value.(alpha))
printx("Q_conductor  ",value.(Q_conductor))
printx("")
printx("beta      ",value.(beta))
printx("P_s       ",value.(P_s).*value.(beta))
printx("Q_s       ",value.(Q_s).*value.(beta))
printx("V_re      ",value.(V_re))
printx("V_im      ",value.(V_im))
printx("V_squared ",V_squared)
printx("constraint 5 LHS  ", value.(P_s) - value.(P_demand))
printx("constraint 6 LHS  ", value.(Q_s) - value.(Q_demand))
