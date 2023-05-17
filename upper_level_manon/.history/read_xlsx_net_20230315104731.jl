import XLSX
import DataFrames

#XLSX_FILE_PATH = "example_short.xlsx"

XLSX_FILE_PATH = "model_2S2H.xlsx"
CONDUCTOR_COST_PER_MM2_PER_KM = 200

df_line = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "line")...)
L = size(df_line)[1]
line_ends = Dict(l => (df_line.from_bus[l]+1, df_line.to_bus[l]+1) for l in 1:L)
line_length = df_line.length_km # km

df_bus = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "bus")...)
N = size(df_bus)[1]

Omega_sending = Dict(n => [] for n in 1:N)
Omega_receiving = Dict(n => [] for n in 1:N)
for l in 1:L
    push!(Omega_sending[line_ends[l][1]], l)
    push!(Omega_receiving[line_ends[l][2]], l)
end


df_conductor = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "line_std_types")...)

K = 3 # TODO select conductors to consider ?
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


df_load = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "load")...)
P_demand = pushfirst!(df_load.p_mw * 1e3 / BASE_POWER, 0.0) # Active power demand at each bus
Q_demand = pushfirst!(df_load.q_mvar * 1e3 / BASE_POWER, 0.0) # Reactive power demand at each bus

n_s_init = 1 # Should be <= n_s
n_s = 1 # susbstation nodes are numbered to correspond to the first n_s nodes in the network
S_rating_init = 0 .* ones(n_s) ./ BASE_POWER
S_rating_max  = 200 .* ones(n_s) ./ BASE_POWER
# substation_cost = [1, 1, 1]*1e0 # EUR small value => shoulf build substations and no line: OK.
substation_fixed_cost = 1e5 .* [10, 1]   # Fixed construction or reinforcement cost of substations [€] 
substation_op_cost    = ones(n_s)   # Substations operation cost [€/kVAh^2]   