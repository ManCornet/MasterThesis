using DataFrames, XLSX, Statistics, Plots

if (!@isdefined LAUNCH_SENSITIVITY_ANALYSIS) || (LAUNCH_SENSITIVITY_ANALYSIS != true)
    include("parameters_key_BFM_1P.jl")
end

# # OTHER PARAMETERS

# User related costs
GENSET_COST = 200 # EUR/kVA
FUEL_COST = 0.5 # EUR/kWh of electricity produced (embeds the conversion)
CONVERTER_COST = 200 # EUR/kVA converter installation
EXP_ELECTRICITY_DSO_COST = IMP_ELECTRICITY_DSO_COST # EUR/kWh, cost of electricity injected in DS

# DSO related costs
SUB_INSTALL_COST = 1e2 # [kEUR/MVA]

# Amortization periods
AMORTIZATION_CONVERTER = 10  # [years] amortization period for users converter
AMORTIZATION_PV = 25  # [years] amortization period for users PV installation
DSO_INTEREST_RATE = 0.06 # Yearly
AMORTIZATION_DSO = 50 # [years] amortization period for lines & substations installation

# Voltages
MIN_VOLTAGE = 0.95  # [pu]
MAX_VOLTAGE = 1.05  # [pu]


# # DIRECTORIES

root_dir = splitdir(@__DIR__)[1]
profiles_data_dir = joinpath(root_dir, "Manchester_data", "LCT_profiles")
network_data_dir = joinpath(root_dir, "network_models")

XLSX_FILE_PATH = joinpath(network_data_dir, "network_Nahman_Peric_2S23H.xlsx")
XLSX_SUMMER_LOAD_PATH = joinpath(profiles_data_dir, "Summer_Load_Profiles.xlsx")
XLSX_WINTER_LOAD_PATH = joinpath(profiles_data_dir, "Winter_Load_Profiles.xlsx")
XLSX_PV_PATH = joinpath(profiles_data_dir, "Summer_PV_Profiles.xlsx")
XLSX_EV_PATH = joinpath(profiles_data_dir, "Winter_EV_Profiles.xlsx")
XLSX_HP_PATH = joinpath(profiles_data_dir, "Winter_EHP_Profiles.xlsx")


# # PER UNITS
BASE_VOLTAGE = 34.5                           # [kV]
BASE_POWER = 1                                # [MVA]
BASE_CURRENT = BASE_POWER / BASE_VOLTAGE      # [kA]
BASE_ADMITTANCE = BASE_CURRENT / BASE_VOLTAGE # [S]
BASE_IMPEDANCE = 1 / BASE_ADMITTANCE          # [Ohm]
BASE_MONEY = 1                                # [kEUR]
PEAK_POWER /= BASE_POWER


# # NETWORK PARAMETERS
function process_conductors(
    df_cond::DataFrames.DataFrame,
    len_lines::Vector{Float64},
    nb_lines::Integer,
)

    nb_cond = DataFrames.nrow(df_cond)
    MAX_CURRENT = Array{Float64}(undef, nb_lines, nb_cond) # absolute, [pu]
    r = Array{Float64}(undef, nb_lines, nb_cond) # absolute, [pu]
    x = Array{Float64}(undef, nb_lines, nb_cond) # absolute, [pu]
    g = Array{Float64}(undef, nb_lines, nb_cond) # absolute, [pu]
    b = Array{Float64}(undef, nb_lines, nb_cond) # absolute, [pu]

    line_cost = Array{Float64}(undef, nb_lines, nb_cond) # [EUR/km]

    # Only take the first conductors of the list in the file
    for l in 1:nb_lines
        for k in 1:nb_cond
            MAX_CURRENT[l, k] = df_cond.max_i_ka[k] ./ BASE_CURRENT

            r[l, k] = len_lines[l] * df_cond.r_ohm_per_km[k] ./ BASE_IMPEDANCE
            x[l, k] = len_lines[l] * df_cond.x_ohm_per_km[k] ./ BASE_IMPEDANCE
            y = 1 / (r[l, k] + im * x[l, k])

            g[l, k] = real(y)
            b[l, k] = abs(imag(y))

            line_cost[l, k] = df_cond.cost_kdollars_per_km[k] * len_lines[l]
        end
    end
    return MAX_CURRENT, line_cost, r, x, g, b
end

df_bus = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "bus"))
df_line = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "line"))
df_cond = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "conductor"))

L_size = DataFrames.nrow(df_line)                                   # Number of lines in the network
L = 1:L_size                                                   # Line set
line_ends = [(df_line.from_bus[l], df_line.to_bus[l]) for l in L]   # Indices of the line extremities
len_lines = convert(Vector{Float64}, df_line.length_km)             # Line lengths [km]

K_size = DataFrames.nrow(df_cond)   # Number of conductor types
K = 1:K_size                   # Set of conductors
MAX_CURRENT, line_cost, R, X, G, B, = process_conductors(df_cond, len_lines, L_size)

N_size = DataFrames.nrow(df_bus)            # Number of buses in the network
N = 1:N_size                           # Buses set
Ns_size = sum(df_bus.type .== "substation") # Number of substation buses
Nu_size = sum(df_bus.type .== "user")       # Number of load nodes
Ns = 1:Ns_size                              # Set of substation buses
Nu = (1:Nu_size) .+ Ns_size                 # Set of load buses
nodes_location = [df_bus.x, df_bus.y]

S_rating_init = convert(Vector{Float64}, df_bus.S_G_init_mva[Ns]) ./ BASE_POWER # [pu]
S_rating_max = convert(Vector{Float64}, df_bus.S_G_max_mva[Ns]) ./ BASE_POWER # [pu]
SUB_INSTALL_COST /= BASE_MONEY   # [pu]

Omega_sending = Dict(n => [] for n in N)
Omega_receiving = Dict(n => [] for n in N)
for l in L
    push!(Omega_sending[line_ends[l][1]], l)
    push!(Omega_receiving[line_ends[l][2]], l)
end


# # ALGORITHMIC PARAMETERS

T = 1:Int(24 / TIME_STEP)
T_size = size(T)[1] # number of periods per day
DAYS_A_YEAR = 365
NB_PROFILES = 1

load_base = (XLSX.readxlsx(XLSX_SUMMER_LOAD_PATH)[1][:])'[1:Nu_size, :] .* df_bus.S_D_mva[Nu]
load_summer = copy(load_base)
load_EV = (XLSX.readxlsx(XLSX_EV_PATH)[1][:])'[1:Nu_size, :] .* df_bus.S_D_mva[Nu]
ELECTRIC_VEHICLES && (load_summer .+= load_EV)
load_size = size(load_base)[2]
divide_by = load_size ÷ T_size # integer divide
averaged_load = [mean(load_summer[j, i:i+divide_by-1]) for j in 1:size(load_base)[1], i in 1:divide_by:T_size*divide_by]
time_steps_lost = load_size - T_size * divide_by
(time_steps_lost > 0) && (println("$time_steps_lost \"5min steps\" have been lost due to TIME_STEP value"))

if WINTER
    load_winter = (XLSX.readxlsx(XLSX_WINTER_LOAD_PATH)[1][:])'[1:Nu_size, :] .* df_bus.S_D_mva[Nu]
    load_HP = (XLSX.readxlsx(XLSX_HP_PATH)[1][:])'[1:Nu_size, :] .* df_bus.S_D_mva[Nu]
    load_base = hcat(load_base, load_winter)
    ELECTRIC_VEHICLES && (load_winter .+= load_EV)
    HEAT_PUMPS && (load_winter .+= load_HP)
    averaged_load_winter = [mean(load_winter[j, i:i+divide_by-1]) for j in 1:size(load_winter)[1], i in 1:divide_by:T_size*divide_by]
    averaged_load = hcat(averaged_load, averaged_load_winter)
    T_size *= 2
    T = 1:T_size
    NB_PROFILES = 2
end

scaling_factor = maximum(sum(load_base, dims=1))
averaged_load *= PEAK_POWER / scaling_factor
S_CONSUMPTION = [zeros(Ns_size, T_size); averaged_load[1:Nu_size, :]]
P_CONSUMPTION = S_CONSUMPTION * COS_PHI
Q_CONSUMPTION = S_CONSUMPTION * sin(acos(COS_PHI))


# PV related data
pv = (XLSX.readxlsx(XLSX_PV_PATH)[1][:])'
pv_size = size(pv)[2]
# divide_by = pv_size ÷ T_size # integer divide
divide_by = pv_size ÷ Int(T_size / NB_PROFILES) # integer divide
# averaged_pv = [mean(pv[j, i:i+divide_by-1]) for j in 1:size(pv)[1], i in 1:divide_by:T_size*divide_by]
averaged_pv = [mean(pv[j, i:i+divide_by-1]) for j in 1:size(pv)[1], i in 1:divide_by:Int(T_size / NB_PROFILES)*divide_by]
averaged_pv ./= maximum(averaged_pv, dims=2)
PV_PRODUCTION = zeros(Float64, N_size, T_size)# [%] so no unit scaling
PV_PRODUCTION[Ns_size+1:end, 1:Int(T_size / NB_PROFILES)] = averaged_pv[1:Nu_size, :]
(NB_PROFILES > 1) && (PV_PRODUCTION[Ns_size+1:end, Int(T_size / NB_PROFILES)+1:end] = averaged_pv[1:Nu_size, :] * PV_SCALE_SUMMER_WINTER)


M = MAX_VOLTAGE^2
Ns_init = [1]
Ns_notinit = setdiff(Ns, Ns_init)
