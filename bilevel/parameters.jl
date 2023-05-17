using DataFrames, XLSX, Statistics, Plots

# Sources
# [1] Aten, M., & Ferris, R. (2009, June). Analysis of distribution losses and life cycle CO 2 emissions. In CIRED 2009-20th International Conference and Exhibition on Electricity Distribution-Part 1 (pp. 1-4). IET.

# User related costs
GENSET_COST = 200 # EUR/kVA
FUEL_COST = 0.5 # EUR/kWh of electricity produced (embeds the conversion)
CONVERTER_COST = 200 # EUR/kVA converter installation
PV_COST = 50 # EUR/kWc PV installation
IMP_ELECTRICITY_ENRG_COST = 0.3 # EUR/kWh, cost of electricity bought to retailer
IMP_ELECTRICITY_DSO_COST = 0.15 # EUR/kWh, cost of electricity imported from DS
EXP_ELECTRICITY_ENRG_COST = 0.15 # EUR/kWh, cost of electricity injected in DS
EXP_ELECTRICITY_DSO_COST = 0.15 # EUR/kWh, cost of electricity injected in DS
GRID_CONNECTION_COST = 300 # EUR/kVA_year

# User related data
PV_MAX_Q = 0.3 # max Q consumed/produced wrt P_peak
PV_SLOPE = -1 # linear slope of the upper bound from the PQ diagram of PV
POWER_FACTOR = 0.9 # P_CONSUMPTION = POWER_FACTOR * S_CONSUMPTION
AMORTIZATION_CONVERTER = 10  # [years] amortization period for users converter
AMORTIZATION_PV = 25  # [years] amortization period for users PV installation

# CO2 emissions -> TODO check all numbers and get reasonable values
CO2_IMP_SUBSTATION = 0.2 # kgCO2/kWh of electrical energy imported from the transmission grid (https://app.electricitymaps.com/map)
CO2_FUEL = 0.4 # kgCO2/kWh of electrical energy from genset
CO2_GENSET = 50 # kgCO2/kW of genset capacity
CO2_PV = 0 # kgCO2/kWp of PV capacity
CO2_CONVERTER = 0 # kgCO2/kVA of converter capacity
CO2_LINE = 100 # kgCO2/km of line constructed
CO2_SUBSTATION = 100 # kgCO2/kVA of substation constructed
# CO2_LOSS # kgCO2/kWh of wasted electrical energy

# DSO related
DSO_INTEREST_RATE = 0.01 # Yearly
SUBSTATION_INSTALLATION_COST = 500 # EUR/kVA
CONDUCTOR_COST_PER_MM2_PER_KM = 200
AMORTIZATION_DSO = 50 # [years] amortization period for lines & substations installation

# Algorithmic parameters
MAX_CO2_BUDGET = -1 # [kgCO2] negative number means no limit
# START_PERIOD = 1
# END_PERIOD = 24
# T = START_PERIOD:TIME_STEP:END_PERIOD # time set
T = 1:12
T_size = size(T)[1] # number of periods per day
TIME_STEP = 24 ÷ T_size # [hour] Duration of the period of time over which we assume powers are constant
DAYS_A_YEAR = 365
# TIME_SPAN = 1 # [years] simulation time span 

# Topology
# XLSX_FILE_PATH = splitdir(pwd())[1]*"\\network_models\\model_2S2H.xlsx"
# XLSX_FILE_PATH = splitdir(pwd())[1]*"\\bilevel\\network_3S7U_2023-03-29_16-42-08.xlsx"
# XLSX_FILE_PATH = splitdir(pwd())[1]*"\\bilevel\\network_8S22U_2023-03-29_16-36-34.xlsx"
XLSX_FILE_PATH = splitdir(pwd())[1] * "\\network_models\\network_Nahman_Peric_1S24H.xlsx"
# XLSX_FILE_PATH = splitdir(pwd())[1]*"\\network_models\\network_Nahman_Peric_2S23H.xlsx"
XLSX_LOAD_PATH = splitdir(pwd())[1] * "\\Manchester_data\\LCT_profiles\\Summer_Load_Profiles.xlsx"
XLSX_PV_PATH = splitdir(pwd())[1] * "\\Manchester_data\\LCT_profiles\\Summer_PV_Profiles.xlsx"


# extracts the lines and the ending nodes
df_line = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "line")...)
L_size = size(df_line)[1]
L = 1:L_size
line_ends = Dict(l => (df_line.from_bus[l], df_line.to_bus[l]) for l in L)
line_length = df_line.length_km # km
LINE_COST = line_length * 50 * CONDUCTOR_COST_PER_MM2_PER_KM # TO BE MODIFIED
S_MAX_LINE = 33 * ones(L_size)

# extracts the buses (nodes) 
df_bus = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "bus")...)
N_size = size(df_bus)[1]
N = 1:N_size
Ns_size = sum(df_bus.type .== "substation")
Ns = 1:Ns_size
Nu_size = sum(df_bus.type .== "user")
Nu = (1:Nu_size) .+ Ns_size
nodes_location = [df_bus.x, df_bus.y] # scale parameter for plot

Omega_sending = Dict(n => [] for n in N)
Omega_receiving = Dict(n => [] for n in N)
for l in L
    push!(Omega_sending[line_ends[l][1]], l)
    push!(Omega_receiving[line_ends[l][2]], l)
end

# df_load = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "load")...)
# # P_CONSUMPTION=[zeros(Ns_size); df_load.p_mw] # TO MODIFY FOR MORE TIME STEPS
# P_CONSUMPTION=[zeros(Ns_size); df_load.s_mw] # TO MODIFY FOR MORE TIME STEPS
S_MAX = 200 * ones(Ns_size)

# Consumptions related data
load = (XLSX.readxlsx(XLSX_LOAD_PATH)[1][:])'
load_size = size(load)[2]
divide_by = load_size ÷ T_size # integer divide
averaged_load = [mean(load[j, i:i+divide_by-1]) for j in 1:size(load)[1], i in 1:divide_by:T_size*divide_by]
P_CONSUMPTION = [zeros(Ns_size, T_size); averaged_load[1:Nu_size, :]]
S_CONSUMPTION = P_CONSUMPTION / POWER_FACTOR

# PV related data
pv = (XLSX.readxlsx(XLSX_PV_PATH)[1][:])'
pv_size = size(load)[2]
divide_by = pv_size ÷ T_size # integer divide
averaged_pv = [mean(pv[j, i:i+divide_by-1]) for j in 1:size(pv)[1], i in 1:divide_by:T_size*divide_by]
averaged_pv ./= maximum(averaged_pv, dims=2)
PV_PRODUCTION = [zeros(Ns_size, T_size); averaged_pv[1:Nu_size, :]]

# # Nahman_Peric_1S24H
# S_CONSUMPTION = [zeros(Ns_size); df_bus.load_kVA[Nu]] # kVA
# df_load_data = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "load_data")...)
# POWER_FACTOR = df_load_data.power_factor[1]
# P_CONSUMPTION = POWER_FACTOR * S_CONSUMPTION
# load_factor = df_load_data.load_factor[1]
# c_i = df_load_data.c_i[1]
# c_l = df_load_data.c_l[1]
# c_k = df_load_data.c_k[1]
# g = df_load_data.g[1]
# beta = 0.15 * load_factor + 0.85 * load_factor .^ 2
# df_line_data = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "line_data")...)
# repair_duration_h = df_line_data.repair_duration_h
# failure_rate_fl_per_yr = df_line_data.failure_rate_fl_per_km_yr .* line_length


PF_square = POWER_FACTOR^2
Q_CONSUMPTION = P_CONSUMPTION * sqrt((1 - PF_square) / PF_square)