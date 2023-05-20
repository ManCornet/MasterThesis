using DataFrames, XLSX, Statistics, Plots

# # KEY PARAMETERS

# ALGORITHMIC PARAMETERS
ACTIVATE_RADIALITY_CONSTRAINTS = true
TIME_STEP = 1 # [hour] Duration of the period of time over which we assume powers are constant
WEIGHT_OP_LIMITS = 1e-2

# CONSUMPTION
WINTER = true # extend the summer profile
ELECTRIC_VEHICLES = false # added to summer & winter profiles
HEAT_PUMPS = false # added to winter profile only
PEAK_POWER = 7 # [MVA] maximum power of all the apparent demands

# PV
MAX_PV_CAPACITY_PER_NODE = 0.4 # MVA 
PV_SCALE_SUMMER_WINTER = 0.1 # winter PV production wrt summer production
PV_MAX_Q = 0.3 # max Q consumed/produced wrt P_peak
PV_SLOPE = -1 # linear slope of the upper bound from the PQ diagram of PV
COS_PHI = 0.95 # P_CONSUMPTION = COS_PHI * S_CONSUMPTION

# DSO related costs
LOSS_COST = 0.1 # 0.7 * IMP_ELECTRICITY_ENRG_COST   # [kEUR/MWh]

# User related costs
PV_COST = 500 # EUR/kWc PV installation
IMP_ELECTRICITY_ENRG_COST = 0.3 # EUR/kWh, cost of electricity bought to retailer
IMP_ELECTRICITY_DSO_COST = 0.1 # EUR/kWh, cost of electricity imported from DS
EXP_ELECTRICITY_ENRG_COST = 0.1 # EUR/kWh, cost of electricity injected in DS
# EXP_ELECTRICITY_DSO_COST = IMP_ELECTRICITY_DSO_COST, see parameters_BFM_1P.jl
GRID_CONNECTION_COST = 80 # EUR/kVA_year


