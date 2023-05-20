#-----------------------------------------------------------------------------
#
#                           - TFE : Bilevel DNEP - 
#                             University of Liege

#-----------------------------------------------------------------------------
# Created By  : Manon Cornet
# Created Date: Saturday March 18 2023
#
# parametes:
#   File with the parameters required to define the model
#
# =============================================================================
#                                   Imports
# =============================================================================
import XLSX
import DataFrames
using Dates
include("building_profiles.jl")

# =============================================================================
#                         Definition of global constants
# =============================================================================

const BASE_VOLTAGE    = 15                                 # [kV]
const BASE_POWER      = 1                                  # [kVA]
const BASE_CURRENT    = BASE_POWER / BASE_VOLTAGE          # [A]
const BASE_ADMITTANCE = BASE_CURRENT / BASE_VOLTAGE        # [S]
const MIN_VOLTAGE     = 0.95                               # [pu]
const MAX_VOLTAGE     = 1.05                               # [pu]
const DAYS_PER_YEAR   = 2
const GRANULARITY     = 5                                  # [min]
const MIN_PER_DAY     = 24 * 60                            # [min]


# =============================================================================
#                                   Functions
# =============================================================================

function process_conductors(df_conductor, line_length, max_current, 
                            conductance, susceptance, line_cost, 
                            conductor_idx, line_idx, conductor_cost)

    max_i = df_conductor.max_i_ka[conductor_idx]*1e3
    max_current[conductor_idx][line_idx] = max_i / BASE_CURRENT

    r = line_length[line_idx] * df_conductor.r_ohm_per_km[conductor_idx]
    x = line_length[line_idx] * df_conductor.x_ohm_per_km[conductor_idx]
    y = 1/(r+im*x) / BASE_ADMITTANCE

    conductance[conductor_idx][line_idx] = real(y) 
    susceptance[conductor_idx][line_idx] = imag(y)

    section = df_conductor.q_mm2[conductor_idx]
    line_cost[conductor_idx][line_idx] = (section * line_length[line_idx] 
                                                  * conductor_cost
                                         )
    return
end

# =============================================================================
#                         Definition of the parameters
# =============================================================================


# ======================== DN parameter definitions =======================

# ---- Path of the file containing the DN topology ----

XLSX_FILE_PATH = joinpath(splitdir(@__DIR__)[1], "network_models/model_2S2H.xlsx")

# ---- Line parameters ----

df_line     = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "line"))
L_size      = size(df_line)[1]     # number of lines
L           = 1:L_size             # Set of lines
L_init      = 0                    # Set of existing lines
line_length = df_line.length_km    # line lengths [km]
line_ends   = Dict(l => (df_line.from_bus[l], df_line.to_bus[l]) for l in L)

# ---- Bus parameters ----

df_bus = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "bus"))

N_size = size(df_bus)[1]  
N = 1:N_size                     # set of buses 
Ns_size = sum(df_bus.type .== "substation") # nb of substation buses
Nu_size = sum(df_bus.type .== "user")       # nb of load nodes

Ns = 1:Ns_size                              # set of substation buses
Nu = (1:Nu_size) .+ Ns_size                 # set of load buses          

# ---- Substation parameters ----

Ns_init       = 0  # Set of existing substations
S_rating_init = 0 .* ones(Ns_size) ./ BASE_POWER
S_rating_max  = 200 .* ones(Ns_size) ./ BASE_POWER
sub_install_cost = 700  # [€/kVA]
sub_expan_cost   = 500  # [€/kVA]

# ---- Link btw lines and nodes ----

Omega_sending   = Dict(n => [] for n in N)
Omega_receiving = Dict(n => [] for n in N)
for l in L
    push!(Omega_sending[line_ends[l][1]], l)
    push!(Omega_receiving[line_ends[l][2]], l)
end

# ---- Line properties ---- 

df_conductor = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "line_std_types"))

K_size = 3           # number of conductor types
K      = 1:K_size    # set of conductors
conductor_cost = 200 # [€/mm^2/km]

max_current = Dict(k => [0.0 for l in L] for k in K) # absolute, [pu]
conductance = Dict(k => [0.0 for l in L] for k in K) # absolute, [pu]
susceptance = Dict(k => [0.0 for l in L] for k in K) # absolute, [pu]
line_cost   = Dict(k => [0.0 for l in L] for k in K) # [€/km]

for k in K
    for l in L
        process_conductors(df_conductor, line_length, max_current, 
                            conductance, susceptance, line_cost, 
                            k, l, conductor_cost)
    end
end
# =========================== Time-dependent parameters =======================

#---- Planning horizon used in the simulation  ----

# idea here is to enter the period we want to cover, then we weight in the 
# objective function as a function of the time covered by the planning horizon

# User enters the simulation data: start date of the simulation and end date
start_date = DateTime(today() + Day(1))
end_date   = DateTime(start_date + Day(1))

# period covered
T = start_date:Minute(GRANULARITY*period):(end_date-Second(1))


nb_years   = 1 # nb of years in the planning horizon
Y          = 1:nb_years # Set containing the years taken into account
period     = 24 # nb of 5min. time steps to aggregate
delta_t    = GRANULARITY * period # [min]

# ---- Demand at buses ---- 

# Fix the PF
cos_phi = 0.98 # inductive power factor for all loads (see paper)
tan_phi = tan(acos(cos_phi)) # links P and Q

# Build active demand profiles
PROF_PATH = joinpath(splitdir(@__DIR__)[1], "Manchester_data/LCT_profiles")

load_files = [joinpath(PROF_PATH, "Summer_Load_Profiles.xlsx"),
              joinpath(PROF_PATH, "Winter_Load_Profiles.xlsx")
              ]

LCT_files = [joinpath(PROF_PATH, "Summer_PV_Profiles.xlsx"),
             joinpath(PROF_PATH, "Summer_EHP_Profiles.xlsx"),
             joinpath(PROF_PATH, "Winter_EV_Profiles.xlsx"),
             joinpath(PROF_PATH, "Winter_uCHP_Profiles.xlsx")
             ]

nb_daily_steps  = 288 ÷ period
nb_yearly_steps = nb_daily_steps * DAYS_PER_YEAR
T = 1:nb_yearly_steps

yearly_load_profiles = zeros(Float64, nb_years, nb_yearly_steps, N_size)
# Build an array with columns = year and row = period
for y in Y
    for d in 1:DAYS_PER_YEAR
        daily_profile = build_daily_load_profiles(load_files, nb_loads = Nu_size, period = period, seed = 1, summer = true)
        yearly_load_profiles[y, 1+(d-1)*nb_daily_steps:d*nb_daily_steps, Nu] = daily_profile
    end
end
#println(yearly_load_profiles)
#println(size(yearly_load_profiles))

# ---- PV profiles ----

# ======================== Objective function parameters =======================

# ---- DSO related parameters ----

DSO_INTEREST_RATE = 0.01      # Yearly interest
MAX_CO2_BUDGET    = -1        # [kgCO2] (negative number means no limit)
losses_cost       = 0.5*1e-2  # [€/kWh]

# ---- CO2 emissions ----

CO2_IMP_SUBSTATION = 0.2 # kgCO2/kWh of electrical energy imported from the transmission grid
CO2_FUEL           = 0.4 # kgCO2/kWh of electrical energy from genset
CO2_GENSET         = 50  # kgCO2/kW of genset capacity
CO2_LINE           = 100 # kgCO2/km of line constructed
CO2_SUBSTATION     = 100 # kgCO2/kVA of substation contructed

