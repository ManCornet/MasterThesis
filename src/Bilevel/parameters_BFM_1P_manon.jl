# -- IMPORTATION OF THE PACKAGES --
using DataFrames, XLSX, Statistics, Plots

#==============================================================================#
#                                   FUNCTIONS                                  #
#==============================================================================#


# -- LOADING THE EXCEL FILE CONTAINING THE NETWORK TOPOLOGY --

XLSX_FILE_PATH = joinpath(splitdir(@__DIR__)[1], "network_models/network_Nahman_Peric_2S23H.xlsx")

# -- FUNCTION THAT PROCESSES THE LINE PROPRETIES  --

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

    line_cost = Array{Float64}(undef, nb_lines, nb_cond) # [€/km]

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

# -- FUNCTIONS FOR COMPUTING THE NPV --

function PV_coeff(tau, lambda)
    return (1 - 1 / (1 + tau)^lambda) / tau
end

function CRF(tau, n)
    return (tau * (1 + tau)^n) / ((1 + tau)^n - 1)
end

# -- FUNCTION THAT NORMALIZES AN ARRAY BY MAX VALUE ALONG ONE OF THE DIMENSIONS --

function norm_col(array::AbstractArray; dims=1)
    return array ./ maximum(array, dims=dims)
end

# -- FUNCTION FOR PROCESSING THE TIME STEPS --

function process_time_steps(nb_time_steps_day::Integer; data_granularity::Integer=5)

    nb_min_day = 24 * 60
    max_nb_time_steps = nb_min_day ÷ data_granularity

    @assert nb_time_steps_day < max_nb_time_steps

    time_step_size = nb_min_day ÷ nb_time_steps_day
    nb_agg_time_steps = time_step_size ÷ data_granularity

    return max_nb_time_steps, time_step_size, nb_agg_time_steps 
end

# -- FUNCTION FOR AGGREGATING THE LOAD PROFILES TO HAVE TO HAVE LESS TIME STEPS --

function change_granularity(time_serie::AbstractArray; nb_agg_periods::Integer = 1)

    # ---- Checking arguments ----
    nb_periods = size(time_serie)[1]
    1 <= nb_agg_periods <= nb_periods || throw(
    ArgumentError("[change_granularity] granularity not in [1, $nb_periods]: $nb_agg_periods"))

    # ---- Changing granularity ----
    nb_rows, nb_columns = size(time_serie)
    new_time_serie = Matrix{Float64}(undef, 0, nb_columns)

    mod = nb_rows % nb_agg_periods 
    for i in 1:nb_agg_periods:(nb_rows - mod)
        new_time_serie = [new_time_serie; sum(time_serie[i:(i+nb_agg_periods-1), :], dims=1)]
    end
    
    return new_time_serie./nb_agg_periods
end

# -- FUNCTION FOR BUILDING DAILY PROFILES FROM MANCHESTER DATA  --

function build_daily_profiles(PROF_PATH::String;
                              nb_profiles::Integer = 1, seed = nothing)

    # In the load profiles, check the resolution of the data
    # The profiles are in kW

    # ---- Checking arguments ----

    1 <= nb_profiles <= 100 || throw(
    ArgumentError("[build_daily_load_profiles] nb_profiles not in [1, 100]: $nb_profiles"))

    1 <= nb_profiles <= 288 || throw(
    ArgumentError("[build_daily_load_profiles] number of aggregated periods not in [1, 288]: $period"))

    # ---- Creating the load profiles ----

    # Fixing the seed to choose randomly the load profiles selected
    if seed !== nothing 
        Random.seed!(seed)
        id_loads = Random.rand(1:100, nb_profiles)
    else 
        id_loads = 1:nb_profiles
    end

    # Load profiles
    profiles = XLSX.readdata(PROF_PATH, "Sheet1", "A1:CV288")
    
    nb_rows, ~ = size(profiles)
    # We only fetch a specific number of profiles < 100 (for nb of load nodes)
    load_profiles = reshape(profiles[:, id_loads], (nb_rows, length(id_loads)))

    return convert(Matrix{Float64}, load_profiles)
end

#==============================================================================#
#                            PART 1 : Upper-Level Parameters                   #
#==============================================================================#

# -- DEFINITION OF THE PER UNIT BASIS --

BASE_VOLTAGE = 34.5                             # [kV]
BASE_POWER = 1                                  # [MVA] check MANCHESTER_SCALING parameter
BASE_CURRENT = BASE_POWER / BASE_VOLTAGE        # [kA]
BASE_ADMITTANCE = BASE_CURRENT / BASE_VOLTAGE   # [S]
BASE_IMPEDANCE = 1 / BASE_ADMITTANCE            # [Ohm]
BASE_MONEY = 1                                  # [k€]
# MANCHESTER_SCALING = 0.001 # Manchester data in kW -> base power in MW
PV_SCALE_SUMMER_WINTER = 0.1 # winter PV production wr. summer production

# -- FETCH THE DATA FROM THE EXCEL SHEET --

df_bus = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "bus"))
df_line = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "line"))
df_cond = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "conductor"))

# -- LINE PARAMETERS DEFINITION --

L_size = DataFrames.nrow(df_line)                                   # Number of lines in the network
L = 1:L_size                                                        # Line set
line_ends = [(df_line.from_bus[l], df_line.to_bus[l]) for l in L]   # Indices of the line extremities
len_lines = convert(Vector{Float64}, df_line.length_km)             # Line lengths [km]
K_size = DataFrames.nrow(df_cond)                                   # Number of conductor types
K = 1:K_size                                                        # Set of conductors

# -- DEFINITION OF THE PHYSICAL QUANTITIES ASSOCIATED TO NETWORK LINES --

MAX_CURRENT, line_cost, R, X, G, B, = process_conductors(df_cond, len_lines, L_size)

# -- BUS PARAMETERS DEFINITION --

N_size = DataFrames.nrow(df_bus)                # Number of buses in the network
N = 1:N_size                                    # Buses set
Ns_size = sum(df_bus.type .== "substation")     # Number of substation buses
Ns_init_size = sum(df_bus.S_G_init_mva .> 0)    # Number of existing substations
Nu_size = sum(df_bus.type .== "user")           # Number of load nodes
Ns = 1:Ns_size                                  # Set of substation buses
Ns_init = 1:Ns_init_size                        # Set of exisitng substations
Ns_notinit = setdiff(Ns, Ns_init)               # Set of potential substations
Nu = (1:Nu_size) .+ Ns_size                     # Set of load buses
nodes_location = [df_bus.x, df_bus.y]

# -- DEFINITION OF THE PHYSICAL QUANTITIES ASSOCIATED TO NETWORK BUSES --

MIN_VOLTAGE = 0.9  # [pu]
MAX_VOLTAGE = 1.1  # [pu]

# -- SUBSTATION PARAMETERS DEFINITION --

S_rating_init = convert(Vector{Float64}, df_bus.S_G_init_mva[Ns]) ./ BASE_POWER # [pu]
S_rating_max = convert(Vector{Float64}, df_bus.S_G_max_mva[Ns]) ./ BASE_POWER   # [pu]
sub_install_cost = 1e3 / BASE_MONEY                                             # [pu]
# sub_op_cost = 0.1 * 1e-3 # k$/kVah^2

# -- LINK BTW LINES AND NODES --

Omega_sending = Dict(n => [] for n in N)
Omega_receiving = Dict(n => [] for n in N)
for l in L
    push!(Omega_sending[line_ends[l][1]], l)
    push!(Omega_receiving[line_ends[l][2]], l)
end

# -- DSO COSTS --

DSO_INTEREST_RATE = 0.05 # Yearly
AMORTIZATION_DSO = 50 # [years] amortization period for lines & substations installation
M = MAX_VOLTAGE^2

# nb_years_planning = 1
# tau = 0.1
line_loss_factor = 0.35     # phi_l : loss factor of lines
# sub_loss_factor = 0.35     # phi_s : cost per energy lost
# K_l = CRF(tau, 1)           # Capital recovery rate of line constructions
# K_s = CRF(tau, 1)           # Capital recovery rate of substation construction or reinforcement
loss_cost = 0.05       # [k€/MWh]
# tau_l = tau             # tau_l : interest rate of circuits
# tau_s = tau             # tau_s : interest rate of substations
# f_l = PV_coeff(tau_l, nb_years_planning)
# f_s = PV_coeff(tau_s, nb_years_planning)

#==============================================================================#
#                            PART 2 : Lower-Level Parameters                   #
#==============================================================================#

# -- FECTHING THE PROFILE FILES --

XLSX_SUMMER_LOAD_PATH = joinpath(splitdir(@__DIR__)[1], "Manchester_data/LCT_profiles/Summer_Load_Profiles.xlsx")
XLSX_WINTER_LOAD_PATH = joinpath(splitdir(@__DIR__)[1], "Manchester_data/LCT_profiles/Winter_Load_Profiles.xlsx")
XLSX_PV_PATH = joinpath(splitdir(@__DIR__)[1], "Manchester_data/LCT_profiles/Summer_PV_Profiles.xlsx")

# -- CREATING THE LOAD PROFILES --

# Definition of time constants
HOURS_PER_YEAR = 8760
DAYS_A_YEAR = 365
GRANULARITY = 5 # [min.] (granularity of the dataset that you use)

# Definition of the time_step of simulation
NB_SIMULATED_DAYS = 2
NB_TIME_STEPS_DAY = 24 # Nb of time steps per day
TOTAL_NB_TIME_STEPS = NB_SIMULATED_DAYS * NB_TIME_STEPS_DAY
T = 1:TOTAL_NB_TIME_STEPS # Set containing all time steps of simulation

# Get the max nb of time steps per day allowed by the data, the time step size and the number of 
# time steps to aggregate in the dataset to get the time step wanted
# ATTENTION, time_step_size is in minutes
max_nb_time_steps, time_step_size, nb_agg_time_steps =
                            process_time_steps(NB_TIME_STEPS_DAY; data_granularity=GRANULARITY)


# Get the max power demand repartition among the nodes in MVA

POWER_FACTOR = 0.9
S_D = convert(Vector{Float64}, df_bus.S_D_mva)[Nu]      # [MVA]
P_D = S_D * POWER_FACTOR                                # [MW]

# Get the load 
PEAK_POWER = sum(P_D) / BASE_POWER
summer_load = build_daily_profiles(XLSX_SUMMER_LOAD_PATH; nb_profiles=Nu_size) .* P_D'
winter_load = build_daily_profiles(XLSX_WINTER_LOAD_PATH; nb_profiles=Nu_size) .* P_D'

agg_summer_load = change_granularity(summer_load; nb_agg_periods=nb_agg_time_steps)'
agg_winter_load = change_granularity(winter_load; nb_agg_periods=nb_agg_time_steps)'

total_load = hcat(agg_summer_load, agg_winter_load)
load = vcat(summer_load, winter_load)'
scaling_factor = maximum(sum(load, dims=1))

total_load *= PEAK_POWER / scaling_factor

P_CONSUMPTION = [zeros(Ns_size, TOTAL_NB_TIME_STEPS); total_load]
Q_CONSUMPTION = P_CONSUMPTION .* tan(acos(POWER_FACTOR))
S_CONSUMPTION = sqrt.(P_CONSUMPTION.^2 + Q_CONSUMPTION.^2)


# -- CREATING THE PV PROFILES --

# User related data
PV_MAX_Q = 0.3 # max Q consumed/produced wrt P_peak
PV_SLOPE = -1 # linear slope of the upper bound from the PQ diagram of PV
POWER_FACTOR = 0.9 # P_CONSUMPTION = POWER_FACTOR * S_CONSUMPTION
AMORTIZATION_CONVERTER = 10  # [years] amortization period for users converter
AMORTIZATION_PV = 25  # [years] amortization period for users PV installation

# PV related data
PV_profiles = build_daily_profiles(XLSX_PV_PATH; nb_profiles=Nu_size)
PV_profiles = change_granularity(PV_profiles; nb_agg_periods=nb_agg_time_steps)'
PV_profiles ./= maximum(PV_profiles, dims=2)   # We suppose the peak in summer is the peak capacity installed

PV_PRODUCTION = zeros(Float64, N_size, TOTAL_NB_TIME_STEPS) # [%] so no unit scaling
PV_SCALES = [1, 0.1] # Scale for each days simulated
for i in 1:NB_SIMULATED_DAYS
    start_time = 1 + NB_TIME_STEPS_DAY*(i - 1)
    end_time  = start_time + NB_TIME_STEPS_DAY - 1
    PV_PRODUCTION[Nu, start_time:end_time] = PV_profiles * PV_SCALES[i]
end


# -- USER COSTS IN k€ --

GENSET_COST = 200 # EUR/kVA
FUEL_COST = 0.5 # EUR/kWh of electricity produced (embeds the conversion)
CONVERTER_COST = 200 # EUR/kVA converter installation
PV_COST = 500 # EUR/kWc PV installation

IMP_ELECTRICITY_ENRG_COST = 0.3 # EUR/kWh, cost of electricity bought to retailer
IMP_ELECTRICITY_DSO_COST = 0.15 # EUR/kWh, cost of electricity imported from DS
EXP_ELECTRICITY_ENRG_COST = 0.15 # EUR/kWh, cost of electricity injected in DS
EXP_ELECTRICITY_DSO_COST = 0.15 # EUR/kWh, cost of electricity injected in DS
GRID_CONNECTION_COST = 80 # EUR/kVA_year

TIME_STEP = time_step_size / 60


# Plot load consumptions


