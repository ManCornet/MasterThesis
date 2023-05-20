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
import XLSX, DataFrames, Dates
include("building_profiles.jl")

# Maybe use the package "MetaGraph.jl"
# =============================================================================
#                         Definition of global constants
# =============================================================================
# Idea: automatic tests with several parameters that are changed and modified

const BASE_VOLTAGE    = 15                                 # [kV]
const BASE_POWER      = 100                                # [kVA] 
const BASE_CURRENT    = BASE_POWER / BASE_VOLTAGE          # [A]
const BASE_ADMITTANCE = BASE_CURRENT / BASE_VOLTAGE*1e3    # [S]
const MIN_VOLTAGE     = 0.95                               # [pu]
const MAX_VOLTAGE     = 1.05                               # [pu]
const GRANULARITY     = 5                                  # [min]
const MIN_PER_DAY     = 24 * 60                            # [min]

# =============================================================================
#                         Definition of structures
# =============================================================================
""" Bus

    Fields:
    -------
        - id : identificator of the node
        - x  : x-coordinate of the node
        - y  : y-coordinate of the node
"""


mutable struct Bus
    id::Int
    coordinates::Union{Nothing, Coord}
    x::Union{Nothing, Float64}
    y::Union{Nothing, Float64}            #
    V_magn::Union{Nothing, Float64}       # Voltage magnitude
    base_voltage::Union{Nothing, Float64} # in kV
    voltage_limits::Tuple{}
    V_max::Float64
    V_min::Float64
    
    function Bus(id::Int64, x::Float64, y::Float64) 
        if id < 0 
            throw(DomainError("[Node]: The id of a node cannot be negative")) 
        end

        return new(id, x, y) 
    end
end


""" Substation 

Fields:
-------
    - rating_init: initial rating of a substation   [kVA]
    - rating_max : maximum rating of the substation [kVA]
    - location   : node       
    - V_must be fixed (reference voltage)                      
"""
mutable struct Substation <: Bus
    id ::Int64
    x  ::Float64
    y  ::Float64
    V  ::Float64 
    P_G::Float64 
    Q_D::Float64
    Q_G::Float64
    Q_G::Float64
    base_kV::Float64
    V_max::Float64
    V_min::Float64
end

""" Line (mutable because the "built field" can change !)

    Fields:
    -------
        _ id            : indentificator of the line
        - built         : Boolean indicating if this line has been built or not
        - line_ends     : tuple containing the sending and receiving ends of the line
        - line_length   : length of the line [km]
        - conductor     : the conductor of the line
        - line_cost     : the total cost of the line [€]
        - max_current   : maximum current that can flow through the line [A]
        - conductance   : conductance of the line [S]
        - susceptance   : susceptance of the line [S]
"""
mutable struct Line
    id          ::Int64
    built       ::Bool
    from_bus    ::Bus 
    to_bus      ::Bus
    line_length ::Float64
    conductor   ::Conductor
    max_i       ::Float64
    r           ::Float64 
    x           ::Float64


    function Line(  id::Int64, line_ends::Tuple{Node, Node}, 
                    line_length::Float64, conductor::Conductor
                ) 
        if id < 0 || line_length < 0 
            throw(DomainError("""[Line]: The id of a line/line 
                                         length cannot be negative."""))
        end
        line_cost = conductor.q_mm2 * line_length * conductor.cost

        max_current = conductor.max_i 

        r = line_length * conductor.r_per_km 
        x = line_length * conductor.x_per_km
        y = 1/(r+im*x)

        conductance = real(y)
        susceptance = imag(y)
       
        return new( id, false, line_ends, line_length, conductor, line_cost,
                    max_current, conductance, susceptance) 
    end
end

""" Conductor

    Fields:
    -------
        _ name      : name of the conductor
        - r_per_km  : resistance [ohm/km]
        - x_per_km  : reactance  [ohm/km]
        - max_i     : maximum current [kA]
        - cost      : the cost of the conductor per mm2 per km [€/mm^2/km]
"""
struct Conductor 
    name    ::String            
    r_per_km::Float64
    x_per_km::Float64
    max_i   ::Float64
    q_mm2   ::Float64
    cost    ::Float64 
    function Conductor( name::String, r_per_km::Float64, x_per_km::Float64, 
                        max_i::Float64, q_mm2::Float64)
        if qmm2 < 0 || r_ohm_per_km < 0 || x_ohm_per_km < 0 || max_i_ka < 0
            throw(DomainError("""[Conductor]: The section of a 
                                              conductor/resistance/reactance/max current 
                                              cannot be negative."""))
        end

        return new(name, r_per_km, x_per_km, max_i, q_mm2, 200) 
    end
end


# Solution is a constructor o



""" Distribution network structure

    Fields:
    -------
        - nb_lines : The number of lines in the network
        - nb_buses : The number of buses in the network
        - nb_subs  : The number of substations in the network
        - lines    : The lines of the network 
        - nodes    : The nodes of the network
        - subs     : The substations of the network
"""
struct DistributionNetwork
    nb_lines::Int64 
    nb_nodes::Int64
    nb_subs ::Int64
    lines   ::Array{Line}
    nodes   ::Vector{Node}
    subs    ::Vector{Substation}
    function DistributionNetwork(nb_lines::Int64, nb_nodes::Int64, nb_subs::Int64,
                    lines::Vector{Line}, nodes::Vector{Node}, 
                    subs::Vector{Substation})
        if nb_lines < 0 || nb_nodes < 0 || nb_subs < 0
            throw(DomainError("""[Network]: The number of lines/nodes/substations 
                                            in the network cannot be negative"""))
        elseif length(lines) != nb_lines || length(nodes) != nb_nodes || length(subs) != nb_subs
            throw(DomainError("""[Network]: The size of the array containing the 
                                            lines/nodes/substations must be equal to 
                                            the nb od lines/nodes/substations"""))
        end

        return new(nb_lines, nb_nodes, nb_subs, lines, nodes, subs) 
    end
end

# =============================================================================
#                         Definition of the parameters
# =============================================================================

# ---- Path of the file containing the DN topology ----
XLSX_FILE_PATH = joinpath(splitdir(@__DIR__)[1], "network_models/network_Nahman_Peric_1S24H.xlsx")

# ======================== DN parameter definitions =======================
# ---- Bus parameters ----

df = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "bus"))
N  = size(df_bus)[1] # set of buses 
Ns = sum(df_bus.type .== "substation") # nb of substation buses

nodes = [i => Node(i, df.x, df.y) for i in 1:N]

# ---- Line parameters ----
df = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "line"))
L  = size(df)[1]

# Fetch the lines 
lines = Array{Line}(undef, L, K)
for l in 1:L
    for k in 1:K
        lines[l, k] = Line(l, (nodes[df.from_bus[l]], nodes[df.to_bus[l])], )
    end
end

line_ends = [(df.from_bus[l], df.to_bus[l]) for l in 1:L]  

Omega_sending   = Dict(n => [] for n in N)
Omega_receiving = Dict(n => [] for n in N)
for l in L
    push!(Omega_sending[line_ends[l][1]], l)
    push!(Omega_receiving[line_ends[l][1]], l)
end








# ---- Link btw lines and nodes ----





# ---- Conductor properties ---- 
# Idea: the number of conductor is something to enter in the command line
K = 3

df = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "line_std_types"))
conductors = [Conductor(df.Name, 
                        df.r_ohm_per_km, 
                        df.x_ohm_per_km, 
                        df.max_i_ka, 
                        df.q_mm2) 
            for k in 1:K]









# ---- Substation parameters ----

Ns_init       = 0  # Set of existing substations
S_rating_init = 0 .* ones(Ns_size) ./ BASE_POWER
S_rating_max  = 200 .* ones(Ns_size) ./ BASE_POWER

substations = [Substation(S_rating_init[i], S_rating_max[i]) for i in 1:Ns]
sub_install_cost = 700  # [€/kVA]
sub_expan_cost   = 500  # [€/kVA]



# Do it in another way : 



# =========================== Time-dependent parameters =======================

#---- Planning horizon used in the simulation  ----

# User enters the simulation data: start date of the simulation and end date
#start_date = DateTime(today() + Day(1))
#end_date   = DateTime(start_date + Day(1))

# period covered
#T = start_date:Minute(GRANULARITY*period):(end_date-Second(1))

period     = 24 # nb of 5min. time steps to aggregate
nb_daily_steps  = (24*60 / 5) ÷ period
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

load_profiles = zeros(Float64, nb_daily_steps, N_size)
# Build an array with columns = year and row = period
load_profiles[:, Nu] = build_daily_load_profiles(load_files, nb_loads = Nu_size, period = period, seed = 1, summer = true)
T = 1:nb_daily_steps

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

