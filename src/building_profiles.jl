#-----------------------------------------------------------------------------
#
#                           - TFE : Bilevel DNEP - 
#                             University of Liege
#
#-----------------------------------------------------------------------------
# Created By  : Manon Cornet
# Created Date: Saturday March 18 2023
#
# building_profiles:
#   File building the load and PV profiles coming from Manchester datasets
#
# =============================================================================
#                                   Imports
# =============================================================================
import XLSX, Random
using Statistics

# =============================================================================
#                                  Functions
# =============================================================================
function change_granularity(time_serie::AbstractArray; nb_agg_periods::Integer = 1)

    # ---- Checking arguments ----
    nb_periods = size(time_serie)[1]
    1 <= period <= nb_periods || throw(
    ArgumentError("[change_granularity] granularity not in [1, $nb_periods]: $period"))

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
                              nb_profiles::Integer = 1, nb_agg_periods::Integer = 1,
                              seed = nothing)

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

    # here change the granularity 
    if period > 1
        return change_granularity(load_profiles; nb_agg_periods)
    end

    return convert(Matrix{Float64}, load_profiles)
end

function build_daily_LCT_profiles(PROF_PATHS::Vector{String}, nb_LCT::Vector{Integer}; 
                                  period::Integer = 1, seed::Integer = 1, 
                                  PV::Bool = true, EHP::Bool = false,
                                  EV::Bool = false, uCHP::Bool = false
                                 )

    # ---- Checking arguments ----

    for i in nb_LCT
        1 <= i <= 100 || throw(
        ArgumentError("[build_daily_LCT_profiles] nb_LCT not in [1, 100]: $i"))
    end 

    1 <= period <= 288 || throw(
    ArgumentError("[build__daily_LCT_profiles] granularity not in [1, 288]: $period" ))
    
    # ---- Creating the LCT profiles ----

    # Fixing the seed to choose randomly the load profiles selected
    id_LCT = Vector{Integer}[]
    for i in nb_LCT
        Random.seed!(seed+i)
        rng           = 1:100
        id_LCT     = push!(id_LCT, Random.rand(rng, i))
        println("id_loads: $id_loads")
    end

    LCT_profiles = Matrix{Float64}[]
    if PV
        profiles = XLSX.readdata(PROF_PATHS[1], "Sheet1", "A1:CV288")
        LCT_profiles = push!(LCT_profiles, 
                             reshape(profiles[:, id_LCT[1]], (nb_rows, length(id_LCT[1])))
                            ) 

    end
    if EHP
        profiles = XLSX.readdata(PROF_PATHS[2], "Sheet1", "A1:CV288")
        LCT_profiles = push!(LCT_profiles,
                             reshape(profiles[:, id_LCT[2]], (nb_rows, length(id_LCT[2])))
                            )
    end  
    if EV 
        profiles = XLSX.readdata(PROF_PATHS[3], "Sheet1", "A1:CV288")
        LCT_profiles = push!(LCT_profiles,
                             reshape(profiles[:, id_LCT[3]], (nb_rows, length(id_LCT[3])))
                            )
    end
    if uCHP 
        profiles = XLSX.readdata(PROF_PATHS[4], "Sheet1", "A1:CV288")
        LCT_profiles = push!(LCT_profiles, 
                            reshape(profiles[:, id_LCT[4]], (nb_rows, length(id_LCT[4])))
                            )
    end

    return LCT_profiles
end

