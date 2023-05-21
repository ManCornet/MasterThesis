#-----------------------------------------------------------------------------
#
#                           - TFE : Bilevel DNEP - 
#                             University of Liege
#
#-----------------------------------------------------------------------------
# Created By  : Manon Cornet
# Created Date: Saturday March 20 2023
#
# profiles:
#   File containing the functions required to build load and PV profiles
#
# =============================================================================
#                                   Imports
# =============================================================================
import Random 
import XLSX 

""" process_time_steps 

    Arguments:
    ----------
        - delta_t           : size of time series time step
        - data_granularity  : granularity of the dataset
    
    Return value:
    -------------
        - nb_agg_time_steps : nb of time steps to aggregate
"""
function process_time_steps(;delta_t::Integer=5, data_granularity::Integer=5)
    # Ensuring that the time_step size is greater or equal than the granularity
    @assert delta_t >= data_granularity
    @assert delta_t % data_granularity == 0

    nb_agg_time_steps = delta_t ÷ data_granularity

    return nb_agg_time_steps 
end


# -- FUNCTION FOR AGGREGATING THE LOAD PROFILES TO HAVE TO HAVE LESS TIME STEPS --

function change_granularity(time_serie::AbstractArray; nb_agg_periods::Integer = 1)
    # We suppose the time series is of size (time, load_profiles)

    # ---- Checking arguments ----
    nb_periods = size(time_serie)[1]
    1 <= nb_agg_periods <= nb_periods || throw(
    ArgumentError("[change_granularity] nb_agg_periods =  $nb_agg_periods not in [1, $nb_periods]"))

    # ---- Changing granularity ----
    nb_rows, nb_columns = size(time_serie)
    new_time_serie = Matrix{Float64}(undef, 0, nb_columns)

    mod = nb_rows % nb_agg_periods 
    for i in 1:nb_agg_periods:(nb_rows - mod)
        new_time_serie = [new_time_serie; sum(time_serie[i:(i+nb_agg_periods-1), :], dims=1)]
    end
    
    return new_time_serie ./ nb_agg_periods
end

# String with significative days
# -- FUNCTIONS TO CREATE THE LOAD PROFILES --
# Function that creates the base load profiles
function create_load_profiles(  PROFILE_DIR::String,
                                nb_profiles::Int64, 
                                peak_power::Float64,
                                peak_distribution::Vector{Float64};
                                winter::Bool=true, 
                                EV::Bool=false, 
                                EHP::Bool=false,
                                seed=nothing
                                )

    

    # ---- Fetch base load profiles ---
    SUMMER_LOAD_PATH = joinpath(PROFILE_DIR, "Summer_Load_Profiles.xlsx")
    base_load_profiles = convert(Matrix{Float64}, XLSX.readxlsx(SUMMER_LOAD_PATH)[1][:])
    summer_load_profiles = deepcopy(base_load_profiles)
    nrows, ncols = size(base_load_profiles)
    if EV 
        EV_LOAD_PATH = joinpath(PROFILE_DIR, "Winter_EV_Profiles.xlsx")
        base_load_EV = convert(Matrix{Float64}, XLSX.readxlsx(EV_LOAD_PATH)[1][:])
        summer_load_profiles .+= base_load_EV
    end

    if EHP 
        EHP_LOAD_PATH = joinpath(PROFILE_DIR, "Winter_EHP_Profiles.xlsx")
        base_load_EHP = convert(Matrix{Float64}, XLSX.readxlsx(EHP_LOAD_PATH)[1][:])
        base_load_profiles .+= base_load_EHP
    end

    if winter 
        WINTER_LOAD_PATH = joinpath(PROFILE_DIR, "Winter_Load_Profiles.xlsx")
        base_winter_profiles = convert(Matrix{Float64}, XLSX.readxlsx(WINTER_LOAD_PATH)[1][:])
        nrows_winter, ncols_winter = size(base_winter_profiles)

        @assert nrows == nrows_winter
        @assert ncols == ncols_winter

        base_load_profiles = vcat(base_load_profiles, base_winter_profiles)
    end

    # Fixing the seed to choose randomly the load profiles selected 
    if !isnothing(seed)
        Random.seed!(seed)
        id_loads = Random.rand(1:ncols, nb_profiles)
    else 
        id_loads = 1:nb_profiles
    end

    # Base load profiles for the users
    base_load_profiles = reshape(base_load_profiles[:, id_loads], (nrows, length(id_loads)))

    # ---- Scale load profiles ---


     

    # ---- Checking arguments ----
    nb_min_day = 24 * 60
    data_granularity = nb_min_day ÷ nrows
    nb_agg_time_steps = process_time_steps(delta_t=delta_t, data_granularity=data_granularity)

    1 <= nb_agg_time_steps <= nrows || throw(
    ArgumentError("[create_load_profiles] number of aggregated periods not in [1, $nrows]"))

    1 <= nb_profiles <= ncols || throw(
    ArgumentError("[create_load_profiles] nb_profiles not in [1, $ncols]"))

    # ---- Rescale according to peak_power before aggregating ----
    
    if nb_agg_time_steps > 1

    end
return




# -- FUNCTIONS TO CREATE THE PV PROFILES --
function create_PV_profiles(   PROFILE_DIR::String; 
    winter::Bool=true, 
    EV::Bool=false, 
    HP::Bool=false,
    peak_power::Float64=7.0,
    delta_t::Int64=60
    )


return