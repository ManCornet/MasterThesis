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
import Plots
using LaTeXStrings

"""
    Idea of this part:
    - Creating the profiles, (think about ways to add additionnal significative days)
    - Think about ways to make the code modular to add additionnal significative days
    - Create PV profiles
    - When it is done add profiles etc to the structures
    - Then compute the model + tests 

"""

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

    # ---- Checking arguments ----
    nrows, ncols = size(time_serie)
    1 <= nb_agg_periods <= nrows || throw(
    ArgumentError("[change_granularity] nb_agg_periods =  $nb_agg_periods not in [1, $nrows]"))

    # ---- Changing granularity ----
    lost_time_steps = nrows % nb_agg_periods 
    new_time_serie  =   [   sum(time_serie[i:(i+nb_agg_periods-1), j] / nb_agg_periods) 
                            for i in 1:nb_agg_periods:(nrows - lost_time_steps), 
                                j in 1:ncols
                        ]
    
    return new_time_serie, lost_time_steps
end

function check_profiles_dim(profiles::Vector{Matrix{Float64}})
    sizes = [size(p) for p in profiles]
    return all(y->y==first(sizes), sizes)
end

function build_profiles(daily_profiles::Vector{Matrix{Float64}};
                        scaling_factor::Float64=1.0)

    # -- Checking argument : all profiles must have same dimensions --
    check_profiles_dim(daily_profiles) || throw(
        ArgumentError("[build_load_profiles] the daily profiles don't have the same dimensions"))

    # -- Merge all daily load profiles and find peak values --
    profiles = vcat(daily_profiles...) * scaling_factor
    peak_value, peak_time_step = findmax(vec(sum(profiles, dims=2)))

    return profiles, peak_value, peak_time_step
end

function build_daily_load_profiles( PROFILE_PATH::String, 
                                    nb_profiles::Int64;
                                    nb_agg_periods::Int64=1,
                                    EV::Bool=false,
                                    EHP::Bool=false,
                                    EV_PATH::Union{Nothing,String}=nothing,
                                    EHP_PATH::Union{Nothing,String}=nothing,
                                    scaling_EV::Float64=1.0,
                                    scaling_EHP::Float64=1.0,
                                    seed::Union{Nothing,Int64}=nothing
    )


    # --- Fetching the load profiles data ---
    load_profiles = convert(Matrix{Float64}, XLSX.readxlsx(PROFILE_PATH)[1][:])
    nrows, ncols = size(load_profiles)

    # ---- Checking arguments ----
    1 <= nb_profiles <= ncols || throw(
    ArgumentError("""   [build_daily_load_profiles]: 
                        nb_profiles not in [1, $ncols]: nb_profiles = $nb_profiles"""))

    1 <= nb_agg_periods <= nrows || throw(
    ArgumentError("""   [build_daily_load_profiles]: 
                        number of aggregated periods not in [1, $nrows]: nb_agg_periods = $period"""))

    # ---- Building load profiles ----
    if !isnothing(seed)
        Random.seed!(seed)
        id_loads = Random.rand(1:ncols, nb_profiles)
    else 
        id_loads = 1:nb_profiles
    end
    # select only a subpart of the profiles
    load_profiles = reshape(load_profiles[:, id_loads], (nrows, nb_profiles))

    # ---- Adding EV profiles to consumption ----
    if EV
        @assert !isnothing(scaling_EV)

        EV_profiles = convert(Matrix{Float64}, XLSX.readxlsx(EV_PATH)[1][:])
        EV_profiles = reshape(load_profiles[:, id_loads], (nrows, nb_profiles))
        load_profiles .+= scaling_EV * EV_profiles
    end 

    # ---- Adding EHP profiles to consumption ----
    if EHP
        @assert !isnothing(scaling_EHP)

        EHP_profiles = convert(Matrix{Float64}, XLSX.readxlsx(EHP_PATH)[1][:])
        EHP_profiles = reshape(load_profiles[:, id_loads], (nrows, nb_profiles))
        load_profiles .+= scaling_EHP * EHP_profiles 
    end

    # ---- Aggregating profiles to change granularity ----
    lost_time_steps = 0
    if nb_agg_periods > 1 
        load_profiles, lost_time_steps = change_granularity(load_profiles; nb_agg_periods=nb_agg_periods)
    end
    return load_profiles, lost_time_steps
end


# -- FUNCTIONS TO CREATE THE PV PROFILES --
function build_daily_PV_profiles(   PV_PROFILE_PATH::String,
                                    nb_profiles::Integer; 
                                    scaling_PV::Float64=1.0,
                                    nb_agg_periods::Integer=1,
                                    seed::Union{Nothing,Int64}=nothing
                                )

    # --- Fetching the PV profiles data ---
    PV_profiles = convert(Matrix{Float64}, XLSX.readxlsx(PV_PROFILE_PATH)[1][:])
    nrows, ncols = size(PV_profiles)

    # ---- Checking arguments ----
    1 <= nb_profiles <= ncols || throw(
    ArgumentError("""   [build_daily_PV_profiles]: 
                        nb_profiles not in [1, $ncols]: nb_profiles = $nb_profiles"""))

    1 <= nb_agg_periods <= nrows || throw(
    ArgumentError("""   [build_daily_PV_profiles]: 
                        number of aggregated periods not in [1, $nrows]: nb_agg_periods = $period"""))

    # ---- Building PV profiles ----
    if !isnothing(seed)
        Random.seed!(seed)
        id_profiles = Random.rand(1:ncols, nb_profiles)
    else 
        id_profiles = 1:nb_profiles
    end

    # select only a subpart of the profiles
    PV_profiles = reshape(PV_profiles[:, id_profiles], (nrows, nb_profiles))

    lost_time_steps = 0
    if nb_agg_periods > 1 
        PV_profiles, lost_time_steps = change_granularity(PV_profiles;nb_agg_periods=nb_agg_periods)
    end

    return (PV_profiles ./ maximum(PV_profiles, dims=1)) * scaling_PV, lost_time_steps, id_profiles
end


function print_load_profiles(fig_name::String,
                        base_profiles::Matrix{Float64},
                        profiles::Matrix{Float64}; 
                        base_granularity::Int64=5,
                        delta_t::Int64=5,
                        EV::Bool=false,
                        EHP::Bool=false)

    colors = ["#bb5b46", "#b37e6e", "#bc9780", "#dfa878", "#5f8d8d", "#3f7e84", "#096f7b"]
    base_time_steps = vec(1:size(base_profiles)[1]) * base_granularity
    time_steps = vec(1:size(profiles)[1]) * delta_t
    label = "Base load"
    fig = Plots.plot(base_time_steps, sum(base_profiles, dims=2), label=label, xlabel = "Time [min.]", ylabel="Power Consumption [MVA]", linewidth=1, formatter=:latex, color=colors[1], linealpha=0.8)
    label = "Aggregated load"
    label = EV ? (label * " + EV") : label 
    label = EHP ? (label * " + EHP") : label 
    Plots.plot!(fig, time_steps, sum(profiles, dims=2), label=label,linewidth=1.5, formatter=:latex, tickfontsize=10, color=colors[7])
    Plots.savefig(fig, fig_name)
    return
end

function print_PV_profiles( fig_name::String,
                            profiles::Matrix{Float64};
                            delta_t::Integer=5,
                            id_profiles)

    colors = ["#bb5b46", "#b37e6e", "#bc9780", "#dfa878", "#5f8d8d", "#3f7e84", "#096f7b"]
    time = vec(1:size(profiles)[1]) * delta_t
    # label = vec(["PV installation of user $i" for i in id_profiles])
    # fig = Plots.plot(time, profiles, label=hcat(id_profiles...), xlabel = "Time [min.]$", ylabel="Power Production [% Peak Power]", linewidth=1.5, formatter=:latex)
    fig = Plots.plot(time, profiles[:, 1], label="", xlabel = "Time [min.]", ylabel="Power Production [% Peak Power]", linewidth=1.5, formatter=:latex, tickfontsize=10, color=colors[1])
    Plots.savefig(fig, fig_name)
    return
end