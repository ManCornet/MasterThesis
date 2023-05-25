#-----------------------------------------------------------------------------
#
#                           - TFE : Bilevel DNEP - 
#                             University of Liege
#
#-----------------------------------------------------------------------------
# Created By  : Manon Cornet
# Created Date: Saturday May 20 2023
#
# main:
#   Main file
#
# =============================================================================
#                                   Imports
# =============================================================================
# Activating the julia environement
# Path: I must add in terminal julia --project -- src/main.jl --

using ArgParse, PrettyTables
using Plots, JSON3
using StructTypes

include("structs.jl")
include("read_network_data.jl")
include("profiles.jl")
include("utils.jl")

# =============================================================================
#                                   Functions
# =============================================================================
function parse_commandline(;as_symbols::Bool=false)

    s = ArgParseSettings()

    @add_arg_table s begin

        "--verbose"
            help = "Print the details of the simulation"
            arg_type = Bool
            default = true

        # ------- Choice of the model -------
        "--bilevel"
            help = """model = - bilevel if bilevel = true
                              - singlelevel if bilevel = false
                    """
            arg_type = Bool
            default = false

        # ------- Profiles parameters -------
        "--days"
            help = "Number of significative days [day]"
            arg_type = Integer
            default = 2
        
        "--delta_t"
            help = """Granularity of the simulation [min.] 
                    [Requirements]: must be a multiple of the granularity of the 
                                    dataset containing the load & PV profiles 
                    """
            arg_type = Int64
            default = 5
        
        "--PP"
            help = "Load peak power"
            arg_type = Float64
            default = 7.0

        "--EV"
            help = "EVs taken into account in residential consumption profiles [-]"
            arg_type = Bool
            default = false
           
        "--EHP"
            help = "EHPs taken into account in residential consumption profiles [-]"
            arg_type = Bool
            default = false
        
        # ------- Users parameters -------
        "--PV_pen"
            help = "Pv penetration in [0, 1] is the % of users having PV)"
            arg_type = Float64
            default = 1.0
        
        "--PV_CAPA"
            help = "Maximum PV capacity per load bus [MVA]"
            arg_type = Float64
            default = 0.4

        "--PVC"
            help = "PV cost [€/kWp]"
            arg_type = Float64
            default = 500.0

        "--PVCC"
            help = "PV converter cost [€/kVA]"
            arg_type = Float64
            default = 200.0
       
        "--AMORT_PV_CONV"
            help = "Duration of the amortization period of a PV converter [year]"
            arg_type = Int64
            default = 10
        
        "--AMORT_PV"
            help = "Duration of the amortization period of a PV installation [year]"
            arg_type = Int64
            default = 25

        "--EIC"
            help = "Imported electricity cost [€/kWh]"
            arg_type = Float64
            default = 0.3

        "--EEC"
            help = "Exported electricity cost [€/kWh]"
            arg_type = Float64
            default = 0.1

        "--DSOEC"
            help = "DSO electricity cost [€/kWh]"
            arg_type = Float64
            default = 0.1

        "--GCC"
            help = "Grid connection cost [€/kVA/y]"
            arg_type = Float64
            default = 80.0

        # ------- DSO parameters -------

        "--SUB_COST"
            help = "Substation cost [k€/MVA]"
            arg_type = Float64
            default = 1e3

        "--AMORT_DSO"
            help = "Duration of the amortization of the DSO [year]"
            arg_type = Int64
            default = 50

        "--IR_DSO"
            help = "Interest rate DSO [%/year]"
            arg_type = Int64
            default = 50
    end
    return parse_args(s; as_symbols=as_symbols)
end


# =============================================================================
#                                Main function
# =============================================================================


function main()
    # ================= Parsing arguments of main command line ================
    parsed_args = parse_commandline(as_symbols=false)
    verbose     = parsed_args["verbose"]

    # ------- Choice of the model -------
    bilevel     = parsed_args["bilevel"]

    # ------- Profiles parameters -------
    nb_days     = parsed_args["days"]
    delta_t     = parsed_args["delta_t"]
    peak_power  = parsed_args["PP"]
    EV          = parsed_args["EV"]
    EHP         = parsed_args["EHP"]

    # ------- Users parameters -------
    PV_pen        = parsed_args["PV_pen"]
    PV_capa       = parsed_args["PV_CAPA"]
    PV_cost       = parsed_args["PVC"]
    PV_conv_cost  = parsed_args["PVCC"]
    amort_PV_conv = parsed_args["AMORT_PV_CONV"]
    amort_PV      = parsed_args["AMORT_PV"]
    EIC           = parsed_args["EIC"]
    EEC           = parsed_args["EEC"]
    DSOEC         = parsed_args["DSOEC"]
    GCC           = parsed_args["GCC"]

    # ------- DSO parameters -------
    substation_cost   = parsed_args["SUB_COST"]
    amort_DSO         = parsed_args["AMORT_DSO"]
    interest_rate_DSO = parsed_args["IR_DSO"]
    
    # ================= Printing the parameters of the simulation =============
    # Idea do the tables for each type of parameters !!!
    if verbose 
        print_header()
        print_title("Running a simulation with the following characteristics:")
        header = (  ["Model Type", "EV", "EHP","Delta_t", "PV_CAPA", "EIC", 
                    "EEC", "DSOEC", "GCC"],
                    ["[-]", "[-]","[-]", "[min]", "[MVA]", "[€/kWh]", 
                    "[€/kWh]", "[€/kWh]", "[€/kVA/y]"])
        model = bilevel ? "bilevel" : "singlelevel"
        data = [model EV EHP delta_t PV_capa EIC EEC DSOEC GCC]
        pretty_table(   data; header = header, 
                        header_crayon = crayon"yellow bold", tf = tf_unicode_rounded)
        print_segment("-")
    end

    # =========================== Network data  ===============================

    # -- Fetching the path of the main directories --
    root_dir = normpath(joinpath(@__FILE__,"..","..",".."))
    plot_dir = joinpath(root_dir, "plots")
    network_data_dir  = joinpath(root_dir, "NetworkModels")
    profiles_data_dir = joinpath(root_dir, "ManchesterData", "LCT_profiles")

    # -- Loading the excel file containing the network topology --
    NETWORK_PATH = joinpath(network_data_dir, "network_Nahman_Peric_2S23H.xlsx") 

    # -- Fetching the network data --
    network, network_topology = get_network_data(NETWORK_PATH)
    # print_network_topology(network_topology)
    # save(network_topology, "network_topology.json")
    # save(network_data, "network_data.json")

    # =========================== Load profiles  ==============================

    # -- Loading the excel file containing the data for the load profiles --
    SUMMER_LOAD_PATH = joinpath(profiles_data_dir, "Summer_Load_Profiles.xlsx")
    WINTER_LOAD_PATH = joinpath(profiles_data_dir, "Winter_Load_Profiles.xlsx")
    EV_PATH  = EV ? joinpath(profiles_data_dir, "Winter_EV_Profiles.xlsx") : nothing
    EHP_PATH = EHP ? joinpath(profiles_data_dir, "Winter_EHP_Profiles.xlsx") : nothing
  
    # -- Building the base load profile on which to scale --
    PROFILE_PATHS = [SUMMER_LOAD_PATH, WINTER_LOAD_PATH]

    base_daily_profiles = Vector{Matrix{Float64}}()
    for path in PROFILE_PATHS 
        base_daily_profile, _ = build_daily_load_profiles(path, get_nb_load_bus(network))
        push!(base_daily_profiles, base_daily_profile * 1e-3 / network.pu_basis.base_power)
    end
  
    base_load_profiles, base_peak, _ = build_profiles(base_daily_profiles)
    scaling_factor = peak_power / base_peak
    base_load_profiles *= scaling_factor

    # -- Building the load profiles with desired characteristics --
    nb_agg_periods = process_time_steps(delta_t=delta_t, data_granularity=5)
    println(nb_agg_periods)
    scaling_EHP = [0.0, 1.0]
    scaling_EV  = [1.0, 1.0]

    @assert length(scaling_EHP) == nb_days
    @assert length(scaling_EV) == nb_days

    daily_profiles = Vector{Matrix{Float64}}()
    for (index, path) in enumerate(PROFILE_PATHS)
        daily_profile, _ = build_daily_load_profiles(  path, 
                                                    get_nb_load_bus(network);
                                                    nb_agg_periods=nb_agg_periods,
                                                    EV=EV, 
                                                    EHP=EHP,
                                                    EV_PATH=EV_PATH, 
                                                    EHP_PATH=EHP_PATH,
                                                    scaling_EHP=scaling_EHP[index],
                                                    scaling_EV=scaling_EV[index])
        println(typeof(daily_profile))
        push!(daily_profiles, daily_profile * 1e-3 / network.pu_basis.base_power)
    end

    load_profiles, _ , _ = build_profiles(daily_profiles; scaling_factor=scaling_factor)
    
    print_load_profiles( joinpath(plot_dir, "load_profiles.pdf"), 
                        base_load_profiles, 
                        load_profiles;
                        base_granularity=5, 
                        delta_t=delta_t, 
                        EV=EV, 
                        EHP=EHP
                    )

    # -- Add load profiles to network structure -- 
    add_load_profiles!(network, load_profiles; delta_t=delta_t)

    # =========================== PV profiles  ==============================
    PV_PATH = joinpath(profiles_data_dir, "Summer_PV_Profiles.xlsx")
    scaling_PV = [1.0, 0.1]
    nb_PV_profiles = floor(Int, PV_pen * get_nb_load_bus(network))
    @assert length(scaling_PV) == nb_days

    daily_PV_profiles = Vector{Matrix{Float64}}()
    ids_profiles = []
    for d in 1:nb_days
        daily_PV_profile, _, id_profiles = build_daily_PV_profiles(PV_PATH,
                                                                nb_PV_profiles; 
                                                                scaling_PV=scaling_PV[d],
                                                                nb_agg_periods=nb_agg_periods,
                                                                seed=nothing)
        push!(daily_PV_profiles, daily_PV_profile)
        push!(ids_profiles, id_profiles)
    end

    PV_profiles, _ , _ = build_profiles(daily_PV_profiles)

    print_PV_profiles(joinpath( plot_dir, "PV_profiles.pdf"), PV_profiles;
                                delta_t=delta_t, id_profiles=ids_profiles[1])

    # -- Add PV profiles to network structure -- 
    PQ_diagram = (max_q = 0.3, slope=-1.0)
    add_PV_profiles!(network, PV_profiles, ids_profiles[1]; 
                    capa_max = PV_capa, PQ_diagram = PQ_diagram, delta_t=delta_t)

    # =========================== Costs definition  ===========================
    money_basis = 1.0

    DSO_costs  = DSOCosts(substation_cost, 0.7 * EIC, amort_DSO, interest_rate_DSO, money_basis)
    User_costs = UserCosts(PV_cost, PV_conv_cost, EIC, EEC, DSOEC, DSOEC, GCC, amort_PV, amort_PV_conv, money_basis)

    # =========================== Model  ===========================
    # -- Running the model -- 
    
    bilevel!(network_data, DSO_costs, User_costs; BILEVEL=bilevel, RADIALITY=true)


    # -- Run a simulation of the model --
    # 
    # this call to the function automatically updates the structures 
    #load_profiles = hcat([b.load_profile.time_serie for b in network.load_buses]...)
    #print_load_profiles(joinpath( plot_dir, "test.pdf"), 
    #                    base_load_profiles, 
    #                    load_profiles;
    #                    base_granularity=5, 
    #                    delta_t=delta_t, 
    #                    EV=EV, 
    #                    EHP=EHP)
    
    #PV_profiles = hcat([b.PV_installation.profile.time_serie for b in network.load_buses]...)
    #print_PV_profiles(joinpath( plot_dir, "PV_profiles.pdf"), PV_profiles;
    #                    delta_t=delta_t, id_profiles=ids_profiles[1])
    # Test for PV profiles 
end

main()
