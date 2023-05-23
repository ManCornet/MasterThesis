#-----------------------------------------------------------------------------
#
#                           - TFE : Bilevel DNEP - 
#                             University of Liege
#
#-----------------------------------------------------------------------------
# Created By  : Manon Cornet
# Created Date: Saturday May 20 2023
#
# new_main:
#   Main file
#
# =============================================================================
#                                   Imports
# =============================================================================
# Activating the julia environement
# Path: I must add in terminal julia --project -- src/main.jl --

using ArgParse, PrettyTables, Printf
using Plots, JSON3
using StructTypes

include("structs.jl")
include("read_network_data.jl")
include("profiles.jl")

# =============================================================================
#                                   Functions
# =============================================================================

function print_title(title::String)
    println(title)
    for _ in 1:length(title) @printf("-") end
    @printf("\n\n")
end

function print_segment()
    println("#-------------------------------------------------------------------------------------------")
end

function print_header()
    println("""#-------------------------------------------------------------------------------------------
#
#                                       Bilevel DNEP
#
#                                      Graduation work
#
#-------------------------------------------------------------------------------------------
# @ Manon Cornet
\n""")
end

function save(R::DistributionNetworkTopology, filename::String)
    f = open(filename,"w")
    JSON3.write(f, R)
    close(f)
end
StructTypes.StructType(::Type{DistributionNetworkTopology}) = StructTypes.Mutable()

function parse_commandline(;as_symbols::Bool=false)

    s = ArgParseSettings()

    @add_arg_table s begin
        "--model"
            help = """model = - bilevel if simulation is run with bilevel model
                              - single if simulation is run with single level model 
                    """
            arg_type = String
            default = "single"
        
        "--delta_t"
            help = """Granularity of the simulation in minutes 
                    [Requirements]: must be a multiple of the granularity of the 
                                    dataset containing the load & PV profiles 
                    """
            arg_type = Number
            default = 5

        "--EV"
            help = "EVs taken into account in residential consumption profiles"
            arg_type = Bool
            default = false
           
        "--EHP"
            help = "EHPs taken into account in residential consumption profiles"
            arg_type = Bool
            default = false
            
        "--PV_CAPA"
            help = "Maximum PV capacity per load bus"
            arg_type = Number
            default = 0.4

        "--IEC"
            help = "Imported electricity cost"
            arg_type = Number
            default = 0.3

        "--EEC"
            help = "Exported electricity cost"
            arg_type = Number
            default = 0.1

        "--DSOEC"
            help = "DSO electricity cost"
            arg_type = Number
            default = 0.1

        "--GCC"
            help = "Grid connection cost"
            arg_type = Number
            default = 80

        "--PVC"
            help = "PV cost"
            arg_type = Number
            default = 500
        
        "--PP"
            help = "Load peak power"
            arg_type = Number
            default = 7
        
        "--verbose"
            help = "Print the details of the simulation"
            arg_type = Bool
            default = true
    end
    return parse_args(s; as_symbols=as_symbols)
end


# =============================================================================
#                                Main function
# =============================================================================


function main()
    # ================= Parsing arguments of main command line ================
    parsed_args = parse_commandline(as_symbols=false)
    model = parsed_args["model"]
    delta_t = parsed_args["delta_t"]
    EV = parsed_args["EV"]
    EHP = parsed_args["EHP"]
    PV_CAPA = parsed_args["PV_CAPA"]
    IEC = parsed_args["IEC"]
    EEC = parsed_args["EEC"]
    DSOEC = parsed_args["DSOEC"]
    GCC = parsed_args["GCC"]
    verbose = parsed_args["verbose"]
    PV_cost = parsed_args["PVC"]
    peak_power = parsed_args["PP"]

    # ================= Printing the parameters of the simulation =============
    if verbose 
        print_header()
        print_title("Running a simulation with the following characteristics:")
        header = (  ["Model Type", "EV", "EHP","Delta_t", "PV_CAPA", "IEC", 
                    "EEC", "DSOEC", "GCC"],
                    ["[-]", "[-]","[-]", "[min]", "[MVA]", "[€/kWh]", 
                    "[€/kWh]", "[€/kWh]", "[€/kVA/y]"])

        data = [model EV EHP delta_t PV_CAPA IEC EEC DSOEC GCC]
        pretty_table(   data; header = header, 
                        header_crayon = crayon"yellow bold", tf = tf_unicode_rounded)
        print_segment()
    end

    # =========================== Network data  ===============================

    # -- Fetching the path of the root directory --
    root_dir = normpath(joinpath(@__FILE__,"..","..",".."))

    # -- Loading the excel file containing the network topology --
    network_data_dir = joinpath(root_dir, "NetworkModels")
    NETWORK_PATH = joinpath(network_data_dir, "network_Nahman_Peric_2S23H.xlsx") 

    # -- Getting the network data --
    network, network_topology = get_network_data(NETWORK_PATH)
    # print_network_topology(network_topology)
    # save(network_topology, "network_topology.json")
    # save(network_data, "network_data.json")

    # =========================== Load profiles  ==============================

    # -- Loading the excel file containing the data for the load profiles --
    profiles_data_dir = joinpath(root_dir, "ManchesterData", "LCT_profiles")
    SUMMER_LOAD_PATH  = joinpath(profiles_data_dir, "Summer_Load_Profiles.xlsx")
    WINTER_LOAD_PATH  = joinpath(profiles_data_dir, "Winter_Load_Profiles.xlsx")
    EV_PATH = EV ? joinpath(profiles_data_dir, "Winter_EV_Profiles.xlsx") : nothing
    EHP_PATH = EHP ? joinpath(profiles_data_dir, "Winter_EHP_Profiles.xlsx") : nothing
  
    # -- Building the base load profile on which to scale --
    PROFILE_PATHS = [SUMMER_LOAD_PATH, WINTER_LOAD_PATH]
    base_daily_profiles = [ build_daily_load_profiles(path, get_nb_load_bus(network)) * 1e-3 / network.pu_basis.base_power 
                            for path in PROFILE_PATHS]
    
    base_load_profiles, base_peak_value, base_peak_time_step = build_load_profiles(base_daily_profiles)
    scaling_factor = peak_power / base_peak_value
    base_load_profiles *= scaling_factor

    # -- Building the load profiles with desired characteristics --
    nb_agg_time_steps = process_time_steps(delta_t=delta_t, data_granularity=5)

    summer_profiles = build_daily_load_profiles( SUMMER_LOAD_PATH, 
                                                get_nb_load_bus(network);
                                                nb_agg_periods=nb_agg_time_steps,
                                                EV=EV, 
                                                EHP=EHP,
                                                EV_PATH=EV_PATH, 
                                                EHP_PATH=EHP_PATH,
                                                scaling_EHP=0.0) * 1e-3 / network.pu_basis.base_power
        
    winter_profiles = build_daily_load_profiles( WINTER_LOAD_PATH, 
                                                get_nb_load_bus(network);
                                                nb_agg_periods=nb_agg_time_steps,
                                                EV=EV, 
                                                EHP=EHP,
                                                EV_PATH=EV_PATH, 
                                                EHP_PATH=EHP_PATH) * 1e-3 / network.pu_basis.base_power

    load_profiles, peak_value, peak_time_step = build_load_profiles([summer_profiles, winter_profiles]; 
                                                                    scaling_factor=scaling_factor)
                   
    
    print_profiles("test1.pdf",base_load_profiles, load_profiles, delta_t, EV=EV, EHP=EHP)

    add_load_profiles!(network, load_profiles, delta_t=delta_t)

    base_load_profiles = hcat([b.load_profile.time_serie for b in network.load_buses]...)
    print_profiles("test2.pdf",base_load_profiles, load_profiles, delta_t, EV=EV, EHP=EHP)
    
    # -- Creating the PV profiles -- 
    #XLSX_PV_PATH = joinpath(profiles_data_dir, "Summer_PV_Profiles.xlsx")

    # -- Adding them to the network data structure -- 

    # -- Run a simulation of the model --
    # bilevel_model(network_data, price_dict)
    # this call to the function automatically updates the structures 

end

main()
