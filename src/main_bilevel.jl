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
#                                   Functions
# =============================================================================
include("./Bilevel/Bilevel.jl")
include("utils.jl")
using .Bilevel
using ArgParse, Dates
using PrettyTables, XLSX, DataFrames

function parse_commandline(;as_symbols::Bool=false)

    s = ArgParseSettings()

    @add_arg_table s begin

        "--verbose"
            help = "Print the details of the simulation"
            arg_type = Bool
            default = true

        "--network_graph_name"
            help = "Name of the network graphs"
            arg_type = String
            default = "network_simu"

        "--plot_file_name"
            help = "Name of the network graphs"
            arg_type = String
            default = "plot_simu"
        
        "--simu_name"
            help = "Simulation name"
            arg_type = String
            default = "simu"

        # ------- Choice of the model -------
        "--bilevel"
            help = """model = - bilevel if bilevel = true
                              - singlelevel if bilevel = false
                    """
            arg_type = Bool
            default = false

        "--storage"
            help = "If storage or not in the model"
            arg_type = Bool
            default = false

        "--network_reconfig"
            help = "If network can be reconfigurated"
            arg_type = Bool
            default = false

        # ------- Profiles parameters -------
        "--days"
            help = "Number of significative days [day]"
            arg_type = Int64
            default = 2
        
        "--delta_t"
            help = """Granularity of the simulation [min.] 
                    [Requirements]: must be a multiple of the granularity of the 
                                    dataset containing the load & PV profiles 
                    """
            arg_type = Int64
            default = 60 # 1h per default
        
        "--PP"
            help = "Load peak power"
            arg_type = Float64
            default = 7.0 # 7.0 CHANGE THAT AFTER

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

        "--Storage_pen"
            help = "Storage penetration in [0, 1] is the % of users having PV)"
            arg_type = Float64
            default = 1.0

        "--Storage_eff"
            help = "Storage efficiency in [0, 1]"
            arg_type = Float64
            default = 0.9

        "--Storage_cost"
            help = "Storage cost"
            arg_type = Float64
            default = 500.0

        "--AMORT_STORAGE"
            help = "Storage amortization"
            arg_type = Int64
            default = 15
        
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
            default = 1.0e3

        "--AMORT_DSO"
            help = "Duration of the amortization of the DSO [year]"
            arg_type = Int64
            default = 50

        "--IR_DSO"
            help = "Interest rate DSO [%/year]"
            arg_type = Float64
            default = 0.06

        "--weight_I"
            help = "Weight of the penalization I [k€/MVA]"
            arg_type = Float64
            default = 1.0e-2

        "--weight_V"
            help = "Weight of the penalization V [k€/MVA]"
            arg_type = Float64
            default = 1.0e-2
        
        "--relax_voltage"
            help = "Weight of the penalization V [k€/MVA]"
            arg_type = Bool
            default = false
    end
    return parse_args(s; as_symbols=as_symbols)
end


# =============================================================================
#                                Main function
# =============================================================================


function main()

    # ================= Parsing arguments of main command line ================
    parsed_args  = parse_commandline(as_symbols=false)
    param_keys   = collect(keys(parsed_args))
    param_values = collect(values(parsed_args))

    verbose     = parsed_args["verbose"]

    # ------- Choice of the model -------
    bilevel     = parsed_args["bilevel"]
    storage     = parsed_args["storage"]

    # ------- Profiles parameters -------
    nb_days     = parsed_args["days"]
    delta_t     = parsed_args["delta_t"]
    peak_power  = parsed_args["PP"]
    EV          = parsed_args["EV"]
    EHP         = parsed_args["EHP"]

    # ------- Users parameters -------
    PV_pen        = parsed_args["PV_pen"]
    storage_pen   = parsed_args["Storage_pen"]
    storage_eff   = parsed_args["Storage_eff"]
    storage_cost  = parsed_args["Storage_cost"]
    amort_storage = parsed_args["AMORT_STORAGE"]
    PV_capa       = parsed_args["PV_CAPA"]
    PV_cost       = parsed_args["PVC"]
    PV_conv_cost  = parsed_args["PVCC"]
    amort_PV_conv = parsed_args["AMORT_PV_CONV"]
    amort_PV      = parsed_args["AMORT_PV"]
    EIC           = parsed_args["EIC"] #old 0.3
    EEC           = parsed_args["EEC"]
    DSOEC         = parsed_args["DSOEC"]
    GCC           = parsed_args["GCC"]
    cos_phi       = 0.95

    # ------- DSO parameters -------
    network_reconfig    = parsed_args["network_reconfig"]
    substation_cost     = parsed_args["SUB_COST"]
    amort_DSO           = parsed_args["AMORT_DSO"]
    interest_rate_DSO   = parsed_args["IR_DSO"]
    weight_I = weight_V = parsed_args["weight_I"]
    money_basis = 1.0
    weight_obj1 = weight_obj2 = 1.0

    # ================= Printing the parameters of the simulation =============
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
    root_dir          = splitdir(@__DIR__)[1]
    simulations_dir   = joinpath(root_dir, "simulations")
    network_data_dir  = joinpath(root_dir, "NetworkModels")
    profiles_data_dir = joinpath(root_dir, "ManchesterData", "LCT_profiles")

    # -- Loading the excel file containing the network topology --

    # Add choice for the test network
    NETWORK_PATH = joinpath(network_data_dir, "network_Nahman_Peric_2S23H.xlsx") 
    #NETWORK_PATH = joinpath(network_data_dir, "model_2S2H.xlsx") 
    pu_basis = define_pu_basis()

    # -- Fetching the network data --
    network, network_topology = Bilevel.get_network_data(NETWORK_PATH; max_pv_capa=PV_capa, pu_basis=pu_basis, cos_phi=cos_phi)
    #print_network_topology(network_topology)


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
        base_daily_profile, _ = build_daily_load_profiles(path, get_nb_loads(network))
        push!(base_daily_profiles, base_daily_profile * 1e-3 / network.pu_basis.base_power)
    end

    base_load_profiles, base_peak, _ = build_profiles(base_daily_profiles)
    scaling_factor = peak_power / base_peak
    base_load_profiles *= scaling_factor

    # -- Building the load profiles with desired characteristics --
    nb_agg_periods = process_time_steps(delta_t=delta_t, data_granularity=5)

    scaling_EHP = [0.0, 1.0]
    scaling_EV  = [1.0, 1.0]

    @assert length(scaling_EHP) == nb_days
    @assert length(scaling_EV) == nb_days

    daily_profiles = Vector{Matrix{Float64}}()
    for (index, path) in enumerate(PROFILE_PATHS)
        daily_profile, _ = build_daily_load_profiles(  path, 
                                                    get_nb_loads(network);
                                                    nb_agg_periods=nb_agg_periods,
                                                    EV=EV, 
                                                    EHP=EHP,
                                                    EV_PATH=EV_PATH, 
                                                    EHP_PATH=EHP_PATH,
                                                    scaling_EHP=scaling_EHP[index],
                                                    scaling_EV=scaling_EV[index])
        push!(daily_profiles, daily_profile * 1e-3 / network.pu_basis.base_power)
    end

    load_profiles, _ , _ = build_profiles(daily_profiles; scaling_factor=scaling_factor)



    # -- Add load profiles to network structure -- 
    add_load_profiles!(network, load_profiles; delta_t=delta_t, pu_basis=pu_basis)
    #save_struct(network, "network_data.json")

    # =========================== PV profiles  ==============================

    if PV_pen > 0

        PV_PATH = joinpath(profiles_data_dir, "Summer_PV_Profiles.xlsx")
        scaling_PV = [1.0, 0.1]
        nb_PV_profiles = floor(Int, PV_pen * get_nb_loads(network))
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


        # -- Add PV profiles to network structure -- 
        PQ_diagram = (max_q = 0.3, slope=-1.0)
        add_PV_profiles!(network, PV_profiles, ids_profiles[1]; 
                        PQ_diagram = PQ_diagram, delta_t=delta_t)

        # -- Add storage -- 
        if storage 
            seed = nothing
        
            nb_storage_units = floor(Int, storage_pen * nb_PV_profiles)
            if !isnothing(seed)
                Random.seed!(seed)
                id_storage = Random.rand(ids_profiles[1], nb_storage_units)
            else 
                id_storage = ids_profiles[1][1:nb_storage_units]
            end

            add_storage!(network, storage_eff, id_storage)
        end
    end

    # =========================== Costs definition  ===========================
                               
    DSO_costs = DSOCosts(   
                            substation_cost, 
                            0.7 * EIC, 
                            amort_DSO, 
                            interest_rate_DSO, 
                            weight_I, 
                            weight_V, 
                            money_basis, 
                            weight_obj1
                        )

    User_costs = UserCosts(
                            PV_cost, 
                            PV_conv_cost, 
                            storage_cost, 
                            EIC, 
                            EEC, 
                            DSOEC, 
                            DSOEC, 
                            GCC, 
                            amort_PV, 
                            amort_PV_conv, 
                            amort_storage, 
                            money_basis, 
                            weight_obj2
                        )

    # =========================== Creating the simulation  ===========================

    if network_reconfig 
        topology_choice = Bilevel.ReconfigAllowed()
    else 
        topology_choice = Bilevel.OneConfig()
    end

    formulation = Bilevel.Formulation(  
                                        powerflow = Bilevel.BFM(),
                                        topology_choice = topology_choice,
                                        graph_type = Bilevel.Undirected(),
                                        radiality = Bilevel.MultiCommodityFlow(),
                                        convexity = Bilevel.Convex(),
                                        v_constraints = Bilevel.StrongVoltages(),
                                        i_constraints = Bilevel.RelaxedCurrents()
                                    )

    nb_sign_days = length(PROFILE_PATHS)
    simulation   = Bilevel.Simulation(
                                        network, 
                                        network_topology, 
                                        DSO_costs, 
                                        User_costs, 
                                        nb_sign_days, 
                                        bilevel, 
                                        storage, 
                                        formulation)

    # =========================== Solving the model  ===========================
    simu_path = joinpath(simulations_dir, parsed_args["simu_name"])
    !isdir(simu_path) && mkdir(simu_path)

    # --- Solving input data --
    input_data_path = joinpath(simu_path, "input_data")
    !isdir(input_data_path) && mkdir(input_data_path)

    print_load_profiles( joinpath(input_data_path, "load_profiles.pdf"), 
                        base_load_profiles, 
                        load_profiles;
                        base_granularity=5, 
                        delta_t=delta_t, 
                        EV=EV, 
                        EHP=EHP
                    )

    print_PV_profiles(joinpath( input_data_path, "PV_profiles.pdf"), PV_profiles;
    delta_t=delta_t, id_profiles=ids_profiles[1])

    # --- Solving model --
    failure = false
    result = nothing
    try
        model = Bilevel.build_model(simulation; set_names=true)
        return_value = Bilevel.solve_model(model, formulation.powerflow)
        result = Bilevel.printed_tables(model)
    catch e
        print(e)
        failure = true
        result =  ["failure", ""]
    end

    param_table = Vector()
    push!(param_table, param_keys)
    push!(param_table, param_values)


    XLSX_PATH = joinpath(simu_path, "results.xlsx")
    isfile(XLSX_PATH) ? mode = "rw" : mode = "w"

    XLSX.openxlsx(XLSX_PATH, mode = mode) do xf
        #mode == "rw" && XLSX.addsheet!(xf)
        #sheet_names = XLSX.sheetnames(xf)
        #current_sheet = sheet_names[end]
        sheet_names = XLSX.sheetnames(xf)
        current_sheet = sheet_names[1]
        current_working_sheet = xf[current_sheet]
        XLSX.writetable!(current_working_sheet, param_table, ["", ""])
        df = DataFrame(result[2:end, :], result[1, :])
        df = permutedims(df, 1)
        XLSX.writetable!(current_working_sheet, df, anchor_cell=XLSX.CellRef("D1"))
    end

    if !failure
        # -- SAVE JSON FILE STRUCTURES --
        json_path = joinpath(simu_path, "json")
        !isdir(json_path) && mkdir(json_path)

        save_struct(network_topology, joinpath(json_path, "network_topology.json"))
        save_struct(network, joinpath(json_path, "network_data.json"))

        # -- PRINTING THE NETWORK AT EACH TIME STEP --
        network_plot_path = joinpath(simu_path, "network_plots")
        !isdir(network_plot_path) && mkdir(network_plot_path)
       
        for t in 1:model[:time_steps]
            Bilevel.print_network_tikz(model[:network_data], t, 5.5, 5; 
            dir=network_plot_path, filename=parsed_args["network_graph_name"] * "_timestep_$t", display=false,reshape=true)
        end

        # -- DIFFERENT PLOTS RESULTING FROM A SIMULATION --
        figures_path = joinpath(simu_path, "figures")
        !isdir(figures_path) && mkdir(figures_path)

        Bilevel.plot_results(model; dir=figures_path, filename=parsed_args["plot_file_name"], pgfplot=true)
    end

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
