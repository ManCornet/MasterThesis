#-----------------------------------------------------------------------------
#
#           - Planning & Operation of electric power and energy systems - 
#                   Homework : Implementation of the DNEP 
#
#                             University of Liege
#
#-----------------------------------------------------------------------------
# Created By  : Geoffrey Bailly, Manon Cornet
# Created Date: Tuesday November 21 2023
#
# main:
#   Main file
#
# =============================================================================
#                                   Imports
# =============================================================================

import XLSX
import DataFrames
using JuMP, Gurobi
using ArgParse, Printf


include("models.jl")

# =============================================================================
#                         Definition of global constants
# =============================================================================

const global HOURS_IN_A_YEAR = 8760                       
const global BASE_VOLTAGE    = 15                                 # [kV]
const global BASE_POWER      = 1e3                                # [kVA]
const global BASE_CURRENT    = BASE_POWER / BASE_VOLTAGE          # [A]
const global BASE_ADMITTANCE = 1e-3 * BASE_CURRENT / BASE_VOLTAGE # [S]
const global MIN_VOLTAGE     = 0.95                               # [pu]
const global MAX_VOLTAGE     = 1.05                               # [pu]
const global CONDUCTOR_COST_PER_MM2_PER_KM = 200                  # [€/mm^2/km]

# =============================================================================
#                                   Functions
# =============================================================================

# Main argument parser
function parse_commandline()

    s = ArgParseSettings()

    @add_arg_table s begin
        "--formulation"
            help = "The formulation to use"
            default = 1
            arg_type = Int
    end
    return parse_args(s)
end

# Prints a basic section title in terminal
function section(title)

    # Number of letters to determine section size
    # title_size = length(title)
    boundary  = "="
    for i in 1:(30)
        boundary *= "="
    end
    
    # Printing section
    @printf("\n \n")
    @printf("%s", boundary)
    @printf(" %s ", title)
    @printf("%s", boundary)
    @printf("\n \n")
end

function print_table(table, table_name, unity, complex)

    end_string = 12 - length(table_name)
    @printf("%s", table_name)
    for i in 1:end_string
        @printf(" ")
    end 
    @printf(": [ ")
    
    for i in table
        if complex
            @printf("%.4f + i %.4f", real(i), imag(i)) 
        else
            @printf("%.4f ", i) 
        end
    end
    @printf("]")
    @printf("    [%s]\n", unity)
end

function process_conductors(df_conductor, line_length, max_current, 
                            conductance, susceptance, line_cost, 
                            conductor_idx, line_idx
                            )

    max_i = df_conductor.max_i_ka[conductor_idx]*1e3
    max_current[conductor_idx][line_idx] = max_i / BASE_CURRENT
    
    r = line_length[line_idx] * df_conductor.r_ohm_per_km[conductor_idx]
    x = line_length[line_idx] * df_conductor.x_ohm_per_km[conductor_idx]
    y = 1/(r+im*x) / BASE_ADMITTANCE
    
    conductance[conductor_idx][line_idx] = real(y) 
    susceptance[conductor_idx][line_idx] = imag(y)

    section = df_conductor.q_mm2[conductor_idx]
    line_cost[conductor_idx][line_idx] = (section * line_length[line_idx] 
                                           *  CONDUCTOR_COST_PER_MM2_PER_KM)

    return
end

# =============================================================================
#                                   Main
# =============================================================================

function main()

    # ================= Parsing arguments of main command line ================
    parsed_args = parse_commandline()
    args = [arg for  (arg, val) in parsed_args]
    vals = [val for (arg, val) in parsed_args]

    formulation = vals[1]
   
    println("Execution of the formulation $formulation")

    # ======================== DN parameter definitions =======================

    # ---- Path of the file containing the DN topology ----

    XLSX_FILE_PATH = "model_2S2H.xlsx"

    # ---- Line parameters ----

    df_line     = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "line"))
    L           = size(df_line)[1]     # number of lines
    line_length = df_line.length_km    # line lengths [km]
    line_ends   = Dict(l => (df_line.from_bus[l], df_line.to_bus[l]) for l in 1:L)

    # ---- Bus parameters ----

    df_bus = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "bus"))
    N      = size(df_bus)[1]    # number of buses 

    # ---- Link btw lines and nodes ----

    Omega_sending   = Dict(n => [] for n in 1:N)
    Omega_receiving = Dict(n => [] for n in 1:N)
    for l in 1:L
        push!(Omega_sending[line_ends[l][1]], l)
        push!(Omega_receiving[line_ends[l][2]], l)
    end

    # ---- Line properties ---- 

    df_conductor = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "line_std_types"))

    K = 3                                  # number of conductor types

    max_current = Dict(k => [0.0 for l in 1:L] for k in 1:K) # absolute, [pu]
    conductance = Dict(k => [0.0 for l in 1:L] for k in 1:K) # absolute, [pu]
    susceptance = Dict(k => [0.0 for l in 1:L] for k in 1:K) # absolute, [pu]
    line_cost   = Dict(k => [0.0 for l in 1:L] for k in 1:K) # [€/km]
    
    for k in 1:K
        for l in 1:L
            process_conductors(df_conductor, line_length, max_current, 
                               conductance, susceptance, line_cost, 
                               k, l)
        end
    end

    # ---- Substation parameters ----

    n_s_init = 0     # Should be <= n_s
    n_s      = 2     # Susbstation nodes are numbered to correspond to the first n_s nodes in the network
    substation_fixed_cost = 1e5 * [10, 1]   # Fixed construction or reinforcement cost of substations [€] 
    substation_op_cost    = ones(n_s)   # Substations operation cost [€/kVAh^2]   
    S_rating_init = 0 * ones(n_s) / BASE_POWER
    S_rating_max  = 200 * ones(n_s) / BASE_POWER

    # ---- Demand at buses ----  

    df_load  = DataFrames.DataFrame(XLSX.readtable(XLSX_FILE_PATH, "load"))
    P_demand = [zeros(n_s); df_load.p_mw * 1e3 / BASE_POWER]
    Q_demand = [zeros(n_s); df_load.q_mvar * 1e3 / BASE_POWER]

    # ---- Objective function related parameters ----  

    K_l = 1     # Capital recovery rate of line constructions
    K_s = 1     # Capital recovery rate of substation construction or reinforcement
    
    line_loss                = 0.01     # phi_l : loss factor of lines
    cost_unit_loss           = 1        # c_l   : loss factor of lines
    substation_utilization   = 0.01     # phi_s : cost per energy lost [€/kWh]
    interest_rate_losses     = 0.02     # tau_l : interest rate for the cost of power losses
    interest_rate_substation = 0.02     # tau_s : interest rate for the substation operation cost

    # ======================== MINLP formulation =======================
    if formulation == 1
        section("Formulation 1: MINLP")

        I_squared, V_re, V_im, x, alpha, beta,
        P_s, Q_s, P_conductor, Q_conductor, 
        obj, time = MINLP_formulation(  N, n_s, L, K ,
                                        K_l, K_s, 
                                        substation_utilization,
                                        S_rating_init, S_rating_max,
                                        conductance, susceptance, max_current, 
                                        line_cost, line_ends, line_length, line_loss, 
                                        interest_rate_losses, interest_rate_substation,
                                        substation_op_cost, substation_fixed_cost,cost_unit_loss,  
                                        Omega_sending, Omega_receiving, 
                                        P_demand, Q_demand
                                    )


        # Refaire les print
        @printf("-------------------------------------------------------------\n")
        @printf("-------------------------- Results --------------------------\n")
        @printf("-------------------------------------------------------------\n")

        alpha_int = [int(i) for i in alpha]
        beta_int = [int(i) for i in beta]
        x_int = [int(i) for i in x]
        
        I = [sqrt(i) for i in I_squared]
        V = [(v_real + im*v_imag)*BASE_VOLTAGE for (v_real, v_imag) in zip(V_re, V_im)]
        P_sub = [p*BASE_POWER for p in P_s]
        Q_sub = [q*BASE_POWER for q in Q_s]
        P_conductor = [p*BASE_POWER for p in P_conductor]
        Q_conductor = [q*BASE_POWER for q in Q_conductor]

        @printf("-> Binary variables: \n")
        print_table(alpha_int, "Alpha" , "-", false) # 2D dimensions
        print_table(x_int, "x", "-", false)
        print_table(beta_int, "Beta", "-", false)
        @printf("-> Continuous variables: \n")
        print_table(I, "I" , "A", false)
        print_table(V, "V", "kV", true)
        print_table(P_s, "P_s", "kW", false)
        print_table(Q_s, "Q_s", "kVar", false)
        print_table(P_conductor, "P_conductor", "kW", false)
        print_table(Q_conductor, "Q_conductor", "kVar", false)

       

    # ======================== MISOCP formulation =======================
    elseif formulation == 2
        section("Formulation 2: MISOCP")

        I_squared, V_re, V_im, x, alpha, beta,
        P_s, Q_s, P_conductor, Q_conductor, 
        obj, time = MISOCP_formulation(  N, n_s, L, K ,
                                        K_l, K_s, 
                                        substation_utilization,
                                        S_rating_init, S_rating_max,
                                        conductance, susceptance, max_current, 
                                        line_cost, line_ends, line_length, line_loss, 
                                        interest_rate_losses, interest_rate_substation,
                                        substation_op_cost, substation_fixed_cost,cost_unit_loss,  
                                        Omega_sending, Omega_receiving, 
                                        P_demand, Q_demand
                                    )

    
    

    end
end

main()