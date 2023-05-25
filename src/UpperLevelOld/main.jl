#-----------------------------------------------------------------------------
#
#                           - TFE : Bilevel DNEP - 
#                             University of Liege
#
#-----------------------------------------------------------------------------
# Created By  : Manon Cornet
# Created Date: Saturday March 18 2023
#
# new_main:
#   Main file
#
# =============================================================================
#                                   Imports
# =============================================================================
import XLSX
import DataFrames
using Dates

using JuMP, Gurobi
using ArgParse

# Formulations
include("upper_level_Jabr.jl")
include("DNEP_time_formulation.jl")

# Parameters
include("post_process.jl")
include("parameters.jl")


# =============================================================================
#                                   Functions
# =============================================================================

# Main argument parser
function parse_commandline()

    s = ArgParseSettings()

    @add_arg_table s begin
        "--formulation"
            help = "The formulation to use"
            default = "Jabr"
            arg_type = String
    end
    return parse_args(s)
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
    
    # ================ Models ==================
    if formulation == "Jabr"
        sets  = (N, Ns, K, L, Omega_sending, Omega_receiving)
        costs = (substation_fixed_cost, substation_op_cost, line_cost, cost_unit_loss)
        substation_param = (S_rating_init, S_rating_max, substation_utilization, interest_rate_substation)
        conductor_param  = (max_current, conductance, susceptance)
        losses           = (line_loss, interest_rate_losses)
        others           = K_s, K_l
        demand           = (P_demand, Q_demand)
    
        var_values, var_sets = UL_Jabr(sets, costs, substation_param, conductor_param, losses, others, demand)

        XLSX_PATH  = "output.xlsx"
        # period covered
        for (key, value) in var_values
            if ndims(value) <= 1 || key == "I_squared" || key == "alpha"
                add_var_to_XLSX(XLSX_PATH, value, key, var_sets[key])
            elseif ndims(value) == 2
                processed_var = process2D_variable(value)
                add_var_to_XLSX(XLSX_PATH, processed_var, key, var_sets[key])
            elseif ndims(value) == 3
                X_i_ij = process_X_i_ij(value)
                add_var_to_XLSX(XLSX_PATH, X_i_ij, "X_i_ij", var_sets[key])
            end
        end

    elseif formulation == "time"
        sets  = (N, Ns, Ns_init, K, L, L_init, T, Omega_sending, Omega_receiving)
        costs = (sub_expan_cost, sub_install_cost, line_cost, losses_cost, DSO_INTEREST_RATE)
        substation_param = (S_rating_init, S_rating_max)
        conductor_param  = (max_current, conductance, susceptance)
        demand           = (load_profiles, tan_phi)
      
        var_values, var_sets = time_dependent_formulation(sets, costs, substation_param, conductor_param, 
                                                          demand, delta_t)

        start_date = DateTime(today() + Day(1))
        end_date   = DateTime(start_date + Day(1))
        XLSX_PATH  = "output.xlsx"
        # period covered
        date_range = start_date:Minute(GRANULARITY*period):(end_date-Second(1))
        for (key, value) in var_values
            if ndims(value) <= 2
                add_var_to_XLSX(XLSX_PATH, value, key, var_sets[key], date_range = date_range)
        
            elseif ndims(value) == 3
                processed_var = process3D_variable(value)
                add_var_to_XLSX(XLSX_PATH, processed_var, key, var_sets[key], date_range = date_range)
        
            elseif ndims(value) == 4
                X_i_ij = process_X_i_ij_time(value)
                add_var_to_XLSX(XLSX_PATH, X_i_ij, "X_i_ij", var_sets[key], date_range = date_range)
            end
        end
    end
end

main()