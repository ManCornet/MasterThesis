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

using JuMP, Gurobi
using ArgParse, Printf

include("building_profiles.jl")
include("parameters.jl")
include("upper_level_models.jl")

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
    
    # ================ Paths of the load and PV profiles files ==================

    if formulation == "Jabr"
        
    elseif formulation == "time"
        sets  = (N, Ns, Ns_init, K, L, L_init, Y, T, Omega_sending, Omega_receiving)
        costs = (sub_expan_cost, sub_install_cost, line_cost, losses_cost, DSO_INTEREST_RATE)
        substation_param = (S_rating_init, S_rating_max)
        conductor_param  = (max_current, conductance, susceptance)
        demand           = (yearly_load_profiles, tan_phi)
     
        results = time_dependent_formulation(sets, costs, substation_param, conductor_param, 
                                             demand, delta_t)

        write_results(results)

    elseif formulation == "radiality"
    
    end

end

main()