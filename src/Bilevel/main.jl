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

using ArgParse
using Plots
include("structs.jl")

# =============================================================================
#                                   Functions
# =============================================================================
# Structure Definition
# =============================================================================
#                         Definition of structures
# =============================================================================
# Function fielnames to get the fields of a structure

# Main argument parser

function parse_commandline()

    s = ArgParseSettings()

    @add_arg_table s begin
        "--EV"
            help = "The model to use"
            default = false
            arg_type = Bool
            
        "--HP"
            help = "The model to use"
            default = false
            arg_type = Bool
        
        "--PV_CAPA"
            help = "Maximum PV capacity per load bus"
            arg_type = Float64
            default = 0.4

        "--IMP_ELECTRICITY_ENRG_COST"
            help = "Cost of the energy that is imported"
            arg_type = Float64
            default = 0.3

        "--EXP_ELECTRICITY_ENRG_COST"
            help = "Cost of the energy that is exported"
            arg_type = Float64
            default = 0.1

        "--IMP_ELECTRICITY_DSO_COST"
            help = "DSO cost of the energy that is imported/exported"
            arg_type = Float64
            default = 0.1

        "--GRID_CONNECTION_COST"
            help = "DSO cost of the energy that is imported/exported"
            arg_type = Float64
            default = 80
    
    end
    return parse_args(s)
end


println()