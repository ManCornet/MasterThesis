module UpperLevel

import StructTypes, JSON3
import DataFrames
import Random 
import XLSX
using JuMP, Gurobi
using Logging, Printf
using Plots
using Graphs, GraphRecipes
using PrettyTables, Latexify

export DSOCosts, UserCosts, Simulation, Formulation, build_model 
export get_nb_loads, get_nb_substations, get_nb_conductors
export get_nb_buses, get_nb_lines, get_nb_nodes, get_nb_time_steps
export get_network_data, define_pu_basis
export build_daily_PV_profiles, build_daily_load_profiles, build_profiles
export add_load_profiles!, add_PV_profiles!, process_time_steps
export build_model
export save_struct
export Network, NetworkTopology 
export print_load_profiles, print_PV_profiles, print_network_tikz


include("structs.jl")
include("read_network_data.jl")
include("profiles.jl")
#export Network
#export Bus
include("formulation/structs.jl")
#export Formulation
include("formulation/variables.jl")
include("formulation/constraints.jl")
include("formulation/objective.jl")
include("build_model.jl")
include("plot_results.jl")


# -- Here I test the functionalities --


end