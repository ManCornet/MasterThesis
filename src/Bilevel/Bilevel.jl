module Bilevel

import StructTypes, JSON3
import DataFrames
import Random 
import XLSX
using JuMP, Gurobi, BilevelJuMP
using Logging, Printf
import Plots
using Graphs, GraphRecipes
using PrettyTables, Latexify, LaTeXStrings
#using PGFPlotsX, 
using Statistics
ENV["CPLEX_STUDIO_BINARIES"] = "/Applications/CPLEX_Studio221/cplex/bin/x86-64_osx/"
import Pkg
Pkg.add("CPLEX")
Pkg.build("CPLEX")
import CPLEX

export DSOCosts, UserCosts, Simulation, Formulation, build_model 
export get_nb_loads, get_nb_substations, get_nb_conductors
export get_nb_buses, get_nb_lines, get_nb_nodes, get_nb_time_steps
export get_network_data, define_pu_basis
export build_daily_PV_profiles, build_daily_load_profiles, build_profiles
export add_load_profiles!, add_PV_profiles!, add_storage!, process_time_steps
export save_struct
export Network, NetworkTopology, Storage 
export print_load_profiles, print_PV_profiles, print_network_tikz



#export Network
#export Bus
include("formulation/structs.jl")
include("structs.jl")
include("read_network_data.jl")
include("profiles.jl")
#export Formulation
include("formulation/variables.jl")
include("formulation/constraints.jl")
include("formulation/objective.jl")
include("build_model.jl")
include("plot_results.jl")


end