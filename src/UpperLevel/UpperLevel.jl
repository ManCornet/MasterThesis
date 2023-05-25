module UpperLevel

import JuMP, Gurobi
import StructTypes
using Logging, Printf

include("structs.jl")
#export COORD, VLIM, PU_BASIS, Node, Bus
#export Node 
#export Network
#export Bus
include("formulation/structs.jl")
export Formulation
include("formulation/variables.jl")
include("formulation/constraints.jl")
include("formulation/objective.jl")
include("build_model.jl")


# -- Here I test the functionalities --


end