# Choice of the formulation
abstract type PowerFlowFormulation end      # can be either Jabr or distflow
abstract type TypeofProdFormulation end     # No_DG or DG generation 
abstract type NetworkGraphFormulation end   # directed or undirected
abstract type RadialityFormulation end      # simple - single commodity flow - multi-commodity - spanning tree constraints
#=  Idea: two additional config that can be added:
    ---------------------------------------------
    abstract type TimeFormulation end               => time dependent formulation or time independent 
    abstract type SubstationCapacityFormulation end => Jabr or my formulation 
=#

# Inspired from: 
# https://github.com/ANL-CEEESA/UnitCommitment.jl/blob/dev/src/model/formulations/base/structs.jl
"""
    struct Formulation
        powerflow::PowerFlowFormulation
    end

    Struct provided to `build_model` that holds various formulation components.
# Fields
- `powerflow`: Formulation for the production decision variables
"""


# Creating a module for each formulation

struct Formulation
    powerflow::PowerFlowFormulation 
    production::TypeofProdFormulation
    radiality::RadialityFormulation
    networkgraph::NetworkGraphFormulation
    

    function Formulation(;
        powerflow::PowerFlowFormulation = Jabr2012.PowerFlow(), # 2 formulations for now
        production::TypeofProdFormulation = ,
        radiality::RadialityFormulation = ,                     # choice among 5 radiality formulations
        networkgraph::NetworkGraphType = )                      # only two either directed or undirected formulation


        return new(powerflow, production, radiality, networkgraph)
    end
end
