# -- CONFIGURATION STRUCTURE --
abstract type PowerFlowFormulation end      # can be either Jabr or distflow
abstract type TypeofProdFormulation end     # No_DG or DG generation (boolean) => replace by a "if" 
abstract type RadialityFormulation end      # simple - single commodity flow - multi-commodity - spanning tree constraints
abstract type TopologyChoiceFormulation end      # simple - single commodity flow - multi-commodity - spanning tree constraints


struct BIM <: PowerFlowFormulation end
struct BFM <: PowerFlowFormulation end

struct NoDG <: TypeofProdFormulation end
struct DG <: TypeofProdFormulation end

struct SimpleRadiality <: RadialityFormulation end
struct SingleCommodityFlow <: RadialityFormulation end
struct MultiCommodityFlow <: RadialityFormulation end
struct SpanningTree <: RadialityFormulation end

struct ReconfigAllowed <: TopologyChoiceFormulation end
struct OneConfig <: TopologyChoiceFormulation end


struct Formulation
    powerflow::PowerFlowFormulation 
    production::TypeofProdFormulation
    radiality::RadialityFormulation
    choice_topology::TopologyChoiceFormulation
    
    function Formulation(;
        powerflow::PowerFlowFormulation = BFM(),
        production::TypeofProdFormulation = NoDG(),
        radiality::RadialityFormulation = SimpleRadiality(),
        choice_topology::TopologyChoiceFormulation = OneConfig()
        )                     
        return new(powerflow, production, radiality, choice_topology)
    end
end