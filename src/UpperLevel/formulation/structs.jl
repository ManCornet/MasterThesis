# -- CONFIGURATION STRUCTURE --
abstract type PowerFlowFormulation end      # can be either Jabr or distflow
abstract type TypeofProdFormulation end     # No_DG or DG generation (boolean) => replace by a "if" 
abstract type RadialityFormulation end      # simple - single commodity flow - multi-commodity - spanning tree constraints
abstract type TopologyChoiceFormulation end      # simple - single commodity flow - multi-commodity - spanning tree constraints
abstract type TypeOfGraph end
abstract type ConvexityFormulation end
abstract type VoltageConstraintsFormulation end
abstract type CurrentConstraintsFormulation end


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

struct Undirected <: TypeOfGraph end
struct Directed <: TypeOfGraph end

struct NonConvex <: ConvexityFormulation end
struct Convex <: ConvexityFormulation end

struct RelaxedVoltages <: VoltageConstraintsFormulation end
struct StrongVoltages <: VoltageConstraintsFormulation end

struct RelaxedCurrents <: CurrentConstraintsFormulation end
struct StrongCurrents <: CurrentConstraintsFormulation end

struct Formulation
    powerflow::PowerFlowFormulation 
    production::TypeofProdFormulation
    radiality::RadialityFormulation
    topology_choice::TopologyChoiceFormulation
    graph_type::TypeOfGraph
    convexity::ConvexityFormulation
    v_constraints::VoltageConstraintsFormulation
    i_constraints::CurrentConstraintsFormulation


    function Formulation(;
        powerflow::PowerFlowFormulation = BFM(),
        production::TypeofProdFormulation = NoDG(),
        radiality::RadialityFormulation = SingleCommodityFlow(),
        topology_choice::TopologyChoiceFormulation = OneConfig(),
        graph_type::TypeOfGraph = Undirected(),
        convexity::ConvexityFormulation = Convex(),
        v_constraints::VoltageConstraintsFormulation = StrongVoltages(),
        i_constraints::CurrentConstraintsFormulation = RelaxedCurrents(),
        )                     
        return new(powerflow, production, radiality, topology_choice, graph_type, 
                    convexity, v_constraints, i_constraints)
    end
end