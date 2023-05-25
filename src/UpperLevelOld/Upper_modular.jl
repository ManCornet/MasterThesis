# -- IMPORTATION OF THE PACKAGES --
import XLSX, DataFrames
using JuMP, Gurobi, Printf, Dates
include("../structs.jl")
include("../read_network_data.jl")
include("../profiles.jl")
include("../utils.jl")
#include("../Bilevel/parameters_BFM_1P.jl")

# ---------------------------------------------------------------------------- #
#                           Structure for the formulation                      #
# ---------------------------------------------------------------------------- #

# -- CONFIGURATION STRUCTURE --
abstract type PowerFlowFormulation end      # can be either Jabr or distflow
abstract type TypeofProdFormulation end     # No_DG or DG generation (boolean) => replace by a "if" 
abstract type RadialityFormulation end      # simple - single commodity flow - multi-commodity - spanning tree constraints

struct BIM <: PowerFlowFormulation end
struct BFM <: PowerFlowFormulation end

struct NoDG <: TypeofProdFormulation end
struct DG <: TypeofProdFormulation end

struct SimpleRadiality <: RadialityFormulation end
struct SingleCommodityFlow <: RadialityFormulation end
struct MultiCommodityFlow <: RadialityFormulation end
struct SpanningTree <: RadialityFormulation end

struct Formulation
    powerflow::PowerFlowFormulation 
    production::TypeofProdFormulation
    radiality::RadialityFormulation
    networkgraph::NetworkGraphFormulation
    condvars::CondVarsFormulation

    function Formulation(;
        powerflow::PowerFlowFormulation = BFM(),
        production::TypeofProdFormulation = NoDG(),
        radiality::RadialityFormulation = SimpleRadiality(),  
        )                     
        return new(powerflow, production, radiality)
    end
end



# ---------------------------------------------------------------------------- #
#                           Function to add pieces of code                     #
# ---------------------------------------------------------------------------- #
 
function _add_BusVariables!(model::JuMP.AbstractModel, ::NoDG)::Nothing 

    N   = model[:network_data][:bus][1]
    Ns  = model[:network_data][:sub_bus][1]
    MIN_VOLTAGE = model[:network_data][:bus][4]
    MAX_VOLTAGE = model[:network_data][:bus][5]

    @variables( model,   
                begin 
                MIN_VOLTAGE^2 <= V_sqr[N] <= MAX_VOLTAGE^2 , (container=Array)
                P_sub[N] >= 0                                , (container=Array)
                Q_sub[N]                                     , (container=Array)
                S_sub[N]                                     , (container=Array)
                beta[Ns]                                   , (container=Array, binary=true)    
                end
            )
            
    for i in Nu
        fix(P_G[i], 0.0; force=true) 
        fix(Q_G[i], 0.0)
        fix(S_G[i], 0.0)
    end
    
    return
end

function _add_BusVariables!(model, ::DG)::Nothing 

    N   = model[:network_data][:bus][1]
    Ns  = model[:network_data][:sub_bus][1]
    MIN_VOLTAGE = model[:network_data][:bus][4]
    MAX_VOLTAGE = model[:network_data][:bus][5]

    @variables( model,   
                begin 
                MIN_VOLTAGE^2 <= V_sqr[N] <= MAX_VOLTAGE^2 , (container=Array)
                P_G[N] >= 0                                , (container=Array)
                Q_G[N]                                     , (container=Array)
                S_G[N]                                     , (container=Array)
                beta[Ns]                                   , (container=Array, binary=true)    
                end
            )
    return
end

# -- BRANCH VARIABLES --
# Here maybe think to aggregate or not in variables, another option 
function _add_BranchVariables!(model, ::DirectedGraph)::Nothing 
    # Here maybe try to test another formulations without variables
    L = model[:network_data][:line][1]
    K = model[:network_data][:conductor]

    JuMP.@variables( model,   
                begin 
                P_ij_k[L, K]        , (container=Array) 
                Q_ij_k[L, K]        , (container=Array)
                I_sqr_k[L, K] >= 0  , (container=Array)
                end
    )

    return
end

function _add_BranchVariables!(model, ::UnDirectedGraph, ::NonAggrCondVars)::Nothing 
    # Here maybe try to test another formulations without variables
    L   = model[:network_data][:line][1]
    K   = model[:network_data][:conductor]

    JuMP.@variables( model,   
                begin 
                P_ij_k[L, K]        , (container=Array) 
                P_ji_k[L, K]        , (container=Array) 
                Q_ij_k[L, K]        , (container=Array)
                Q_ji_k[L, K]        , (container=Array)
                I_sqr_k[L, K] >= 0  , (container=Array)
                end
    )
    return
end

function _add_BranchVariables!(model, ::DirectedGraph, ::AggrCondVars)::Nothing 
    # Here maybe try to test another formulations without variables
    L   = model[:network_data][:line][1]
    K   = model[:network_data][:conductor]

    JuMP.@variables( model,   
                begin 
                P_ij_k[L, K]        , (container=Array) 
                Q_ij_k[L, K]        , (container=Array)
                I_sqr_k[L, K] >= 0  , (container=Array)
                P_ij[L]             , (container=Array) 
                Q_ij[L]             , (container=Array) 
                I_sqr[L] >= 0       , (container=Array)
                end
    )

    return
end

function _add_BranchVariables!(model, ::UnDirectedGraph, ::AggrCondVars)::Nothing 
    # Here maybe try to test another formulations without variables
    L   = model[:network_data][:line][1]
    K   = model[:network_data][:conductor]

    JuMP.@variables( model,   
                begin 
                P_ij_k[L, K]        , (container=Array) 
                P_ji_k[L, K]        , (container=Array) 
                Q_ij_k[L, K]        , (container=Array)
                Q_ji_k[L, K]        , (container=Array)
                I_sqr_k[L, K] >= 0  , (container=Array)
                P_ij[L]             , (container=Array)
                P_ji[L]             , (container=Array) 
                Q_ij[L]             , (container=Array)
                Q_ji[L]             , (container=Array)  
                I_sqr[L] >= 0       , (container=Array)
                end
    )

    return
end

# -- REFERENCE VOLTAGE AT SUBSTATIONS -- 

function _add_RefVoltages!(model)::Nothing

    Ns      = model[:network_data][:sub_bus][1]
    Ns_init = model[:network_data][:sub_bus][2]
    Ns_notinit = setdiff(Ns, Ns_init)
    MIN_VOLTAGE = model[:network_data][:bus][4]
    MAX_VOLTAGE = model[:network_data][:bus][5]

    JuMP.@constraint(model, 
                ref_voltage_sub_init[i=Ns_init], 
                model[:V_sqr][i] == 1
    )
    
    JuMP.@constraint(model, 
                ref_voltage_sub_notinit1[i=Ns_notinit], 
                model[:V_sqr][i] - 1 <= (MAX_VOLTAGE^2 - 1)*(1 - model[:beta][i])
    )

    JuMP.@constraint(model, 
                ref_voltage_sub_notinit2[i=Ns_notinit], 
                model[:V_sqr][i] - 1 >= (MIN_VOLTAGE^2 - 1)*(1 - model[:beta][i])
    ) 
    
    return 
end


# -- LOAD OVERSATISFACTION CONSTRAINTS --

function _add_LoadOverSatisfaction!(model)::Nothing

    N   = model[:network_data][:bus][1]
    P_D = model[:network_data][:load_bus][2]
    Q_D = model[:network_data][:load_bus][3]
    
    JuMP.@constraint(model, 
                    LoadOverSatisfaction_active, 
                    sum(model[:P_G][i] - P_D[i] for i in N) >= 0
    )

    JuMP.@constraint(model, 
                    LoadOverSatisfaction_reactive, 
                    sum(model[:Q_G][i] - Q_D[i] for i in N) >= 0
    )

    return 
end

# Add the functionaility to design for 1 significative days or a snapshot
# 
# Function that builds the model
function build_model(network_data::Network,
                    DSO_costs::DSOCosts;
                    formulation = Formulation()
                    )::JuMP.Model

    @info "Building model..."
    time_model = @elapsed begin
        model = JuMP.Model(Gurobi.Optimizer)

        # ====================== Set up the Gurobi solver =====================
        set_optimizer_attribute(model, "TimeLimit", 200)
        set_optimizer_attribute(model, "Presolve", 0)

        model[:network_data] = network_data
        model[:DSO_costs]    = DSO_costs

        # =========================== Build the model =========================
        # -- Add the variables of the model --
        _add_BusVariables!(model, formulation.production)
        _add_BranchVariables!(model, )

        #_add_BranchVariables!(model, formulation.networkgraph, formulation.condvars)
        _add_RefVoltages!(model)
        _add_LoadOverSatisfaction!(model)
        #_add_PowerFlowEqs!(model, formulation.powerflow, formulation.networkgraph, formulation.condvars)
            
    end

    @info @sprintf("Built model in %.2f seconds", time_model)


    return model
end