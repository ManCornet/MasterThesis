#-----------------------------------------------------------------------------
#
#                           - TFE : Bilevel DNEP - 
#                             University of Liege
#
#-----------------------------------------------------------------------------
# Created By  : Manon Cornet
# Created Date: Wednesday May 24 2023
#
# bilevel_model:
#   file containing the bilevel model
#
# =============================================================================
#                                   Imports
# =============================================================================
# Activating the julia environement
# Path: I must add in terminal julia --project -- src/main.jl --

using JuMP, BilevelJuMP, Gurobi
using GraphRecipes, Graphs, Plots, Dates, DataFrames
if (!@isdefined LAUNCH_SENSITIVITY_ANALYSIS) || (LAUNCH_SENSITIVITY_ANALYSIS != true)
    include("parameters_BFM_1P.jl")
end
include("export_xlsx.jl")

# =============================================================================
#                                   Functions
# =============================================================================
# Several options within this function: storage or not 

# In the upper model there are several possibilities:
# - Test the radiality constraints to use: 4 to test
# - If the configuration of the network can change at each time 
#   step
function _add_upper_variables!(Upper::JuMP.AbstractModel)::Nothing
    @variables(Upper, begin
        MIN_VOLTAGE^2 <= V_sqr[T, N] <= MAX_VOLTAGE^2, (container = Array)
        I_sqr_k[T, L, K] >= 0, (container = Array)
        I_sqr[T, L] >= 0, (container = Array)
        P_sub[T, Ns], (container = Array)
        Q_sub[T, Ns], (container = Array)
        S_sub[T, Ns] >= 0, (container = Array) # Because apparent power
        P_ij_k[T, L, K], (container = Array)
        Q_ij_k[T, L, K], (container = Array)
        P_ij[T, L], (container = Array)
        Q_ij[T, L], (container = Array)
        S_sub_capa[Ns] >= 0, (container = Array)
        Y[L], (container=Array, binary=true)
        Alpha[L, K], (container=Array, binary=true)
        Beta[Ns], (container = Array, binary=true)
        Curr_limit[T, L, K], (container=Array, binary=true)
        Slack[T, L, K], (container = Array)
        DSO_fixed_costs
        DSO_loss_costs
        DSO_op_limits
    end)
    return
end
function _add_upper_constraints!(model::JuMP.AbstractModel)

end

# Add the possibility to have storage 
# Activate or deactivate PQ diagram 
# 
function _add_lower_variables!(model::JuMP.AbstractModel, )

end
function _add_lower_constraints!(model::JuMP.AbstractModel, )

end



function build_bilevel_model(network::Network, 
                            DSO_costs::DSOCosts,
                            User_costs::UserCosts;;
                            formulation = Formulation()
                            )::JuMP.AbstractModel

    @info "Building model..."
    time_model = @elapsed begin
        model = JuMP.Model(Gurobi.Optimizer)
        set_optimizer_attribute(model, "TimeLimit", 200)
        set_optimizer_attribute(model, "Presolve", 0)

        model[:network_data] = network_dict
        model[:obj_data]     = obj_dict

        _add_BusVariables!(model, formulation.production)
        #_add_BranchVariables!(model, formulation.networkgraph, formulation.condvars)
        _add_RefVoltages!(model)
        _add_LoadOverSatisfaction!(model)
        #_add_PowerFlowEqs!(model, formulation.powerflow, formulation.networkgraph, formulation.condvars)
            
    end

    @info @sprintf("Built model in %.2f seconds", time_model)


    return model
end



function bilevel!(network::Network, 
                  DSO_costs::DSOCosts,
                  User_costs::UserCosts; 
                  BILEVEL::Bool=false, RADIALITY::Bool=true)
    

    

end