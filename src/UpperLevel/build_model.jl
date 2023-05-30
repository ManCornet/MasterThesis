using Printf, Logging
import JuMP, Gurobi
include("../structs.jl")
include("./formulation/structs.jl")
include("./formulation/variables.jl")
include("./formulation/constraints.jl")
include("./formulation/objective.jl")

function build_model(   simulation::Simulation;
                        formulation = Formulation()
                    )::JuMP.Model

    @info "Building model..."
    time_model = @elapsed begin
        # ====================== Set up the Gurobi solver =====================
        model = JuMP.Model(Gurobi.Optimizer)
        JuMP.set_optimizer_attribute(model, "TimeLimit", 600)
        JuMP.set_optimizer_attribute(model, "MIPGap", 1e-2)
        JuMP.set_optimizer_attribute(model, "MIPFocus", 1)
        if typeof(formulation.convexity) <: NonConvex
            JuMP.set_optimizer_attribute(model, "NonConvex", 2)  
        end

        model[:network_data] = simulation.network
        model[:network_topology] = simulation.network_topology
        model[:DSO_costs]    = simulation.DSO_costs 
        model[:User_costs]   = simulation.User_costs 
        model[:time_steps]   = 3 #simulation.nb_time_steps ATTENTION PUT BACK TO NB TIME STEPS
        model[:delta_t]      = simulation.delta_t
        
        # =========================== Build the model =========================
        # -- Add the variables of the model --
        _add_BusVariables!(model, formulation.production)
        _add_BranchVariables!(model, formulation.powerflow)
        _add_CondChoiceVariables!(model, formulation.topology_choice, formulation.graph_type)

        _add_RefVoltages!(model)
        if isa(formulation.production, DG)
            _add_PVOperationConstraints!(model)
        end
        _add_LoadOverSatisfaction!(model, formulation.production)
        _add_SubstationConstraints!(model, formulation.convexity)
        _add_CurrentOpConstraints!(model, formulation.topology_choice, formulation.i_constraints)
        _add_VoltageOpConstraints!(model, formulation.v_constraints)
        _add_PowerBalanceConstraints!(model, formulation.production, formulation.powerflow)
        _add_RotatedConicConstraints!(model, formulation.powerflow, formulation.convexity)
        _add_PowerFlowConstraints!(model, formulation.topology_choice,  formulation.graph_type, formulation.powerflow)
        
    end
    @info @sprintf("Built model in %.2f seconds", time_model);

    return model
end