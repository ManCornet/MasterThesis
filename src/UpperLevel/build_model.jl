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
        JuMP.set_optimizer_attribute(model, "TimeLimit", 200)
        JuMP.set_optimizer_attribute(model, "Presolve", 0)

        model[:network_data] = simulation.network
        model[:DSO_costs]    = simulation.DSO_costs 
        model[:User_costs]   = simulation.User_costs 
        model[:time_steps]   = simulation.nb_time_steps
        model[:delta_t]      = simulation.delta_t
        
        # =========================== Build the model =========================
        # -- Add the variables of the model --
        _add_BusVariables!(model, formulation.production)
        #_add_BranchVariables!(model, )

        #_add_BranchVariables!(model, formulation.networkgraph, formulation.condvars)
        #_add_RefVoltages!(model)
        #_add_LoadOverSatisfaction!(model)
        #_add_PowerFlowEqs!(model, formulation.powerflow, formulation.networkgraph, formulation.condvars)
    end

    @info @sprintf("Built model in %.2f seconds", time_model)
    

    return model
end