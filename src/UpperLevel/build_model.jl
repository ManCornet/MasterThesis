# using Printf, Logging
# import JuMP, Gurobi
# include("../structs.jl")
# include("./formulation/structs.jl")
# include("./formulation/variables.jl")
# include("./formulation/constraints.jl")
# include("./formulation/objective.jl")

function build_model(   simulation::Simulation;
                        formulation = Formulation(),
                        TimeLimit = 600,
                        MIPGap = 1e-2,
                        MIPFocus = 1
                    )::Union{Nothing, JuMP.Model}

    @info "Building model..."
    time_model = @elapsed begin
        # ====================== Set up the Gurobi solver =====================
        model = JuMP.Model(Gurobi.Optimizer)
        #JuMP.set_optimizer_attribute(model, "mode", Mode)
        JuMP.set_optimizer_attribute(model, "TimeLimit", TimeLimit)
        JuMP.set_optimizer_attribute(model, "MIPGap", MIPGap)
        JuMP.set_optimizer_attribute(model, "MIPFocus", MIPFocus)

        if typeof(formulation.convexity) <: NonConvex
            JuMP.set_optimizer_attribute(model, "NonConvex", 2)  
        end

        model[:network_data]     = simulation.network
        model[:network_topology] = simulation.network_topology
        model[:DSO_costs]        = simulation.DSO_costs 
        model[:User_costs]       = simulation.User_costs 
        model[:time_steps]       = 1 #simulation.nb_time_steps #ATTENTION PUT BACK TO NB TIME STEPS
        model[:delta_t]          = simulation.delta_t
        model[:nb_sign_days]     = simulation.nb_sign_days
        
        # =========================== Build the model =========================

        # -- Add the variables of the model --
        _add_BusVariables!(model, formulation.production)
        _add_BranchVariables!(model, formulation.powerflow)
        _add_CondChoiceVariables!(model, formulation.topology_choice, formulation.graph_type)

        # -- Add the constraints of the model --
        _add_RefVoltages!(model)
        if isa(formulation.production, DG)
            _add_PVOperationConstraints!(model)
        end
        #_add_LoadOverSatisfaction!(model, formulation.production)
        _add_SubstationConstraints!(model, formulation.convexity)
        _add_CurrentOpConstraints!(model, formulation.topology_choice, formulation.i_constraints)
        _add_VoltageOpConstraints!(model, formulation.v_constraints)
        _add_PowerBalanceConstraints!(model, formulation.production, formulation.powerflow)
        _add_RotatedConicConstraints!(model, formulation.powerflow, formulation.convexity)
        _add_PowerFlowConstraints!(model, formulation.topology_choice,  formulation.graph_type, formulation.powerflow)
        _set_Objective!(model, formulation.i_constraints, formulation.v_constraints)
        
    end
    @info @sprintf("Built model in %.2f seconds", time_model);
    return model
end

function solve_model(model::JuMP.AbstractModel; verbose=false)

    JuMP.optimize!(model)
    JuMP.solution_summary(model, verbose=verbose)

    try
        
    catch


    end
    if JuMP.termination_status(model) == JuMP.MOI.OPTIMAL
        var_values = Dict(  k => JuMP.value.(v) for (k, v) 
                            in JuMP.object_dictionary(model) 
                            if (v isa AbstractArray{JuMP.VariableRef} || 
                                v isa JuMP.VariableRef)
                        )
        
        # var_sets = Dict("V_sqr"     => ["node"], 
        #                 "I_sqr_k"   => ["line", "conductor"],
        #                 "I_sqr"     => ["line"],
        #                 "P_G"       => ["node"],
        #                 "Q_G"       => ["node"],
        #                 "S_G"       => ["node"],
        #                 "P_ij_k"    => ["line", "conductor"],
        #                 "Q_ij_k"    => ["line", "conductor"],
        #                 "P_ij"      => ["line"],
        #                 "Q_ij"      => ["line"],
        #                 "y_send"    => ["line"],
        #                 "y_rec"     => ["line"],
        #                 "alpha"     => ["line", "conductor"],
        #                 "beta"      => ["node"]
        #                 )
    
        return JuMP.objective_value(model), var_values, JuMP.solve_time(model) 
    
    elseif JuMP.termination_status(model) == DUAL_INFEASIBLE
        println("problem unbounded")
        return 
    elseif JuMP.termination_status(model) == JuMP.MOI.INFEASIBLE
        println("problem infeasible")
        return
    end
end 
