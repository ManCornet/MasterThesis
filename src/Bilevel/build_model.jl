# using Printf, Logging
# import JuMP, Gurobi
# include("../structs.jl")
# include("./formulation/structs.jl")
# include("./formulation/variables.jl")
# include("./formulation/constraints.jl")
# include("./formulation/objective.jl")



# Remove for large models the ability to have names inside the model
# look at the UnitCommitment.jl package



function build_model(   simulation::Simulation;
                        TimeLimit::Float64 = 600.0,
                        MIPGap::Float64 = 1e-2,
                        MIPFocus::Int64 = 1, 
                        set_names::Bool = false,
                        bilevel_mode=BilevelJuMP.StrongDualityMode() # other option: BilevelJuMP.SOS1Mode()
                    )::Union{Nothing, JuMP.AbstractModel}

    # ====================== Set up the Gurobi solver =====================

    if simulation.bilevel
        model = BilevelJuMP.BilevelModel(Gurobi.Optimizer, mode=bilevel_mode)
        upper = Upper(model)
        lower = Lower(model)
    else
        model = JuMP.Model(Gurobi.Optimizer)
        JuMP.set_string_names_on_creation(model, set_names)
        upper = lower = model
    end

    #JuMP.set_optimizer_attribute(model, "mode", Mode)
    JuMP.set_optimizer_attribute(model, "TimeLimit", TimeLimit)
    JuMP.set_optimizer_attribute(model, "MIPGap", MIPGap)
    JuMP.set_optimizer_attribute(model, "MIPFocus", MIPFocus)

    formulation = simulation.formulation
    if isa(formulation.convexity, NonConvex)
        JuMP.set_optimizer_attribute(model, "NonConvex", 2)  
    end

    model[:network_data]     = simulation.network
    model[:network_topology] = simulation.network_topology
    model[:DSO_costs]        = simulation.DSO_costs 
    model[:User_costs]       = simulation.User_costs 
    model[:time_steps]       = simulation.nb_time_steps #ATTENTION PUT BACK TO NB TIME STEPS
    model[:bilevel] = simulation.bilevel
    model[:storage] = simulation.storage

    model[:delta_t]          = simulation.delta_t
    model[:nb_sign_days]     = simulation.nb_sign_days

    
    # =========================== Build the model =========================

    @info "Building model..." 
    time_model = @elapsed begin

        # -- Add the variables of the model --
        _add_BusVariables!(model)
        _add_BranchVariables!(upper, formulation.powerflow)
        _add_CondChoiceVariables!(upper, formulation.topology_choice, formulation.graph_type)
        _add_VoltageOpConstraints!(upper, formulation.v_constraints)
        _add_CurrentOpConstraints!(upper, formulation.topology_choice, formulation.i_constraints)
        _add_RadialityVariables!(upper, formulation.topology_choice, formulation.radiality)
        
        # -- Add the constraints of the upper model --
           #_add_LoadOverSatisfaction!(Upper(model))
        _add_RefVoltages!(upper)
        _add_SubstationConstraints!(upper, formulation.convexity)
        _add_PowerBalanceConstraints!(upper, formulation.topology_choice, formulation.powerflow)
        _add_RotatedConicConstraints!(upper, formulation.powerflow, formulation.convexity)
        _add_PowerFlowConstraints!(upper, formulation.topology_choice, formulation.graph_type, formulation.powerflow)
        _add_RadialityConstraints!( upper, 
                                    formulation.graph_type,
                                    formulation.topology_choice,
                                    formulation.radiality)::Nothing

        # -- Add the constraints of the lower model --
        _add_LowerConstraints!(lower)

        _set_Objective!(model, formulation.i_constraints, formulation.v_constraints) 
    
        #_add_RadialityVariables!(model, formulation.topology_choice, formulation.radiality)
        #_add_RadialityConstraints!( model, 
                                    # formulation.graph_type,
                                    # formulation.topology_choice,
                                    # formulation.radiality)::Nothing
        
    end
    @info @sprintf("Built model in %.2f seconds", time_model);
    return model
end

function solve_model(model::JuMP.AbstractModel,
                    power_flow::PowerFlowFormulation;
                    verbose=true)

    try 
        JuMP.optimize!(model)
        JuMP.solution_summary(model, verbose=verbose)  
    catch 
        println("problem unbounded or infeasible")
        return
    finally
        _update_buses!(model)
        _update_lines!(model, power_flow)
        return JuMP.objective_value(model), JuMP.solve_time(model)   
    end
end

function _update_buses!(model::JuMP.AbstractModel)

    network   = model[:network_data]
    buses     = network.buses 
    Nu, Ns, T = network.nb_loads, network.nb_substations, model[:time_steps] 
    N = Ns + Nu 

    # Update substation buses
    # don't forget to add the fact that a substation can be built or not
    for (i,b) in enumerate(buses[1:Ns])
        beta = isapprox(JuMP.value.(model[:Beta][i]), 1; rtol = 1e-4)
        b.built  = (b.built || beta) 
        b.V_magn = sqrt.(JuMP.value.(model[:V_sqr][:, i]))

        if b.built 
            b.S_rating += (beta ? JuMP.value(model[:S_sub_capa][i]) : 0) # Additional capa
            b.P_sup = JuMP.value.(model[:P_sub][:, i])
            b.Q_sup = JuMP.value.(model[:Q_sub][:, i])
        else 
            b.S_rating = 0.0 
            b.P_sup = zeros(Float64, T)
            b.Q_sup = zeros(Float64, T)
        end
       
    end

    # Update load buses 
    for (i, b) in enumerate(buses[(Ns+1):N])
        b.V_magn = sqrt.(JuMP.value.(model[:V_sqr][:, Ns + i]))
        if !isnothing(b.PV_installation)
            PV = b.PV_installation
            PV.capa = JuMP.value.(model[:p_pv_max][i])
            PV.P = JuMP.value.(model[:p_pv][:, i])
            PV.Q = zeros(Float64, T) # for now leave like this
        end
        if !isnothing(b.storage)
            storage = b.storage
            storage.capa = JuMP.value.(model[:storage_capacity][i])
            storage.P = JuMP.value.(model[:p_storage][:, i])
            storage.state = JuMP.value.(model[:storage_state][:, i])
        end

        # storage update will be at this location also
    end

    return
end


function _update_lines!(model::JuMP.AbstractModel, power_flow::PowerFlowFormulation)

    lines      = model[:network_data].lines
    conductors = model[:network_data].conductors

    for l in lines

        i = l.edge.id
        # update the built field 
        l.built = isapprox(JuMP.value(sum(model[:Alpha][i, :])), 1; rtol = 1e-4)

        if l.built 
            # update the conductor field
            cond_idx = findfirst(x->isapprox(x, 1; rtol=1e-4), JuMP.value.(model[:Alpha][i, :]))
            l.conductor = conductors[cond_idx]

            # update the cost field
            l.cost = l.conductor.cost * l.length

            # update the current field
            l.I_magn = sqrt.(vec(JuMP.value.(sum(model[:I_sqr_k][:, i, :], dims=2))))

            # Update the power field + direction of power flow : 1 if from_node 
            # Change edge struct also (normalement ça va se changer aussi dans network_topology)
    
            if isa(power_flow, BIM)
                l.P_send = vec(JuMP.value.(sum(model[:P_ij_k][:, i, :], dims=2)))
                l.P_rec = vec(JuMP.value.(sum(model[:P_ji_k][:, i, :], dims=2)))
                l.Q_send = vec(JuMP.value.(sum(model[:Q_ij_k][:, i, :], dims=2)))
                l.Q_rec = vec(JuMP.value.(sum(model[:Q_ji_k][:, i, :], dims=2)))
        
            elseif isa(power_flow, BFM)
                l.P_send = vec(JuMP.value.(sum(model[:P_ij_k][:, i, :], dims=2)))
                l.P_rec = l.P_send .- l.conductor.r * l.length .* l.I_magn.^2
                l.Q_send = vec(JuMP.value.(sum(model[:Q_ij_k][:, i, :], dims=2)))
                l.Q_rec = l.Q_send .- l.conductor.x * l.length .* l.I_magn.^2
            end
        end

    end
    return
end