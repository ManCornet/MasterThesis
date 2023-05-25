module KnapsackModel

import JuMP
import Gurobi

struct _KnapsackObject 
    profit::Float64 
    weight::Float64

    function _KnapsackObject(profit::Float64, weight::Float64) 
        if weight < 0
            throw(DomainError("Weight of object cannot be negative")) 
        end
    return new(profit, weight) end
end


struct _KnapsackData 
    objects::Dict{String,_KnapsackObject}
    capacity::Float64
end


function _read_data(filename)
    d = JSON.parsefile(filename) 
    return _KnapsackData(
                        Dict(
                        k => _KnapsackObject(v["profit"], v["weight"]) for (k, v) in d["objects"]
                        ),
                        d["capacity"],
    )
end


abstract type _AbstractConfiguration end


"""
    BinaryKnapsackConfig()
    Create a binary knapsack problem where each object can be taken 0 or 1 times. 
"""
struct BinaryKnapsackConfig <: _AbstractConfiguration end

"""
    IntegerKnapsackConfig()
    Create an integer knapsack problem where each object can be taken any number of times.
"""
struct IntegerKnapsackConfig <: _AbstractConfiguration end


function _add_knapsack_variables( model::JuMP.Model, data::_KnapsackData, ::BinaryKnapsackConfig,
)
    return JuMP.@variable(model, x[keys(data.objects)], Bin)
end

function _add_knapsack_variables( model::JuMP.Model, data::_KnapsackData, ::IntegerKnapsackConfig,
)
    return JuMP.@variable(model, x[keys(data.objects)] >= 0, Int)
end

function _add_knapsack_constraints( model::JuMP.Model, data::_KnapsackData, ::_AbstractConfiguration,
)
    x = model[:x]
    JuMP.@constraint(model,capacity_constraint,
                    sum(v.weight * x[k] for (k, v) in data.objects) <= data.capacity, 
                    )
    return 
end

function _add_knapsack_objective( model::JuMP.Model, data::_KnapsackData, ::_AbstractConfiguration,
    )
    x = model[:x]
    JuMP.@objective(model, Max, sum(v.profit * x[k] for (k, v) in data.objects))
return end


function _solve_knapsack( optimizer,
    data::_KnapsackData,
    config::_AbstractConfiguration,
)
    model = JuMP.Model(optimizer) 
    _add_knapsack_variables(model, data, config) 
    _add_knapsack_constraints(model, data, config) 
    _add_knapsack_objective(model, data, config) 
    JuMP.optimize!(model)

    if JuMP.termination_status(model) != JuMP.OPTIMAL
        @warn("Model not solved to optimality")
        return nothing 
    end
    return JuMP.value.(model[:x]) 
end


function solve_knapsack(
                    optimizer,
                    data_filename::String,
                    config::_AbstractConfiguration,
                    )
    return _solve_knapsack(optimizer, _read_data(data_filename), config)
end


end
