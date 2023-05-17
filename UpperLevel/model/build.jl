using JuMP
using Gurobi

# Est-ce que ça sert à quelque chose de faire ça ? OUI SI PLUSIEURS FORMULATIONS A TESTER !!!

"""
    function build_model(;
        optimizer = nothing,
        formulation = Formulation(),
    )::JuMP.Model

Build the JuMP model corresponding to the given unit commitment instance.

Arguments
---------
- `optimizer`:
    the optimizer factory that should be attached to this model (e.g. Cbc.Optimizer).
    If not provided, no optimizer will be attached.
- `formulation`:
    the MIP formulation to use. By default, uses a formulation that combines
    modeling components from different publications that provides good
    performance across a wide variety of instances. An alternative formulation
    may also be provided.
- `variable_names`: 
    if true, set variable and constraint names. Important if the model is going
    to be exported to an MPS file. For large models, this can take significant
    time, so it's disabled by default.

"""

# Use expressions or something like that:
# https://jump.dev/JuMP.jl/stable/manual/expressions/


function build_model(;
                    formulation = Formulation(),
                    variable_names::Bool = false,
)::JuMP.Model

    @info "Building model..."
    time_model = @elapsed begin
        model = Model(Gurobi.Optimizer)
        set_optimizer_attribute(model, "TimeLimit", 200)
        set_optimizer_attribute(model, "Presolve", 0)

        model[:obj] = QuadExpr()

        _add_RefVoltage!(model,  MIN_VOLTAGE, MAX_VOLTAGE)
        _add_LoadOverSatisfaction!(model, formulation.powerflow, )
    end

    @info @sprintf("Built model in %.2f seconds", time_model)
    if variable_names
        _set_names!(model)
    end
     
    return model

end
