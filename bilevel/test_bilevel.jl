using JuMP, BilevelJuMP, Gurobi
model = BilevelModel(Gurobi.Optimizer, mode=BilevelJuMP.SOS1Mode())
@variable(Upper(model), x)
@variable(Lower(model), y)
@objective(Upper(model), Min, x - 4y)
@constraints(Upper(model), begin
    x >= 0
end)
@objective(Lower(model), Min, y)
@constraints(Lower(model), begin
    y >= 0
    -x - y <= -3
    -2x + y <= 0
    2x + y <= 12
    -3x + 2y <= -4
end)
optimize!(model)
objective_value(model) # = 3 * (3.5 * 8/15) + 8/15 
value(x) # = 3.5 * 8/15 
value(y) # = 8/15
