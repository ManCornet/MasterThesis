using XLSX, DataFrames, Plots, JuMP, Gurobi
env=Gurobi.Env()


include("functions.jl")


# PARAMETERS
# arbitrary values, to be modified
cost_PV=1800 # €/kWp
cost_storage=1000 # €/kWhp
grid_voltage = 0.230 # kV
grid=DataFrame(capacity=[16,32,64], cost=[5,12.5,30]) # [A],[EUR/year]
cost_grid_imp=0.5 # €/kWh
cost_grid_exp=0.2 # €/kWh
cost_community=0.3 # €/kWh

eff_PV=0.18
eff_Bc=0.95
eff_Bd=0.95


COMMUNITY_SIZE=1                 # homes in the community
PROFILES=range(1, length=COMMUNITY_SIZE) # must be <= 100

# careful when feeders implemented
PV_RATE=0.6         # PV penetration rate <=1

# half summer (PV), half winter (EV, HP, uCHP)
LOAD_PROFILE, PV_PROFILE=create_profiles(
    summer=true,
    PV=true,
    winter=true,
    EV=true,
    HP=false);
steps_1p=288                # nbr of steps in one profile
steps=size(LOAD_PROFILE,1)  # nbr of steps in total

COMMUNITY = false    # Is energy exchanged allowed among the community ?
TIME_STEP=1/12      # [hour/step]
MONTHS=12*20
days=Int(floor(steps/(24/TIME_STEP)))    # nbr of days in one profile
n_profiles=MONTHS*30/days    # nbr of profiles occuring within the MONTHS period


model=simulate(COMMUNITY);

# list_of_constraint_types(model)
# all_constraints(model, AffExpr, MOI.EqualTo{Float64})
# all_constraints(model, AffExpr, MOI.GreaterThan{Float64})
# all_constraints(model, AffExpr, MOI.LessThan{Float64})
# all_constraints(model, VariableRef, MOI.GreaterThan{Float64})
# all_constraints(model, VariableRef, MOI.LessThan{Float64})
# all_constraints(model, VariableRef, MOI.Integer)

print_results(model)
display_results(model,[2 3 4]) #  1:load+PV  2:powers  3:capacities  4:battery


# cost=[]
# for i in 1:2:40
#     model=simulate(i,true)
#     append!(cost,objective_value(model))
# end
# plot(cost)


# load, pv = create_profiles(EV=true);
# load
# pv
