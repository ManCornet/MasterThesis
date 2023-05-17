using JuMP, BilevelJuMP, Gurobi
using GraphRecipes, Graphs, Plots, Dates, DataFrames
include("parameters.jl")
include("export_xlsx.jl")


# model = BilevelModel(Gurobi.Optimizer, mode=BilevelJuMP.SOS1Mode())
# set_silent(model)
model = BilevelModel(Gurobi.Optimizer, mode=BilevelJuMP.StrongDualityMode())


@variables(Upper(model), begin
    Alpha[l=L], Bin # build line ij (1) or not (0)
    # BETA[i=Ns] # build substation (1) or not (0)
    P_line[l=L, t=T] # active power from i to j
    Q_line[l=L, t=T] # reactive power from i to j
    S_MAX_build_line[l=L] >= 0 # max apparent power for each line when alpha=1
    S_sub[i=Ns] >= 0 # max apparent power transiting at substation i
    P_sub[i=Ns, t=T] # active power from transmission grid
    P_sub_imp[i=Ns, t=T] >= 0 # active power from transmission grid (positive part of P_sub)
    P_sub_exp[i=Ns, t=T] >= 0 # active power to transmission grid (negative part of P_sub)
    Q_sub[i=Ns, t=T] # reactive power from transmission grid 
    Q_sub_imp[i=Ns, t=T] >= 0 # reactive power from transmission grid (positive part of Q_sub)
    Q_sub_exp[i=Ns, t=T] >= 0 # reactive power to transmission grid (negative part of Q_sub)
    CO2_budget >= 0
    DSO_costs >= 0
end)
@variables(Lower(model), begin
    p_imp[i=Nu, t=T] >= 0 # active power imported at time t
    q_imp[i=Nu, t=T]    # reactive power imported at time t
    p_exp[i=Nu, t=T] >= 0 # active power exported at time t
    q_exp[i=Nu, t=T]    # reactive power exported at time t
    s_grid_max[i=Nu] >= 0 # kVA, grid imp/exp capacity
    p_pv[i=Nu, t=T] >= 0 # active power generated at time t
    q_pv[i=Nu, t=T]    # reactive power generated at time t
    s_conv_pv[i=Nu] >= 0 # kVA, PV converter max capacity
    p_pv_max[i=Nu] >= 0 # kWp, PV capacity
    user_costs[i=Nu]
end)

@objective(Upper(model), Min, DSO_costs)
@objective(Lower(model), Min, sum(user_costs))

@constraints(Upper(model), begin
    DSO_costs == sum(Alpha[l] * LINE_COST[l] for l in L) +
                 sum(S_sub[i] * SUBSTATION_INSTALLATION_COST for i in Ns) # max coef: 27500 ALpha[20]
    DSO_costs * (1 + DSO_INTEREST_RATE)^AMORTIZATION_DSO <= AMORTIZATION_DSO * (sum(s_grid_max[i] * GRID_CONNECTION_COST for i in Nu) +
                                                                        sum(p_imp[i, t] * IMP_ELECTRICITY_DSO_COST - p_exp[i, t] * EXP_ELECTRICITY_DSO_COST for i in Nu, t in T) * TIME_STEP * DAYS_A_YEAR) # max coef: GRID_CONNECTION_COST*WEIGHT_DSO_REVENUES s_grid_max[:]
    CO2_budget == sum(P_sub_imp[i, t] * CO2_IMP_SUBSTATION for i in Ns, t in T) +
                  sum(S_sub[i] * CO2_SUBSTATION for i in Ns) +
                  sum(Alpha[l] * CO2_LINE for l in L) +
                  sum(p_pv_max[i] * CO2_PV for i in Nu) +
                  sum(s_conv_pv[i] * CO2_CONVERTER for i in Nu) # max coef: CO2_xxx
    [i in Ns, t in T], P_sub[i, t] == sum(P_line[l, t] for l in Omega_sending[i]) # max coef: 1
    [i in Ns, t in T], Q_sub[i, t] == sum(Q_line[l, t] for l in Omega_sending[i]) # max coef: 1
    [i in Nu, t in T], p_imp[i, t] - p_exp[i, t] == -sum(P_line[l, t] for l in Omega_sending[i]) + sum(P_line[l, t] for l in Omega_receiving[i]) # max coef: 1
    [i in Nu, t in T], q_imp[i, t] - q_exp[i, t] == -sum(Q_line[l, t] for l in Omega_sending[i]) + sum(Q_line[l, t] for l in Omega_receiving[i]) # max coef: 1
    # [l in L, t in T], P_line[l, t] <= Alpha[l] * S_MAX_line[l] # TODO Remove if Alpha in cone constraint works.
    # [l in L, t in T], P_line[l, t] >= -Alpha[l] * S_MAX_line[l]
    # [l in L, t in T], Q_line[l, t] <= Alpha[l] * S_MAX_line[l]
    # [l in L, t in T], Q_line[l, t] >= -Alpha[l] * S_MAX_line[l]
    # [l in L, t in T], [Alpha[l] * S_MAX_line[l], P_line[l, t], Q_line[l, t]] in SecondOrderCone()
    [l in L], S_MAX_build_line[l] <= Alpha[l] * S_MAX_LINE[l]  # max coef: S_MAX_line
    [l in L, t in T], [S_MAX_build_line[l], P_line[l, t], Q_line[l, t]] in SecondOrderCone()
    [i in Ns], S_sub[i] <= S_MAX[i] # Substation capacity
    [i in Ns, t in T], P_sub[i, t] == P_sub_imp[i, t] - P_sub_exp[i, t] # Decomposition of active power exchange at substation in positive and negative parts
    [i in Ns, t in T], Q_sub[i, t] == Q_sub_imp[i, t] - Q_sub_exp[i, t] # Decomposition of reactive power exchange at substation in positive and negative parts
    [i in Ns, t in T], [S_sub[i], P_sub[i, t], Q_sub[i, t]] in SecondOrderCone() # Bounding the apparent power flow at the substation to the max substation capacity
    sum(Alpha) == Nu_size
end)
if MAX_CO2_BUDGET > 0
    @constraint(Upper(model), CO2_budget <= MAX_CO2_BUDGET)
end
@constraints(Lower(model), begin
    [i in Nu], user_costs[i] == sum(p_imp[i, t] * (IMP_ELECTRICITY_ENRG_COST + IMP_ELECTRICITY_DSO_COST) - p_exp[i, t] * (EXP_ELECTRICITY_ENRG_COST - EXP_ELECTRICITY_DSO_COST) for t in T) * TIME_STEP * DAYS_A_YEAR +
                                (s_conv_pv[i] * CONVERTER_COST / AMORTIZATION_CONVERTER + p_pv_max[i] * PV_COST / AMORTIZATION_PV + s_grid_max[i] * GRID_CONNECTION_COST)
    [i in Nu, t in T], p_imp[i, t] - p_exp[i, t] + p_pv[i, t] == P_CONSUMPTION[i, t]
    [i in Nu, t in T], q_imp[i, t] - q_exp[i, t] + q_pv[i, t] == Q_CONSUMPTION[i, t]
    # [i in Nu, t in T], [s_grid_max[i], p_imp-p_exp[i, t], q_imp-q_exp[i, t]] in SecondOrderCone() # see next constraints
    [i in Nu, t in T], p_imp[i, t] <= s_grid_max[i]
    [i in Nu, t in T], q_imp[i, t] <= s_grid_max[i]
    [i in Nu, t in T], -q_imp[i, t] <= s_grid_max[i]
    [i in Nu, t in T], p_exp[i, t] <= s_grid_max[i]
    [i in Nu, t in T], q_exp[i, t] <= s_grid_max[i]
    [i in Nu, t in T], -q_exp[i, t] <= s_grid_max[i]

    [i in Nu, t in T], q_pv[i, t] <= PV_MAX_Q * p_pv_max[i]
    [i in Nu, t in T], -q_pv[i, t] <= PV_MAX_Q * p_pv_max[i]
    [i in Nu, t in T], p_pv[i, t] <= p_pv_max[i] + PV_SLOPE * q_pv[i, t]
    [i in Nu, t in T], p_pv[i, t] <= p_pv_max[i] - PV_SLOPE * q_pv[i, t]
    # [i in Nu, t in T], [s_conv_pv[i], p_pv[i, t], q_pv[i, t]] in SecondOrderCone() # see next constraint
    [i in Nu], p_pv_max[i] <= s_conv_pv[i]

    [i in Nu, t in T], p_pv[i, t] <= PV_PRODUCTION[i, t] * p_pv_max[i]
end)

# m2=BilevelJuMP._build_single_model(model)

optimize!(model)
println(objective_value(model))


# graphical display
for t in T
    falf = [value(P_line[l, t]) for l in L]

    edge_label_dict = Dict{Tuple{Int64,Int64},Float64}()
    edge_width = Dict(line_ends[l] => (value(Alpha[l]) ≈ 1 ? 1 : 0.1) for l = L)

    g = SimpleDiGraph(N_size)
    for l in L
        add_edge!(g, line_ends[l][1], line_ends[l][2])
        if abs(falf[l]) > 1e-6 # TO BE MODIFIED        
            edge_label_dict[line_ends[l]] = round.(falf[l]; digits=3)
        end
    end

    box_size = 8
    colors = [:green, :red] # substations in green, users in red
    graph = graphplot(
        adjacency_matrix(g),
        x=nodes_location[1] / box_size,
        y=nodes_location[2] / box_size,
        curves=false,
        arrow=true,
        edgelabel=edge_label_dict,
        names=round.([[value.(S_sub)[n] for n = Ns]; [S_CONSUMPTION[n, t] for n = Nu]]; digits=3),
        nodeshape=:rect,
        nodecolor=colors[[[1 for _ in Ns]; [2 for _ in Nu]]],
        ew=edge_width,
        size=(1000, 1000)
    )
    display(graph)
end



# log record
dict_input = Dict(
    "S_CONSUMPTION" => Containers.DenseAxisArray(S_CONSUMPTION, N, T),
    "P_CONSUMPTION" => Containers.DenseAxisArray(P_CONSUMPTION, N, T),
    "Q_CONSUMPTION" => Containers.DenseAxisArray(Q_CONSUMPTION, N, T),
    "PV_PRODUCTION" => Containers.DenseAxisArray(PV_PRODUCTION, N, T),
    "S_MAX_LINE" => Containers.DenseAxisArray(S_MAX_LINE, L))
dict_output = Dict(String(k) => value.(v) for (k, v) in object_dictionary(model))
data = merge(dict_output, dict_input)

xlsx_output(data)
