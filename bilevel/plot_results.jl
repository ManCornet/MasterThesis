 # graphical display

# PV understanding
for u in Nu
    plt = Plots.plot()
    plot!(T, P_CONSUMPTION[u, :], label="P_conso", title="Bus $(u)")
    #plot!(plt, T, Q_CONSUMPTION[u, :], label="Q_conso")
    plot!(plt, T, [value.(model[:p_pv_max])[u] for t in T], label="pv max")
    plot!(plt, T, value.(model[:p_pv])[u, :], label="p_pv")
    plot!(plt, T, PV_PRODUCTION[u, :] * value.(model[:p_pv_max])[u], label="max PV prod")
    #plot!(plt, T, value.(model[:q_pv])[u, :], label="q_pv")
    plot!(plt, T, value.(model[:p_exp])[u, :], label="p_exp")
    display(plt)
end

# -- SHAPES OF THE NODES OF THE NETWORK --
house_nodeshape(x_i, y_i, s) = [(x_i + 0.7s * dx, y_i + 0.7s * dy) for (dx, dy) in [(1, 1), (0, 1.6), (-1, 1), (-1, -1), (1, -1), (1, 1)]]
subs_nodeshape(x_i, y_i, s) = [(x_i + 1s * dx, y_i + 1s * dy) for (dx, dy) in [(1, 1), (-1, 1), (-1, -1), (1, -1), (1, 1)]]

# -- FUNCTION THAT PRINTS THE NETWORK --
function print_network(nb_nodes::Int64, 
                        nb_sub_nodes::Int64, 
                        line_ends::Vector{Tuple{Int64,Int64}}, 
                        Alpha, P_ij, I_sqr_k, R, p_pv, P_sub, P_CONSUMPTION)

    
    L = 1:length(line_ends)
    N = 1:nb_nodes
    Ns = 1:nb_sub_nodes
    Nu = setdiff(N, Ns)
    
    chosen_lines = sum(Alpha, dims=2) 
    edge_label_dict = Dict{Tuple{Int64,Int64},Float64}()
  
    # Find values of line powers

    line_power = [(P_ij[l] >= 0) ? P_ij[l] : abs(P_ij[l] - sum(R[l, k] * I_sqr_k[l] for k in K)) for l in L]
    
    # Creating the graph
    edge_width = Dict()
    g = SimpleDiGraph(nb_nodes)
    for l in L
        if P_ij[l] < 0
            add_edge!(g, line_ends[l][2], line_ends[l][1])
            key = (line_ends[l][2], line_ends[l][1])
        else
            add_edge!(g, line_ends[l][1], line_ends[l][2])
            key = line_ends[l]
        end
        if isapprox(chosen_lines[l], 1; rtol=1e-4) # TO BE MODIFIED
            edge_label_dict[key] = abs(round(line_power[l]; digits=4) .* BASE_POWER)
        end
        edge_width[key] = (isapprox(chosen_lines[l], 1; rtol=1e-4) ? 1 : 0.1)
    end

    colors = [colorant"lightseagreen", colorant"orange"]
    house_labels = ["+$(round(value.(p_pv)[n]*BASE_POWER, digits=3)) \n\n -$(round(P_CONSUMPTION[n]*BASE_POWER, digits=3))" for n in Nu]
    subs_labels = ["+$(round(value.(P_sub)[n]*BASE_POWER, digits=3))" for n in Ns]
    node_labels = [subs_labels; house_labels]
    node_shapes = [[subs_nodeshape for _ in Ns]; [house_nodeshape for _ in Nu]]

    # See components.jl from Plots.jl to have all the arguments
    graph = graphplot(adjacency_matrix(g),
        x=[[df_bus.x[n] + 0.5 for n in Ns]
            [df_bus.x[n] for n in Nu]],             # x-coordinate of the nodes
        y=df_bus.y,                                 # y-coordinate of the nodes
        nodesize=0.1,
        nodestrokewidth=0,                                        # coutour line width of the node
        edgestyle=:solid,
        nodealpha=1,                                        # transparency of node color
        names=node_labels,                              # node label
        nodeshape=node_shapes,                              # :circle, :ellipse, :hexagon
        nodecolor=colors[[[1 for _ in Ns]; [2 for _ in Nu]]],
        curves=false,                                    # if an edge is curved or not
        #arrow=Plots.arrow(:closed, 0.8, 0.8),                 # other choices : :open, :closed
        ew=edge_width,
        shorten=0.05,
        edgelabel=edge_label_dict,
        edgelabeloffset=3,
        axis_buffer=0.1,
        fontsize=10,
        size=(1000, 2000)
    )
    return graph
end

function if_small(x)
    for (i, a) in enumerate(x)
        if abs(a) < 1e-6
            x[i] = 0
        end
    end
    return x #round.(x, digits=6)
end

#using PyPlot
for t in T
    pyplot()
    network_plot = print_network(  N_size, 
                    Ns_size, 
                    line_ends, 
                    value.(model[:Alpha]), 
                    value.(model[:P_ij])[:, t],
                    value.(model[:I_sqr_k])[:, :, t],
                    R, 
                    value.(model[:p_pv])[:, t], 
                    value.(model[:P_sub])[:, t],
                    P_CONSUMPTION[:, t]
    )
    display(network_plot)
    #PyPlot.savefig("network_state_time_step_$(t).png")
end

# # log record
# dict_input = Dict(
#     "S_CONSUMPTION" => Containers.DenseAxisArray(S_CONSUMPTION, N, T),
#     "P_CONSUMPTION" => Containers.DenseAxisArray(P_CONSUMPTION, N, T),
#     "Q_CONSUMPTION" => Containers.DenseAxisArray(Q_CONSUMPTION, N, T),
#     "PV_PRODUCTION" => Containers.DenseAxisArray(PV_PRODUCTION, N, T),
#     "MAX_CURRENT" => Containers.DenseAxisArray(MAX_CURRENT, L, K),
#     "MAX_VOLTAGE" => Containers.DenseAxisArray([MAX_VOLTAGE], 1:1),
# )
# dict_output = Dict(String(k) => value.(v) for (k, v) in object_dictionary(model))
# data = merge(dict_output, dict_input)

# xlsx_output(data)

# Plot currents
for l in L
    plt = Plots.plot()
    chosen_cond_index = findfirst(i->isapprox(i, 1; rtol=1e-4), value.(model[:Alpha])[l, :])
    if isnothing(chosen_cond_index)
        continue
    end
    plot!(T,  sqrt.(if_small(value.(model[:I_sqr])[l, :])) * BASE_CURRENT, label="I [kA]", title="Line $(l)")
    plot!(plt, T, MAX_CURRENT[l, chosen_cond_index]*BASE_CURRENT * ones(T_size), label="MAX_CURRENT [kA]")
    display(plt)
end

# Display the current for each conductor
for l in L
    plt = Plots.plot()
    chosen_cond_index = findfirst(i->isapprox(i, 1; rtol=1e-4), value.(model[:Alpha])[l, :])
    title = isnothing(chosen_cond_index) ? "Line $l is not built" : "Line $l is built with conductor $chosen_cond_index"
    for k in K
        plot!(plt, T,  sqrt.(if_small(value.(model[:I_sqr_k])[l, k, :])) * BASE_CURRENT, label="I_ij_$k [kA]", title=title)
    end
    plot!(plt, T,  sum(sqrt.(if_small(value.(model[:I_sqr_k])[l, :, :])) * BASE_CURRENT, dims=1)', label="I_ij [kA]")
    display(plt)
end

# Display the power for each conductor
for l in L
    plt = Plots.plot()
    chosen_cond_index = findfirst(i->isapprox(i, 1; rtol=1e-4), value.(model[:Alpha])[l, :])
    title = isnothing(chosen_cond_index) ? "Line $l is not built" : "Line $l is built with conductor $chosen_cond_index"
    for k in K
        plot!(plt, T,  value.(model[:P_ij_k])[l, k, :]  /  MAX_VOLTAGE * BASE_CURRENT, label="P_ij_$k [kA]", title=title)
    end
    
    if isnothing(chosen_cond_index)
        continue 
    end
    plot!(plt, T, MAX_CURRENT[l, chosen_cond_index]*BASE_CURRENT * ones(T_size), label="+MAX_CURRENT [kA]")
    plot!(plt, T, -MAX_CURRENT[l, chosen_cond_index]*BASE_CURRENT * ones(T_size), label="-MAX_CURRENT [kA]")

    display(plt)
end

# Check voltage limits
for i in N
    plt = Plots.plot()
    plot!(T,  sqrt.(value.(model[:V_sqr])[i, :]), label="V [p.u.]", title="Bus $(i)")
    plot!(plt, T, MIN_VOLTAGE * ones(T_size), label="MIN_VOLTAGE [p.u.]")
    plot!(plt, T, MAX_VOLTAGE * ones(T_size), label="MAX_VOLTAGE [p.u.]")
    display(plt)
end

for i in Ns
    plt = Plots.plot()
    plot!(T,  value.(model[:P_sub])[i, :], label="P_sub", title="Bus $(i)")
    display(plt)
end

# Compute total losses per time step: 
total_losses = [sum(R[:, :] .* value.(model[:I_sqr_k])[:, :, t] .* BASE_POWER) for t in T]
# Check P_sub, Q_sub
p = Plots.plot()
plot!(T,  [sum(value.(model[:P_sub])[:, t]) for t in T] * BASE_POWER, label="P_sub [MW]")
plot!(p, T, [sum(P_CONSUMPTION[:, t]) * BASE_POWER for t in T], label="Total P_CONSUMPTION [MW]")
plot!(p, T, [sum(value.(model[:p_pv])[:, t]) * BASE_POWER for t in T], label="Total p_pv [MW]")
plot!(p, T, total_losses, label="Total Joule Losses [MW]")
plot!(p, T, [sum(value.(model[:p_pv])[:, t]) * BASE_POWER for t in T] .+  [sum(value.(model[:P_sub])[:, t])  for t in T] * BASE_POWER .- total_losses , label="p_pv + P_sub - losses[MW]")
display(p)