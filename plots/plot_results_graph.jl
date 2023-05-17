using GraphRecipes
using Graphs
using Plots

falf = [sum(value(P_conductor_f[k, l]) for k in 1:K) for l in 1:L]
falb = [sum(value(P_conductor_b[k, l]) for k in 1:K) for l in 1:L]
active_losses = falf + falb
frlf = [sum(value(Q_conductor_f[k, l]) for k in 1:K) for l in 1:L]

edge_label_dict = Dict(line_ends[l] => round.(falf[l] * BASE_POWER; digits=3) for l in 1:L)
edge_width = Dict(line_ends[l] => (value(x[l]) ≈ 1 ? 1 : 0.1) for l = 1:L)

x_pos = [0 for n = 1:N]
y_pos = [0 for n = 1:N]

g = SimpleDiGraph(N)
for l in 1:L
    add_edge!(g, line_ends[l][1], line_ends[l][2])
end

next_x = [1 for n = 1:N]
function minimum_ws(g, node, depth)
    x_pos[node] = next_x[depth]
    y_pos[node] = depth
    next_x[depth] += 2
    for c in neighbors(g, node)
        minimum_ws(g, c, depth + 1)
    end
end

minimum_ws(g, 1, 1)


plot!(size=(2000, 2000))

graphplot(
    adjacency_matrix(g),
    x=x_pos,
    y=y_pos,
    curves=false,
    arrow=true,
    edgelabel=edge_label_dict,
    names=[(round.(-P_demand[n] * BASE_POWER; digits=3)) for n = 1:N],
    nodeshape=:rect,
    nodecolor=:red,
    ew=edge_width
)

