using ExcelFiles, DataFrames, SparseArrays, GraphRecipes, Plots, LinearAlgebra

NETWORK=2
FEEDER=3
vertices=DataFrame(load("Manchester_data/LV_network_models/network_$NETWORK/Feeder_$FEEDER/XY_Position.xls", "Sheet1"))
edges=DataFrame(load("Manchester_data/LV_network_models/network_$NETWORK/Feeder_$FEEDER/Feeder_Data.xls", "Sheet1"))

n_edges=size(edges,1)

adjmat=sparse(edges.NodeA, edges.NodeB, ones(n_edges))
adjmat=[adjmat;zeros(length(vertices[!,1]),size(adjmat,2)-size(adjmat,1))']

graphplot(
    adjmat, 
    x=vertices.X, 
    y=vertices.Y, 
    curves=false,
    arrow=false,
    edgecolor=reshape(range(colorant"red", stop=colorant"blue",length=n_edges), 1, n_edges )
    )
