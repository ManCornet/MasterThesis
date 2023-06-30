using Test, Graphs

P_sub = value.(model[:P_sub])
Q_sub = value.(model[:Q_sub])
S_sub = value.(model[:S_sub])
P_ij  = value.(model[:P_ij])
Q_ij  = value.(model[:Q_ij])
V_sqr = value.(model[:V_sqr])
I_sqr = value.(model[:I_sqr])


# -- BFM: OPTIMALITY VERIFICATION --
@testset verbose = true "Checking optimality" begin
    
    @testset "Conic constraint for substation $i at time step $t" for i in Ns, t in T
        @test S_sub[i, t]^2 ≈ P_sub[i, t]^2 + Q_sub[i, t]^2 rtol=1e-3
    end
    @testset "Rotated Conic constraints at line $l and time step $t" for l in L, t in T
        @test V_sqr[line_ends[l][1], t] * I_sqr[l, t] ≈ P_ij[l, t]^2 + Q_ij[l, t]^2 rtol=1e-3
    end
end


# -- BFM: RADIALITY VERIFICATION --

@testset verbose = true "Checking optimality" begin
    nb_nodes = Nu_size + length(Ns_init) + convert(Int64, sum(value.(model[:Beta])[Ns_notinit]))
    g = SimpleGraph(nb_nodes)
    for l in L
        if isapprox(value.(model[:Y])[l], 1; rtol=1e-4) # TO BE MODIFIED
            if line_ends[l][1] > length(Ns_init)
                send_end = line_ends[l][1] - length(Ns_notinit) + convert(Int64, sum(value.(model[:Beta])[Ns_notinit]))
            else 
                send_end = line_ends[l][1]
            end
            if line_ends[l][2] > length(Ns_init)
                rec_end = line_ends[l][2] - length(Ns_notinit) + convert(Int64, sum(value.(model[:Beta])[Ns_notinit]))
            else 
                send_end = line_ends[l][2]
            end
            
            add_edge!(g, send_end, rec_end)
        end
    end
    
    display(graphplot( adjacency_matrix(g), curves=false))
    @test is_tree(g)
end