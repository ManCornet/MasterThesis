function create_profiles(;summer=true, PV=true, winter=true, EV=false, HP=false, uCHP=false)

    load_profile=Matrix{Float64}(undef, 0, COMMUNITY_SIZE) # nbr of columns must be given
    if summer
        summer_load=XLSX.readdata("Manchester_data/LCT_profiles/Summer_Load_Profiles.xlsx","Sheet1","A1:CV288")
        load_profile=[load_profile;summer_load[:,PROFILES]]
        steps_1p=size(summer_load,1) # nbr of steps in one profile
    end
    if winter
        winter_load=XLSX.readdata("Manchester_data/LCT_profiles/Winter_Load_Profiles.xlsx","Sheet1","A1:CV288")
        load_profile=[load_profile;winter_load[:,PROFILES]]
        steps_1p=size(winter_load,1) # nbr of steps in one profile
    end
    load_profile=Array{Float64}(load_profile)
    steps=size(load_profile,1)
    if winter && EV
        winter_EV=XLSX.readdata("Manchester_data/LCT_profiles/Winter_EV_Profiles.xlsx","Sheet1","A1:CV288")
        load_profile[steps-steps_1p+1:steps,:].+=winter_EV[:,PROFILES]
    end
    if winter && HP
        winter_HP=XLSX.readdata("Manchester_data/LCT_profiles/Winter_EHP_Profiles.xlsx","Sheet1","A1:CV288")
        load_profile[steps-steps_1p+1:steps,:]+=winter_HP[:,PROFILES]
    end
    if winter && uCHP
        winter_uCHP=XLSX.readdata("Manchester_data/LCT_profiles/Winter_uCHP_Profiles.xlsx","Sheet1","A1:CV288")
        load_profile[steps-steps_1p+1:steps,:]+=winter_uCHP[:,PROFILES]
    end

    pv_profile=zeros(size(load_profile))
    if summer && PV
        summer_PV=XLSX.readdata("Manchester_data/LCT_profiles/Summer_PV_Profiles.xlsx","Sheet1","A1:CV288")
        pv_profile[1:steps_1p,:]=summer_PV[:,PROFILES]
        pv_profile[:,Int(round(PV_RATE*COMMUNITY_SIZE))+1:COMMUNITY_SIZE].=0
    end
    
    return load_profile, pv_profile
end


function simulate(community)
    
    if community
        islanded_costs=value.(simulate(false)[:costs])

    end

    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "MIPGap", 0.01)
    
    @variables(model, begin
            capacity_PV[1:COMMUNITY_SIZE]>=0
            capacity_storage[1:COMMUNITY_SIZE]
            storage_discharge[1:steps,1:COMMUNITY_SIZE]>=0
            storage_charge[1:steps,1:COMMUNITY_SIZE]<=0
            storage_SOC[1:steps,1:COMMUNITY_SIZE]>=0
            grid_imported[1:steps,1:COMMUNITY_SIZE]>=0
            grid_exported[1:steps,1:COMMUNITY_SIZE]<=0
            grid_capacity[1:size(grid,1),1:COMMUNITY_SIZE]>=0, (Int)
            community_imported[1:steps,1:COMMUNITY_SIZE]>=0
            community_exported[1:steps,1:COMMUNITY_SIZE]<=0
            costs[1:COMMUNITY_SIZE]
            diff_costs[1:COMMUNITY_SIZE]>=0
        end)
    
    for n in 1:COMMUNITY_SIZE
        @constraint(model, LOAD_PROFILE[:,n].==PV_PROFILE[:,n]+storage_discharge[:,n]+storage_charge[:,n]+grid_imported[:,n]+grid_exported[:,n]+community_imported[:,n]+community_exported[:,n])
        @constraint(model, grid_imported[:,n].<=sum(grid_capacity[:,n].*grid.capacity)*grid_voltage)
        @constraint(model, -grid_exported[:,n].<=sum(grid_capacity[:,n].*grid.capacity)*grid_voltage)
        @constraint(model, sum(grid_capacity[:,n])<=1)
        @constraint(model, PV_PROFILE[:,n].<=eff_PV*capacity_PV[n]) # réflexion inverse
        @constraint(model, storage_SOC[:,n].<=capacity_storage[n])

        for k in 0:days-1   # for each day in one profile
            @constraint(model, [i=(k*steps_1p+2):(k+1)*steps_1p],
                storage_SOC[i,n].==storage_SOC[i-1,n]-eff_Bc*storage_charge[i,n]*TIME_STEP-(1/eff_Bd)*storage_discharge[i,n]*TIME_STEP )    
            @constraint(model, # boundary condition
                storage_SOC[(k*steps_1p+1),n].==storage_SOC[(k+1)*steps_1p,n]-eff_Bc*storage_charge[(k*steps_1p+1),n]*TIME_STEP-(1/eff_Bd)*storage_discharge[(k*steps_1p+1),n]*TIME_STEP )
        end
    end
    if community
        @constraint(model, sum(community_imported,dims=2).==-sum(community_exported,dims=2))
        # @constraint(model, costs.<=islanded_costs)
        @constraint(model, diff_costs.>=costs-islanded_costs)
    else
        @constraint(model, community_imported.==0)
        @constraint(model, community_exported.==0)
    end
      
    @constraint(model, costs.== cost_PV*capacity_PV
                                +cost_storage*capacity_storage
                                +cost_grid_imp*som(grid_imported,1)*TIME_STEP*n_profiles
                                +cost_grid_exp*som(grid_exported,1)*TIME_STEP*n_profiles
                                +som(grid.cost.*grid_capacity,1)*MONTHS/12
                                +cost_community*som(community_imported+community_exported,1)*TIME_STEP*n_profiles # somme totale nulle
    )
    @objective(model, Min, sum(costs))#+sum(diff_costs))
    
    
    optimize!(model)

    # GRBgetconstrs(model,...)
    # if community
    #     write_to_file(model, "with_objective.mps")
    #     aaa=hcat(value.(model[:costs]),islanded_costs)
    #     println(aaa)
    #     iii=findall(aaa[:,1].>aaa[:,2])
    #     println(iii)
    #     # @constraint(model, costs[iii].<=islanded_costs[iii])
    #     # optimize!(model)
    #     # aaa=hcat(value.(model[:costs]),islanded_costs)
    #     # iii=findall(aaa[:,1].>aaa[:,2])
    #     # println(iii)

    # end
    # générer fichiers .mps


    # solution_summary(model)
        
    return model
end

function print_results(model)
    grid_net=(som(value.(model[:grid_imported]),1)+som(value.(model[:grid_exported]),1))*TIME_STEP*n_profiles
    grid_net_cost=(cost_grid_imp*som(value.(model[:grid_imported]),1)+cost_grid_exp*som(value.(model[:grid_exported]),1))*TIME_STEP*n_profiles
    comm_net=(som(value.(model[:community_imported]),1)+som(value.(model[:community_exported]),1))*TIME_STEP*n_profiles
    
    println("")
    printx("solve_time ",solve_time(model)) # computational time
    # println("relative_gap ",relative_gap(model)) # not with Gurobi
    # println("simplex_iterations ", simplex_iterations(model))
    println("barrier_iterations ", barrier_iterations(model))
    
    printx("")
    printx("PV        : ", value.(model[:capacity_PV]), " kW installed")
    printx("battery   : ", value.(model[:capacity_storage]), " kWh installed")
    printx("grid capa : ", som(grid.capacity.*value.(model[:grid_capacity]),1), " A")                                    # defined above
    printx("grid cons : ", grid_net, " kWh")                                    # defined above
    printx("comm cons : ", comm_net, " kWh")                                    # defined above
    printx("")
    printx("PV        : ", cost_PV*value.(model[:capacity_PV]), " €")
    printx("battery   : ", cost_storage*value.(model[:capacity_storage]), " €")
    printx("grid capa : ", som(grid.cost.*value.(model[:grid_capacity]),1)*MONTHS/12, " €")                                    # defined above
    printx("grid cons : ", grid_net_cost, " €")                                 # defined above
    printx("comm cons : ", cost_community*comm_net, " €")
    printx("total     : ", value.(model[:costs]), " €")
    # printx("total mean: ", objective_value(model)/COMMUNITY_SIZE, " €")
    printx("total     : ", objective_value(model), " €")
end


function display_results(model, graphs)
    index=range(start=TIME_STEP, step=TIME_STEP, length=steps)
    ymin=minimum(value.(model[:community_exported]))
    ymax=ceil(2*maximum([PV_PROFILE value.(model[:grid_imported]) value.(model[:community_imported])]))/2

    if 1 in graphs
        p1=plot(xlabel = "hours", size=(800,300), ylims=[ymin, ymax])
        for n in eachindex(PROFILES)
            plot!(p1, index, LOAD_PROFILE[:,n], label="$(PROFILES[n]) load")
            plot!(p1, index, PV_PROFILE[:,n], label="$(PROFILES[n]) PV")
        end
        display(p1)
    end

    if 2 in graphs
        for n in eachindex(PROFILES)
            p2=plot(xlabel = "hours", size=(800,300), ylims=[ymin, ymax])
            plot!(p2, index, LOAD_PROFILE[:,n], label="$(PROFILES[n]) load")
            plot!(p2, index, PV_PROFILE[:,n], label="$(PROFILES[n]) PV")
            plot!(p2, index, value.(model[:grid_imported][:,n]), label="$(PROFILES[n]) grid imp")
            plot!(p2, index, value.(model[:grid_exported][:,n]), label="$(PROFILES[n]) grid exp")
            plot!(p2, index, value.(model[:storage_charge][:,n]), label="$(PROFILES[n]) charge")
            plot!(p2, index, value.(model[:storage_discharge][:,n]), label="$(PROFILES[n]) discharge")
            plot!(p2, index, value.(model[:community_imported][:,n]), label="$(PROFILES[n]) community imp")
            plot!(p2, index, value.(model[:community_exported][:,n]), label="$(PROFILES[n]) community exp")
            display(p2)
        end
    end

    if 3 in graphs
        for n in eachindex(PROFILES)
            p3=plot(xlabel = "hours", size=(800,300)) #, ylims=[0, ymax])
            plot!(p3, index, value.(model[:capacity_PV][n])*ones(steps), label="$(PROFILES[n]) PV [kWp]")
            plot!(p3, index, value.(model[:capacity_storage][n])*ones(steps), label="$(PROFILES[n]) battery [kWh]")
            display(p3)
        end
    end

    if 4 in graphs
        for n in eachindex(PROFILES)
            p4=plot(xlabel = "hours", size=(800,300)) #, ylims=[0, ymax])
            plot!(p4, index, value.(model[:storage_charge][:,n]), label="$(PROFILES[n]) charge")
            plot!(p4, index, value.(model[:storage_discharge][:,n]), label="$(PROFILES[n]) discharge")
            plot!(p4, index, value.(model[:storage_SOC][:,n]), label="$(PROFILES[n]) SOC [kWh]")
            plot!(p4, index, value.(model[:capacity_storage][n]).*ones(steps), label="$(PROFILES[n]) battery [kWh]")
            display(p4)
        end
    end

end

# functions to display values without too many digits
using Formatting
function printx(args...)
    for arg in args
        if typeof(arg) == Float64
            _print_float(arg)
        elseif typeof(arg) == Array{Float64} || typeof(arg) == Vector{Float64}
            print("[")
            for ar in eachindex(arg)
                _print_float(arg[ar])
                if ar != length(arg) print("  ") end
            end
            print("]")
        elseif typeof(arg) == Matrix{Float64}
            print("[")
            for ar in eachindex(arg)
                _print_float(arg[ar])
                if ar != length(arg) print("  ") end
            end
            print("]")
        else print(arg)
        end
    end
    println()
end

function _print_float(arg)
    if abs(arg)<1e-5 print("0")
    elseif iszero(arg-floor.(arg)) print(Int.(round.(arg)))
    elseif abs(arg) > 1e3 print(Int.(round.(arg)))
    else
        for i in 2:-1:-8
            if abs(arg) > 10.0^i
                print(sprintf1("%.$(3-i)f",arg))
                break
            end
        end
    end
end

function som(A,dim)
    return dropdims(sum(A,dims=dim),dims=dim)
end
