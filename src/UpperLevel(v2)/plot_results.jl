#-----------------------------------------------------------------------------
#
#                           - TFE : Bilevel DNEP - 
#                             University of Liege
#
#-----------------------------------------------------------------------------
# Created By  : Manon Cornet
# Created Date: Saturday May 20 2023
#
# plot_results:
#   plot results
#
# =============================================================================
#                                   Functions
# =============================================================================
# -----------------------------------------------------------------------------
#                               PRINT NETWORK IN TIKZ 
# -----------------------------------------------------------------------------
function print_network_tikz(network, time_step, x_scale, y_scale;
                        dir=pwd(), filename="graph", display=true, reshape=false)
    
    
    preamble = """\\documentclass{standalone}
    \\usepackage{tikz}
    \\usepackage{amsmath}
    \\usetikzlibrary{graphs, quotes, arrows.meta, positioning}
    \\usepackage{xcolor}
    \\definecolor{Ulg_blue}{RGB}{9, 111, 123}
    \\definecolor{Ulg_red}{RGB}{132, 38, 0}
    \\begin{document}
    \\begin{tikzpicture}[
        every label/.style = {align=center, font=\\tiny, inner sep=2pt},
        every edge quotes/.style = {font=\\scriptsize, text=black, fill=white, inner sep=2pt}
        ]
    \\graph [no placement]
    {\n"""

    postamble = """};
    \\end{tikzpicture}
    \\end{document}"""
 
    style_sub = "rectangle, draw, fill=white, minimum size=1em, inner sep=1pt,  align=center, text width=0.6cm"
    style_load = "circle, draw, fill=white, minimum size=1em, inner sep=1pt, align=center, text width=0.6cm"
    base_power = network.pu_basis.base_power

    # Creating the .tex file
    file_path = joinpath(dir ,"$filename.tex")
    touch(file_path)

    # Writing the .tex file
    open(file_path, "w") do file
        write(file, preamble)
        
        for b in network.buses
            node = b.node
            n = node.id
            x_coord = node.coord.x
            y_coord = node.coord.y
            if reshape
                if n in [1, 2]
                    x_coord += 0.5
                elseif n == 16
                    y_coord += 0.1
                elseif n == 20
                    y_coord -=  0.1
                elseif n in [3, 9]
                    if n == 9 
                        x_coord += 0.5  
                    end
                    y_coord -= 1.8
                elseif n in [4, 5]
                    if n == 4
                        x_coord += 0.8
                    else 
                        x_coord += 0.5
                    end
                    y_coord += 0.4
                end
            end
            x = x_scale * x_coord
            y = y_scale * y_coord
            
            ns = get_nb_substations(network)
            nu = n - ns
            if n <= ns
                P_gen = round((b.P_sup[time_step]); digits = 2)
                col = P_gen >= 0 ? "\\color{Ulg_blue}" : "\\color{Ulg_red}"
                if b.built
                    write(file,
                    "    $n [x=$(x)cm, y=$(y)cm, as={\$\\mathcal{S}_{$n}\$  \\vspace{0.1em}  \\scriptsize \$\\mathbf{$col $P_gen}\$}, $(style_sub), very thick];\n")
                else 
                    write(file,
                    "    $n [x=$(x)cm, y=$(y)cm, as={\$\\mathcal{S}_{$n}\$ \\vspace{0.1em}  \\scriptsize \$\\mathbf{$col $P_gen}\$}, $(style_sub)];\n")
                end
            else 
                P_cons = round(b.load_profile.time_serie[time_step] * b.cos_phi; digits = 2)
                P_gen = isnothing(b.PV_installation) ? 0.0 : b.PV_installation.P[time_step]
                P_gen = round(P_gen; digits=2)
                
                write(file,
                "    $n [x=$(x)cm, y=$(y)cm, as={\$\\mathcal{U}_{$nu}\$  \\vspace{-0.2em} \\scriptsize  \$ \\color{Ulg_blue} \\mathbf{$P_gen}\$ \\vspace{0.1em} \\scriptsize \$ \\color{Ulg_red} \\mathbf{$P_cons}\$}, $(style_load)];\n")
            end
        end

        for l in network.lines
            if l.built
                p_s = l.P_send[time_step]
                p_r = l.P_rec[time_step]
                if  p_s >= 0                   
                    i = l.edge.from_node.id ; j = l.edge.to_node.id
                    p = l.P_send[time_step]
                else 
                    i = l.edge.to_node.id ; j = l.edge.from_node.id
                    p = abs(l.P_rec[time_step])
                end
              
                if l.conductor.name == "Poppy"
                    write(file,
                    "    $(i) ->[\"$(round(p_s, digits=2)), $(round(p_r, digits=2))\"] $(j);\n")
                elseif l.conductor.name == "Oxlip"
                    write(file,
                    "    $(i) ->[\"$(round(p_s, digits=2)), $(round(p_r, digits=2))\", thick] $(j);\n")
                elseif l.conductor.name == "Daisy"
                    write(file,
                    "    $(i) ->[\"$(round(p, digits=2))\", very thick] $(j);\n")
                elseif l.conductor.name == "Tulip"
                    write(file,
                    "    $(i) ->[\"$(round(p, digits=2))\", ultra thick] $(j);\n")
                end
                print(l.conductor.name)
            else 
                i = l.edge.from_node.id ; j = l.edge.to_node.id
                write(file,
                "    $(i) --[gray, dashed] $(j);\n")
            end       
        end  
        write(file, postamble)
    end

    if display
        run(Cmd(`lualatex $filename.tex`, dir="$dir"))
        run(Cmd(`open $filename.pdf`, dir="$dir"))
    end
end

# -----------------------------------------------------------------------------
#                               PLOT NETWORK RESULTS 
# -----------------------------------------------------------------------------

# -- FUNCTION THAT PRINTS THE NETWORK --
function print_network(
    p, Alpha, P_ij, I_sqr_k, p_pv, P_sub, P_CONSUMPTION)

    chosen_lines = sum(Alpha, dims=2)
    edge_label_dict = Dict{Tuple{Int64,Int64},Float64}()

    # Find values of line powers

    line_power = [(P_ij[l] >= 0) ? P_ij[l] : abs(P_ij[l] - sum(p.R[l, k] * I_sqr_k[l] for k in p.K)) for l in p.L]

    # Creating the graph
    edge_width = Dict()
    g = SimpleDiGraph(p.N_SIZE)
    for l in p.L
        if P_ij[l] < 0
            add_edge!(g, p.LINE_ENDS[l][2], p.LINE_ENDS[l][1])
            key = (p.LINE_ENDS[l][2], p.LINE_ENDS[l][1])
        else
            add_edge!(g, p.LINE_ENDS[l][1], p.LINE_ENDS[l][2])
            key = p.LINE_ENDS[l]
        end
        if isapprox(chosen_lines[l], 1; rtol=1e-4) # TO BE MODIFIED
            edge_label_dict[key] = abs(round(line_power[l]; digits=4) .* p.BASE_POWER)
        end
        edge_width[key] = (isapprox(chosen_lines[l], 1; rtol=1e-4) ? 1 : 0.1)
    end

    # -- SHAPES OF THE NODES OF THE NETWORK --
    house_nodeshape(x_i, y_i, s) = [(x_i + 0.7s * dx, y_i + 0.7s * dy) for (dx, dy) in [(1, 1), (0, 1.6), (-1, 1), (-1, -1), (1, -1), (1, 1)]]
    subs_nodeshape(x_i, y_i, s) = [(x_i + 1.2s * dx, y_i + 1.2s * dy) for (dx, dy) in [(1, 1), (-1, 1), (-1, -1), (1, -1), (1, 1)]]


    # colors = [colorant"lightseagreen", colorant"orange"]
    colors = ["#689BAA", "#C2C5DB"]
    house_labels = ["+$(round(value.(p_pv)[n]*p.BASE_POWER, digits=3)) \n\n -$(round(P_CONSUMPTION[n]*p.BASE_POWER, digits=3))" for n in p.NU]
    subs_labels = ["+$(round(value.(P_sub)[n]*p.BASE_POWER, digits=3))" for n in p.NS]
    node_labels = [subs_labels; house_labels]
    node_shapes = [[subs_nodeshape for _ in p.NS]; [house_nodeshape for _ in p.NU]]

    # See components.jl from Plots.jl to have all the arguments
    graph = graphplot(adjacency_matrix(g),
        x=[[p.NODES_LOCATION[1][n] + 0.5 for n in p.NS]
            [p.NODES_LOCATION[1][n] for n in p.NU]],             # x-coordinate of the nodes
        y=p.NODES_LOCATION[2],                                # y-coordinate of the nodes
        nodesize=0.1,
        nodestrokewidth=0,                                        # coutour line width of the node
        edgestyle=:solid,
        nodealpha=1,                                        # transparency of node color
        names=node_labels,                              # node label
        nodeshape=node_shapes,                              # :circle, :ellipse, :hexagon
        nodecolor=colors[[[1 for _ in p.NS]; [2 for _ in p.NU]]],
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

function plot_results(p, model)
    # PV understanding
    for u in p.NU
        plt = Plots.plot()
        plot!(p.T, p.P_CONSUMPTION[u, :], label="P_conso", title="Bus $(u)")
        #plot!(plt, T, Q_CONSUMPTION[u, :], label="Q_conso")
        plot!(plt, p.T, [value.(model[:p_pv_max])[u] for t in p.T], label="pv max")
        plot!(plt, p.T, value.(model[:p_pv])[u, :], label="p_pv")
        plot!(plt, p.T, p.PV_PRODUCTION[u, :] * value.(model[:p_pv_max])[u], label="max PV prod")
        #plot!(plt, T, value.(model[:q_pv])[u, :], label="q_pv")
        plot!(plt, p.T, value.(model[:p_exp])[u, :], label="p_exp")
        display(plt)
    end


    #using PyPlot
    for t in p.T
        pyplot()
        network_plot = print_network(
            p,
            value.(model[:Alpha]),
            value.(model[:P_ij])[:, t],
            value.(model[:I_sqr_k])[:, :, t],
            value.(model[:p_pv])[:, t],
            value.(model[:P_sub])[:, t],
            p.P_CONSUMPTION[:, t]
        )
        display(network_plot)
        #PyPlot.savefig("network_state_time_step_$(t).png")
    end

    # # log record
    # dict_input = Dict(
    #     "S_CONSUMPTION" => Containers.DenseAxisArray(p.S_CONSUMPTION, p.N, p.T),
    #     "P_CONSUMPTION" => Containers.DenseAxisArray(p.P_CONSUMPTION, p.N, p.T),
    #     "Q_CONSUMPTION" => Containers.DenseAxisArray(p.Q_CONSUMPTION, p.N, p.T),
    #     "PV_PRODUCTION" => Containers.DenseAxisArray(p.PV_PRODUCTION, p.N, p.T),
    #     "MAX_CURRENT" => Containers.DenseAxisArray(p.MAX_CURRENT, p.L, p.K),
    #     "MAX_VOLTAGE" => Containers.DenseAxisArray([p.MAX_VOLTAGE], 1:1),
    # )
    # dict_output = Dict(String(k) => value.(v) for (k, v) in object_dictionary(model))
    # data = merge(dict_output, dict_input)

    # xlsx_output(data)

    # Plot currents
    for l in p.L
        plt = Plots.plot()
        chosen_cond_index = findfirst(i -> isapprox(i, 1; rtol=1e-4), value.(model[:Alpha])[l, :])
        if isnothing(chosen_cond_index)
            continue
        end
        plot!(p.T, sqrt.(if_small(value.(model[:I_sqr])[l, :])) * p.BASE_CURRENT, label="I [kA]", title="Line $(l)")
        plot!(plt, p.T, p.MAX_CURRENT[l, chosen_cond_index] * p.BASE_CURRENT * ones(p.T_SIZE), label="MAX_CURRENT [kA]")
        display(plt)
    end

    # Display the power for each conductor
    for l in p.L
        plt = Plots.plot()
        chosen_cond_index = findfirst(i -> isapprox(i, 1; rtol=1e-4), value.(model[:Alpha])[l, :])
        title = isnothing(chosen_cond_index) ? "Line $l is not built" : "Line $l is built with conductor $chosen_cond_index"
        for k in p.K
            plot!(plt, p.T, sqrt.(if_small(value.(model[:I_sqr_k])[l, k, :])) * p.BASE_CURRENT, label="I_ij_$k [kA]", title=title)
        end
        plot!(plt, p.T, sum(sqrt.(if_small(value.(model[:I_sqr_k])[l, :, :])) * p.BASE_CURRENT, dims=1)', label="I_ij [kA]")
        display(plt)
    end


    # Check voltage limits
    for i in p.N
        plt = Plots.plot()
        plot!(p.T, sqrt.(value.(model[:V_sqr])[i, :]), label="V [p.u.]", title="Bus $(i)")
        plot!(plt, p.T, p.MIN_VOLTAGE * ones(p.T_SIZE), label="MIN_VOLTAGE [p.u.]")
        plot!(plt, p.T, p.MAX_VOLTAGE * ones(p.T_SIZE), label="MAX_VOLTAGE [p.u.]")
        display(plt)
    end

    for i in p.NS
        plt = Plots.plot()
        plot!(p.T, value.(model[:P_sub])[i, :], label="P_sub", title="Bus $(i)")
        display(plt)
    end

    # Compute total losses per time step: 
    total_losses = [sum(p.R[:, :] .* value.(model[:I_sqr_k])[:, :, t] .* p.BASE_POWER) for t in p.T]
    # Check P_sub, Q_sub
    plt = Plots.plot()
    plot!(p.T, [sum(value.(model[:P_sub])[:, t]) for t in p.T] * p.BASE_POWER, label="P_sub [MW]")
    plot!(plt, p.T, [sum(p.P_CONSUMPTION[:, t]) * p.BASE_POWER for t in p.T], label="Total P_CONSUMPTION [MW]")
    plot!(plt, p.T, [sum(value.(model[:p_pv])[:, t]) * p.BASE_POWER for t in p.T], label="Total p_pv [MW]")
    plot!(plt, p.T, total_losses, label="Total Joule Losses [MW]")
    plot!(plt, p.T, [sum(value.(model[:p_pv])[:, t]) * p.BASE_POWER for t in p.T] .+ [sum(value.(model[:P_sub])[:, t]) for t in p.T] * p.BASE_POWER .- total_losses, label="p_pv + P_sub - losses[MW]")
    display(plt)
end

function plot_storage(p, model)
    # Plot storage power transfers
    for i in p.NU
        plt = Plots.plot(legend=:bottomleft)
        plot!(p.T, value.(model[:p_storage])[i, :], label="storage power [MW]", title="Bus $(i)")
        plot!(p.T, value.(model[:storage_state])[i, :], label="storage state [MWh]")
        hline!([value(model[:storage_capacity][i])], label="storage capacity [MWh]")
        summer_check = round.([value(model[:storage_state][i, 1]) value(model[:storage_state][i, 24])+value(model[:p_storage][i, 1])* p.STORAGE_EFF], sigdigits=2)
        winter_check = round.([value(model[:storage_state][i, 25]) value(model[:storage_state][i, 48])+value(model[:p_storage][i, 25])* p.STORAGE_EFF], sigdigits=2)
        annotate!(12, 0.9 * value(model[:storage_capacity][i]), "$(summer_check[1]) ?= $(summer_check[2])")
        annotate!(36, 0.9 * value(model[:storage_capacity][i]), "$(winter_check[1]) ?= $(winter_check[2])")
        # ramping_max : fraction of the battery (dis)charged on 1 hour
        ramping_max=round(maximum(abs.(diff([value.(model[:storage_state][i, 1:24]) value.(model[:storage_state][i, 25:48])], dims=1))), sigdigits=2)
        annotate!(24, 0.7 * value(model[:storage_capacity][i]), "Ramping max $(ramping_max)/h")
        display(plt)
    end
end

"""Compute KPIs."""

function sig_round(x)
    for (i, a) in enumerate(x)
        if isa(a, AbstractFloat)
            x[i] = round(a, sigdigits=4)
            if abs(a) < 1e-3
                x[i] = 0
            end
        end
    end
    return x
end

# This one seems okay: but check with results of Geoffrey
function print_case_description(model::JuMP.AbstractModel)
    network = model[:network_data]
    delta_t = model[:delta_t]
    nb_sign_days = model[:nb_sign_days]
    nb_time_steps = model[:time_steps]  
    bilevel = model[:bilevel]
    N = get_nb_buses(network)
    Ns = get_nb_substations(network)
    buses = network.buses
    energy_consumed = 0
    for i in Ns+1:N
        P_cons = sum(buses[i].load_profile.time_serie .* buses[i].cos_phi)
        energy_consumed += P_cons * delta_t * 365 / nb_sign_days
    end
    
    peak_sub_P = value.(sum(model[:P_sub], dims=2))
    peak_sub_max = maximum(peak_sub_P)
    peak_sub_min = minimum(peak_sub_P)
    case_headers = (["Model", "Objective", "Gap", "Solve Time", "Time periods", "Energy consumed", "Peak sub max", "Peak sub min"], ["", "kEUR/year", "%", "sec.","", "MWh/year", "MW on 1 DT", "MW on 1 DT"])
    case_params = sig_round([(bilevel ? "Bilevel" : "Single Level") objective_value(model) relative_gap(model)*100 solve_time(model) nb_time_steps energy_consumed peak_sub_max peak_sub_min])
    pretty_table(case_params, header=case_headers)
    # returns the table as a Matrix
    return vcat([i[j] for i in case_headers, j in 1:length(case_headers[1])], case_params)
end

# This one seems okay: but check with results of Geoffrey
function print_network_results(model; latex=false)
    network = model[:network_data]
    nb_lines_built = 0
    for l in network.lines 
        if l.built 
            nb_lines_built += 1
        end
    end

    Ns = get_nb_substations(network)
    buses = network.buses
    nb_substations = 0
    for i in 1:Ns
        if buses[i].built
            nb_substations += 1
        end
    end
 
    lines_cond1 = [[i for (i, a) in enumerate(value.(model[:Alpha])[:, 1]) if isapprox(a, 1, rtol=1e-2)]]
    lines_cond2 = [[i for (i, a) in enumerate(value.(model[:Alpha])[:, 2]) if isapprox(a, 1, rtol=1e-2)]]
    lines_cond3 = [[i for (i, a) in enumerate(value.(model[:Alpha])[:, 3]) if isapprox(a, 1, rtol=1e-2)]]
    lines_cond4 = [[i for (i, a) in enumerate(value.(model[:Alpha])[:, 4]) if isapprox(a, 1, rtol=1e-2)]]
    kpis_header = (["Nb lines built", "Nb substations", "Nb lines cond 1", "Nb lines cond 2", "Nb lines cond 3", "Nb lines cond 4"], ["-", "-", "-", "-", "-", "-"])
    kpis = [nb_lines_built nb_substations lines_cond1 lines_cond2 lines_cond3 lines_cond4]
    [kpis[i] = string(a) for (i, a) in enumerate(kpis) if (i >= 3) && (!isempty(a))]
    if latex
        pretty_table(kpis, header=kpis_header, backend=Val(:latex))
    else
        pretty_table(kpis, header=kpis_header)
    end
    # returns the table as a Matrix
    return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
end

# function print_cost_results(p, model; latex=true)
#     DSO_fixed_line = value(sum(sum(model[:Alpha] .* p.LINE_COST)))
#     DSO_fixed_sub = value(sum(model[:S_sub_capa]) * p.SUB_INSTALL_COST)
#     DSO_loss = value(model[:DSO_loss_costs])* p.AMORTIZATION_DSO
#     DSO_future_value = (value(model[:DSO_fixed_costs]) * (1 + p.DSO_INTEREST_RATE)^p.AMORTIZATION_DSO) + value(model[:DSO_loss_costs]) * p.AMORTIZATION_DSO
#     DSO_revenues = value(sum(model[:grid_costs])) * p.AMORTIZATION_DSO
#     user_costs = value(sum(model[:user_costs]))
#     PV_costs = value(sum(model[:PV_costs]))
#     grid_costs = value(sum(model[:grid_costs]))
#     energy_costs = value(sum(model[:energy_costs]))
#     energy_revenues = value(sum(model[:energy_revenues]))
#     kpis_header = (["DSO fixed line", "DSO fixed sub", "DSO loss", "DSO future value", "DSO revenues", "User", "PV", "Grid connection", "Energy imported", "Energy exported"],
#         ["kEUR", "kEUR", "kEUR", "kEUR", "kEUR", "kEUR/year", "kEUR/year", "kEUR/year", "kEUR/year", "kEUR/year"])
#     kpis = sig_round([DSO_fixed_line DSO_fixed_sub DSO_loss DSO_future_value DSO_revenues user_costs PV_costs grid_costs energy_costs energy_revenues])
#     if latex
#         pretty_table(kpis, header=kpis_header, backend=Val(:latex))
#     else
#         pretty_table(kpis, header=kpis_header)
#     end
#     # returns the table as a Matrix
#     return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
# end

# function print_case_results(p, model; latex=true)
#     PV_penetration = sum(value.(model[:p_pv_max]))
#     PV_energy = sum(value.(model[:p_pv])) * p.TIME_STEP * p.DAYS_A_YEAR / p.NB_PROFILES / PV_penetration
#     PV_potential = sum(p.PV_PRODUCTION[i, t] * value(model[:p_pv_max][i]) for i in p.NU, t in p.T) * p.TIME_STEP * p.DAYS_A_YEAR / p.NB_PROFILES / PV_penetration
#     energy_to_grid = value(sum(model[:p_exp])) * p.TIME_STEP * p.DAYS_A_YEAR / p.NB_PROFILES
#     energy_from_grid = value(sum(model[:p_imp])) * p.TIME_STEP * p.DAYS_A_YEAR / p.NB_PROFILES
#     substation1_capacity = value(p.S_RATING_INIT[1] + model[:S_sub_capa][1])
#     substation2_capacity = value(p.S_RATING_INIT[2] + model[:S_sub_capa][2])
#     kpis_header = (["PV penetration", "PV production", "PV potential", "Energy bought from grid", "Energy sold to grid", "Substation1 capacity", "Substation2 capacity"], ["MVA", "MWh/MVA/year", "MWh/MVA/year", "MWh/year", "MWh/year", "MVA", "MVA"])
#     kpis = sig_round([PV_penetration PV_energy PV_potential energy_from_grid energy_to_grid substation1_capacity substation2_capacity])
#     if latex
#         pretty_table(kpis, header=kpis_header, backend=Val(:latex))
#     else
#         pretty_table(kpis, header=kpis_header)
#     end
#     # returns the table as a Matrix
#     return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
# end

# function print_case_processed_results(p, model; latex=true)
#     self_sufficiency = mean(min.(value.(model[:p_pv])[p.NU, :], p.P_CONSUMPTION[p.NU, :]) ./ p.P_CONSUMPTION[p.NU, :])
#     PV_potential = (p.PV_PRODUCTION.*value.(model[:p_pv_max]))[p.NU, :]
#     index_non_zeros = PV_potential .!= 0
#     self_consumption = mean(min.(PV_potential, p.P_CONSUMPTION[p.NU, :])[index_non_zeros] ./
#                             PV_potential[index_non_zeros])
#     production_ratio = mean(value.(model[:p_pv][p.NU, :])[index_non_zeros] ./
#                             PV_potential[index_non_zeros])
#     loss = value(sum(p.R .* sum(model[:I_sqr_k], dims=3))) * p.DAYS_A_YEAR * p.TIME_STEP / p.NB_PROFILES
#     pv_cost = value(sum(model[:s_conv_pv])) * p.CONVERTER_COST / p.AMORTIZATION_CONVERTER + value(sum(model[:p_pv_max])) * p.PV_COST / p.AMORTIZATION_PV
#     pv_production = sum(value.(model[:p_pv])) * p.TIME_STEP * p.DAYS_A_YEAR / p.NB_PROFILES
#     LCOE_pv = pv_cost / pv_production
#     grid_cost = (value(sum(model[:p_imp])) * p.IMP_ELECTRICITY_ENRG_COST) * p.TIME_STEP * p.DAYS_A_YEAR / p.NB_PROFILES + #=- value(sum(model[:p_exp])) * EXP_ELECTRICITY_ENRG_COST=#
#                 (value(sum(model[:p_imp])) * p.IMP_ELECTRICITY_DSO_COST) * p.TIME_STEP * p.DAYS_A_YEAR / p.NB_PROFILES + #=+ value(sum(model[:p_exp])) * EXP_ELECTRICITY_DSO_COST=#
#                 value(sum(model[:s_grid_max])) * p.GRID_CONNECTION_COST
#     grid_production = (value(sum(model[:p_imp]))) * p.TIME_STEP * p.DAYS_A_YEAR / p.NB_PROFILES #=- sum(model[:p_exp])=#
#     LCOE_grid = grid_cost / grid_production
#     kpis_header = (["self sufficiency", "self consumption", "production ratio", "loss", "LCOE pv", "LCOE grid"], ["average, [-]", "average, [-]", "average, [-]", "MWh/year", "kEUR/MWh", "kEUR/MWh"])
#     kpis = sig_round([self_sufficiency self_consumption production_ratio loss LCOE_pv LCOE_grid])
#     if latex
#         pretty_table(kpis, header=kpis_header, backend=Val(:latex))
#     else
#         pretty_table(kpis, header=kpis_header)
#     end
#     # returns the table as a Matrix
#     return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
# end

# function print_summary(p, model; latex=true)
#     DSO_capex = value(model[:DSO_fixed_costs]) / p.AMORTIZATION_DSO / 1000
#     DSO_opex = value(model[:DSO_loss_costs]) / 1000
#     PV_costs = value(sum(model[:PV_costs])) / 1000
#     grid_costs = value(sum(model[:grid_costs])) / 1000
#     net_energy_costs = (value(sum(model[:energy_costs])) - value(sum(model[:energy_revenues]))) / 1000
#     self_sufficiency = mean(min.(value.(model[:p_pv])[p.NU, :], p.P_CONSUMPTION[p.NU, :]) ./ p.P_CONSUMPTION[p.NU, :]) *100
#     PV_potential = (p.PV_PRODUCTION.*value.(model[:p_pv_max]))[p.NU, :]
#     index_non_zeros = PV_potential .!= 0
#     self_consumption = mean(min.(PV_potential, p.P_CONSUMPTION[p.NU, :])[index_non_zeros] ./
#                             PV_potential[index_non_zeros]) *100
#     kpis_header = (["CAPEX", "OPEX", "UPVC", "UGCC", "UNEEC", "USS", "USC"],
#         ["M€/y", "M€/y", "M€/y", "M€/y", "M€/y", "%", "%"])
#     kpis = sig_round([DSO_capex DSO_opex PV_costs grid_costs net_energy_costs self_sufficiency self_consumption])
#     if latex
#         pretty_table(kpis, header=kpis_header, backend=Val(:latex))
#     else
#         pretty_table(kpis, header=kpis_header)
#     end
#     # returns the table as a Matrix
#     return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
# end


# function printed_tables(p, model; print_latex=false)
#     return hcat(
#         print_case_description(p, model),
#         print_network_results(p, model, latex=print_latex),
#         print_cost_results(p, model, latex=print_latex),
#         print_case_results(p, model, latex=print_latex),
#         print_case_processed_results(p, model, latex=print_latex),
#         print_summary(p, model, latex=print_latex),
#     )
# end

