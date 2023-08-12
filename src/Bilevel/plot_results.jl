# =============================================================================
#                                   Functions
# =============================================================================
# -----------------------------------------------------------------------------
#                               PRINT NETWORK IN TIKZ 
# -----------------------------------------------------------------------------
"""Draws a graph using TikZ."""

function print_initial_network(network, x_scale, y_scale;
                   dir=pwd(),
                   filename="graph",
                   display=false, reshape=false)

    # Let's insert some boilerplate styling
    # and necessary preamble/postamble
    preamble = """\\documentclass{standalone}
    \\usepackage{tikz}
    \\usepackage{amsmath}
    \\usepackage{xcolor}
    \\definecolor{Ulg_blue}{RGB}{9, 111, 123}
    \\usetikzlibrary{graphs, quotes, arrows.meta, positioning}

    \\begin{document}
    \\Large
    \\begin{tikzpicture}[
        every label/.style = {align=center, font=\\normalsize, inner sep=2pt},
        every edge quotes/.style = {font=\\normalsize, text=black, fill=white, inner sep=2pt}
        ]
    \\graph [no placement]
    {\n"""

    
    postamble = """};
    \\end{tikzpicture}
    \\end{document}"""
 
    #style_sub = "rectangle, draw=Ulg_blue, fill=white, minimum size=1.5em, inner sep=1pt"
    style_sub = "rectangle, draw, fill=white, minimum size=1.5em, inner sep=1pt"
    style_load = "circle, draw, fill=white, minimum size=2em, inner sep=1pt"

    file_path = joinpath(dir ,"$filename.tex")

    touch(file_path)

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
            if n <= ns
                write(file,
                "    $n [x=$(x)cm, y=$(y)cm, as={\$\\mathcal{S}_{$n}\$}, $(style_sub)];\n")
            else 
                write(file,
                "    $n [x=$(x)cm, y=$(y)cm, as={\$\\mathcal{U}_{$(n-ns)}\$}, $(style_load)];\n")
            end
        end

        for l in network.lines
            i = l.edge.from_node.id ; j = l.edge.to_node.id
            write(file,
            "    $(i) --[gray, dashed] $(j);\n")
        end  
        write(file, postamble)
    end

    if display
        run(Cmd(`lualatex $filename.tex`, dir="$dir"))
        run(Cmd(`open $filename.pdf`, dir="$dir"))
    else 
        run(Cmd(`lualatex $filename.tex`, dir="$dir"))
        #run(Cmd(`open $filename.pdf`, dir="$dir"))
    end
end

function print_network_tikz(network, time_step, x_scale, y_scale;
                        dir=pwd(), filename="graph", display=false, reshape=false)
    
    
    preamble = """\\documentclass{standalone}
    \\usepackage{tikz}
    \\usepackage{amsmath}
    \\usetikzlibrary{graphs, quotes, arrows.meta, positioning}
    \\usepackage{xcolor}
    \\definecolor{Ulg_blue}{RGB}{9, 111, 123}
    \\definecolor{Ulg_red}{RGB}{132, 38, 0}
    \\definecolor{Ulg_orange}{RGB}{223, 168, 120}
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
 
    style_sub = "rectangle, draw, fill=white, minimum size=1.7em, inner sep=1pt,  align=center, text width=0.6cm"
    style_load = "circle, draw, fill=white, minimum size=1.2em, inner sep=1pt, align=center, text width=0.6cm"
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
            
                if isnothing(b.PV_installation)
                    # if isnothing(b.storage) 
                    write(file,
                    "    $n [x=$(x)cm, y=$(y)cm, as={\$\\scriptsize \$ \\color{Ulg_red} \\mathbf{$P_cons} \$}, $(style_load)];\n")
                    # else 
                    #     P_stor_discharge = b.storage.P_discharge[time_step]
                    #     P_stor_charge = b.storage.P_charge[time_step]
                    #     write(file,
                    #     "    $n [x=$(x)cm, y=$(y)cm, as={\$\\mathcal{U}_{$nu}\$  \\vspace{-0.2em} \\scriptsize  \$ \\color{Ulg_orange} \\mathbf{$P_stor_discharge}/\\mathbf{$P_stor_charge}\$ \\vspace{0.1em} \\scriptsize \$ \\color{Ulg_red} \\mathbf{$P_cons}\$}, $(style_load)];\n")
                    # end
                else 
                    P_gen = b.PV_installation.P[time_step]
                    P_gen = round(P_gen; digits=2)
                   
                    write(file,
                    "    $n [x=$(x)cm, y=$(y)cm, as={\$\\mathcal{U}_{$nu}\$  \\vspace{-0.2em} \\scriptsize  \$ \\color{Ulg_blue} \\mathbf{$P_gen}\$ \\vspace{0.1em} \\scriptsize \$ \\color{Ulg_red} \\mathbf{$P_cons}\$}, $(style_load)];\n")
                    # else 
                    #     P_stor_discharge = round(b.storage.P_discharge[time_step]; digits=2)
                    #     P_stor_charge = round(b.storage.P_charge[time_step]; digits=2)
                    #     write(file,
                    #     "    $n [x=$(x)cm, y=$(y)cm, as={\$\\mathcal{U}_{$nu}\$  \\vspace{-0.2em} \\scriptsize  \$ \\color{Ulg_blue} \\mathbf{$P_gen}\$ \\vspace{0.1em} \\scriptsize\$ \\color{Ulg_orange} \\mathbf{$P_stor_discharge}/\\mathbf{$P_stor_charge}\$ \\vspace{0.1em} \\scriptsize \$ \\color{Ulg_red} \\mathbf{$P_cons}\$}, $(style_load)];\n")
                    # end
                end
            end
        end

        for l in network.lines
            if l.built
                if  l.P_send[time_step] >= 0                   
                    i = l.edge.from_node.id ; j = l.edge.to_node.id
                    p = l.P_send[time_step]
                else 
                    i = l.edge.to_node.id ; j = l.edge.from_node.id
                    p = abs(l.P_rec[time_step])
                end
              
                if l.conductor.name == "Poppy"
                    write(file,
                    "    $(i) ->[\"$(round(p, sigdigits=2))\"] $(j);\n")
                elseif l.conductor.name == "Oxlip"
                    write(file,
                    "    $(i) ->[\"$(round(p, sigdigits=2))\", thick] $(j);\n")
                elseif l.conductor.name == "Daisy"
                    write(file,
                    "    $(i) ->[\"$(round(p, sigdigits=2))\", very thick] $(j);\n")
                elseif l.conductor.name == "Tulip"
                    write(file,
                    "    $(i) ->[\"$(round(p, sigdigits=2))\", ultra thick] $(j);\n")
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
    else
        run(Cmd(`lualatex $filename.tex`, dir="$dir"))
    end
end

# -----------------------------------------------------------------------------
#                               PRINT PLOTS
# -----------------------------------------------------------------------------

function if_small(x)
    for (i, a) in enumerate(x)
        if abs(a) < 1e-6
            x[i] = 0
        end
    end
    return x #round.(x, digits=6)
end

function plot_results(model::JuMP.AbstractModel; dir=pwd(), filename="plot", pgfplot=true)

    # First: 3 nuances of brown, 3 nuances of 
    colors = ["#bb5b46", "#b37e6e", "#bc9780", "#dfa878", "#5f8d8d", "#3f7e84", "#096f7b"]


    network = model[:network_data]
    buses = network.buses
    lines = network.lines 
    conductors = network.conductors
    T = model[:time_steps]
    L = get_nb_lines(network)
    Ns = get_nb_substations(network)
    Nu = get_nb_loads(network)
    K = get_nb_conductors(network)
    BASE_POWER = network.pu_basis.base_power
    BASE_CURRENT = network.pu_basis.base_current

    P_consumed = [buses[Ns + i].load_profile.time_serie[t] * buses[Ns+i].cos_phi * BASE_POWER  for t in 1:T, i in 1:Nu] 
    PV_prod = [buses[Ns + i].PV_installation.profile.time_serie[t] * BASE_POWER  for t in 1:T, i in 1:Nu]

    time = collect(1:T)
    # PV understanding
    for u in 1:Nu
        current_filename = filename * "pv_bus_$(Ns+u)"

        fig = Plots.plot(xlabel="Time Steps", ylabel="Power [MW]", tex_output_standalone = false, legend = :best, tickfontsize=10, formatter=:latex)

        Plots.plot!(fig, time, P_consumed[:, u], label=L"P_{D}", color=colors[7], linewidth=1.5)

        Plots.plot!(fig, time, [value.(model[:p_pv_max])[u] for t in time], label=L"p_{PV,capa}", color=colors[1], linewidth=1.5, linestyle=:dash)

        Plots.plot!(fig, time, value.(model[:p_pv])[:, u], label=L"p_{PV}", linewidth=1.5, color=colors[4])

        Plots.plot!(fig, time, PV_prod[:, u] .* value.(model[:p_pv_max])[u], label=L"p_{PV,max}", linewidth=1.5, color=colors[3])

        Plots.plot!(fig, time, value.(model[:p_exp])[:, u], label=L"p_{exp}", linewidth=1.5, color=colors[2])

        #display(fig)
        #pgfplot && Plots.savefig(fig, joinpath(dir, "$current_filename.tikz"))
        Plots.savefig(fig, joinpath(dir, "$current_filename.pdf"))
    end

    # Plot currents
    for l in 1:L
        current_filename = filename * "current_line_$(l)"
        #pgfplot && Plots.pgfplotsx()
        fig = Plots.plot(xlabel="Time Steps", ylabel="Current [kA]", tex_output_standalone = false, legend = :best, tickfontsize=10,formatter=:latex)
        chosen_cond_index = findfirst(i -> isapprox(i, 1; rtol=1e-4), value.(model[:Alpha])[l, :])
        if isnothing(chosen_cond_index)
            continue
        end

        Plots.plot!(fig, time, network.conductors[chosen_cond_index].max_i * BASE_CURRENT * ones(T), label=L"(\|I\|)^{max}", linewidth=1.5, color=colors[7], linestyle=:dash)

        Plots.plot!(fig, time, sqrt.(if_small(value.(model[:I_sqr])[:, l])) * BASE_CURRENT, label=L"\|I\|", linewidth=1.5, color=colors[1])

        #pgfplot && Plots.savefig(fig, joinpath(dir, "$current_filename.tikz"))
        Plots.savefig(fig, joinpath(dir, "$current_filename.pdf"))
    end

    # Plot currents conductors
    for l in 1:L
        current_filename = filename * "current_k_line_$(l)"
        #pgfplot && Plots.pgfplotsx()
        fig = Plots.plot(xlabel="Time Steps", ylabel="Current [kA]", tex_output_standalone = false, legend = :best, tickfontsize=10,formatter=:latex)
       
        chosen_cond_index = findfirst(i -> isapprox(i, 1; rtol=1e-4), value.(model[:Alpha])[l, :])
        title = isnothing(chosen_cond_index) ? "Line $l is not built" : "Line $l is built with conductor $chosen_cond_index"
        for k in 1:K
            Plots.plot!(fig, time, sqrt.(if_small(value.(model[:I_sqr_k])[:, l, k])) * BASE_CURRENT, label=L"\|I_{ij}\|"*" (k = $k)", title=title, color=colors[k])
        end
        Plots.plot!(fig, time, sum(sqrt.(if_small(value.(model[:I_sqr_k])[:, l, :])) * BASE_CURRENT, dims=2), label=L"\|I_{ij}\|", color=colors[7])
        
        #pgfplot && Plots.savefig(fig, joinpath(dir, "$current_filename.tikz"))
        Plots.savefig(fig, joinpath(dir, "$current_filename.pdf"))
    end

    # Check voltage limits
    for i in 1:(Ns+Nu)
        current_filename = filename * "voltage_bus_$(i)"
        #pgfplot && Plots.pgfplotsx()
        fig = Plots.plot(xlabel="Time Steps", ylabel="Voltage [p.u.]", tex_output_standalone = false, legend = :best, tickfontsize=10,formatter=:latex)

        Plots.plot!(fig, time, buses[i].V_limits.V_min * ones(T), label=L"V^{min}", color=colors[4], linewidth=1.5, linestyle=:dash)
        Plots.plot!(fig, time, buses[i].V_limits.V_max * ones(T), label=L"V^{max}", color=colors[7], linewidth=1.5, linestyle=:dash)
        Plots.plot!(fig, time, sqrt.(value.(model[:V_sqr])[:, i]), label=L"V", color=colors[1], linewidth=1.5)
       

        #pgfplot && Plots.savefig(fig, joinpath(dir, "$current_filename.tikz"))
        Plots.savefig(fig, joinpath(dir, "$current_filename.pdf"))
    end

    # P_substations
    for i in 1:Ns
        current_filename = filename * "substation_$(i)"
        #pgfplot && Plots.pgfplotsx()
        fig = Plots.plot(xlabel="Time Steps", ylabel="Voltage [kV]", tex_output_standalone = false, legend = :best, tickfontsize=10,formatter=:latex)
       
        Plots.plot!(fig, time, value.(model[:S_sub])[:, i], label=L"S_{S}", color=colors[1], linewidth=1.5)
        Plots.plot!(fig, time, value.(model[:P_sub])[:, i], label=L"P_{S}", color=colors[4], linewidth=1.5)
        Plots.plot!(fig, time, value.(model[:S_sub_capa])[i]*ones(T), label=L"S_{S,capa}", color=colors[7], linestyle=:dash, linewidth=1.5)

        #pgfplot && Plots.savefig(fig, joinpath(dir, "$current_filename.tikz"))
        Plots.savefig(fig, joinpath(dir, "$current_filename.pdf"))
    end

    DELTA_T = model[:delta_t]
    # Power Summary
    total_losses_per_t = [(BASE_POWER * value(sum(lines[l].length * conductors[k].r * model[:I_sqr_k][t, l, k] for l in 1:L, k in 1:K))) for t in time] 
    total_losses_per_t_MWh = total_losses_per_t .* DELTA_T/60

    current_filename = filename * "summary_system_power"
    fig = Plots.plot(xlabel="Time Steps", ylabel="Power [MW]", tex_output_standalone = false, legend = :best, tickfontsize=10,formatter=:latex)

    Plots.plot!(fig, time, [sum(value.(model[:P_sub])[t, :]) for t in time] * BASE_POWER, label=L"P_{S,tot}", linewidth=1.5, color=colors[1])
    Plots.plot!(fig, time, [sum(P_consumed[t, :]) * BASE_POWER for t in time], label=L"P_{D,tot}", linewidth=1.5, color=colors[3])
    Plots.plot!(fig, time, [sum(value.(model[:p_pv])[t, :]) * BASE_POWER for t in time], label=L"p_{PV,tot}", linewidth=1.5, color=colors[5])
    Plots.plot!(fig, time, total_losses_per_t, label=L"P_{L,tot}", linewidth=1.5, color=colors[7])
    # plot!(fig, time, [sum(value.(model[:p_pv])[t, :]) * BASE_POWER for t in time] .+ [sum(value.(model[:P_sub])[t, :]) for t in time] * BASE_POWER .- total_losses, label="p_pv + P_sub - losses[MW]")
    Plots.savefig(fig, joinpath(dir, "$current_filename.pdf"))

    # plot!(fig, time, [sum(value.(model[:p_pv])[t, :]) * BASE_POWER for t in time] .+ [sum(value.(model[:P_sub])[t, :]) for t in time] * BASE_POWER .- total_losses, label="p_pv + P_sub - losses[MW]")

    # Energy
    total_losses_per_t_MWh = total_losses_per_t .* DELTA_T/60
    current_filename = filename * "summary_system_energy"    
    fig = Plots.plot(xlabel="Time Steps", ylabel="Energy [MWh]", tex_output_standalone = false, legend = :best, tickfontsize=10,formatter=:latex)
    Plots.plot!(fig, time, [sum(value.(model[:P_sub])[t, :]) for t in time] * BASE_POWER * DELTA_T/60, label=L"E_{S,tot}", linewidth=1.5, color=colors[1])
    Plots.plot!(fig, time, [sum(P_consumed[t, :]) * BASE_POWER * DELTA_T/60 for t in time], label=L"E_{D,tot}", linewidth=1.5, color=colors[3])
    Plots.plot!(fig, time, [sum(value.(model[:p_pv])[t, :]) * BASE_POWER * DELTA_T/60 for t in time], label=L"e_{PV,tot}", linewidth=1.5, color=colors[5])
    Plots.plot!(fig, time, total_losses_per_t_MWh, label=L"E_{L,tot}", linewidth=1.5, color=colors[7])
    Plots.savefig(fig, joinpath(dir, "$current_filename.pdf"))

    # plot!(fig, time, [sum(value.(model[:p_pv])[t, :]) * BASE_POWER for t in time] .+ [sum(value.(model[:P_sub])[t, :]) for t in time] * BASE_POWER .- total_losses, label="p_pv + P_sub - losses[MW]")

    # Plot storage
    if model[:storage]
        for i in 1:Nu
            current_filename = filename * "storage_bus_$(Ns + i)"
            #pgfplot && Plots.pgfplotsx()
            fig = Plots.plot(xlabel="Time Steps", ylabel="Energy [MWh]", tex_output_standalone = false, legend =:topleft, tickfontsize=10,formatter=:latex)

            Plots.plot!(fig, time, value.(model[:p_storage_charge])[:, i] .* DELTA_T/60, label=L"e_{st}^+", linewidth=1.5, color=colors[1])


            Plots.plot!(fig, time, value.(model[:p_storage_discharge])[:, i] .* DELTA_T/60, label=L"e_{st}^-", linewidth=1.5, color=colors[3])

            Plots.plot!(fig, time, value.(model[:p_pv])[:, i] .* DELTA_T/60, label=L"e_{PV}^-", linewidth=1.5, color=colors[4])

            Plots.plot!(fig, time, value.(model[:storage_state])[:, i], label=L"SOC", color=colors[5], linewidth=1.5)

            Plots.hline!(fig, [value(model[:storage_capacity][i])], label=L"e_{st,capa}", color=colors[7], linestyle=:dash, linewidth=1.5)

            # summer_check = round.([value(model[:storage_state][1, i]) value(model[:storage_state][24, i])+ value(model[:p_storage][1, i])* buses[Ns+i].storage.efficiency], sigdigits=2)
            # winter_check = round.([value(model[:storage_state][25, i]) value(model[:storage_state][48, i])+value(model[:p_storage][25, i])* p.STORAGE_EFF], sigdigits=2)
            # Plots.annotate!(12, 0.9 * value(model[:storage_capacity][i]), "$(summer_check[1]) ?= $(summer_check[2])")
            # Plots.annotate!(36, 0.9 * value(model[:storage_capacity][i]), "$(winter_check[1]) ?= $(winter_check[2])")
            # ramping_max : fraction of the battery (dis)charged on 1 hour
            # ramping_max=round(maximum(abs.(diff([value.(model[:storage_state][1:24, i]) value.(model[:storage_state][25:48, i])], dims=2))), sigdigits=2)
            # Plots.annotate!(24, 0.7 * value(model[:storage_capacity][i]), "Ramping max $(ramping_max)/h")
            #pgfplot && Plots.savefig(fig, joinpath(dir, "$current_filename.tikz"))
            Plots.savefig(fig, joinpath(dir, "$current_filename.pdf"))
        end
    end
end

# -----------------------------------------------------------------------------
#                                   KPIs
# -----------------------------------------------------------------------------

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


# -- PRINT CASE DESCRIPTION FUNCTION --
#
function print_case_description(model::JuMP.AbstractModel)
    
    # Fetch data
    network         = model[:network_data]
    DELTA_T         = model[:delta_t]
    NB_PROFILES     = model[:nb_sign_days]
    T               = model[:time_steps]  
    bilevel         = model[:bilevel]
    Ns              = get_nb_substations(network)
    Nu              = get_nb_loads(network)
    N = Ns + Nu
    buses           = network.buses
    BASE_POWER      = network.pu_basis.base_power

    P_consumed = [buses[Ns + i].load_profile.time_serie[t] * buses[Ns + i].cos_phi for t in 1:T, i in 1:Nu]
    
    # 1. Peak load demand
    peak_demand = maximum(sum(P_consumed, dims=2)) * BASE_POWER

    # 2. Energy consumed
    energy_consumed = sum(P_consumed) * BASE_POWER * DELTA_T / 60 * 365 / NB_PROFILES
    
    # 3. P_sub max and min
    peak_sub_P = value.(sum(model[:P_sub], dims=2)) * BASE_POWER
    peak_sub_max = maximum(peak_sub_P) * BASE_POWER
    peak_sub_min = minimum(peak_sub_P) * BASE_POWER

    # Build table
    case_headers = (
        ["Model", "Objective", "Gap", "Solve Time", "Time periods", "Energy consumed", "Peak Demand", "Peak sub max", "Peak sub min"], ["", "kEUR/year", "%", "sec.","", "MWh/year", "MW", "MW", "MW"])
    case_params = sig_round([(bilevel ? "Bilevel" : "Single Level") objective_value(model) relative_gap(model)*100 solve_time(model) T energy_consumed peak_demand peak_sub_max peak_sub_min])
    pretty_table(case_params, header=case_headers)

    # Returns the table as a Matrix
    return vcat([i[j] for i in case_headers, j in 1:length(case_headers[1])], case_params)
end

# -- PRINT COST RESULTS FUNCTION --
#
function print_cost_results(model; latex=false)
    network = model[:network_data]
    DSO_costs = model[:DSO_costs]

    L = get_nb_lines(network)
    K = get_nb_conductors(network)
    lines = network.lines 
    conductors = network.conductors

    # Investements to build new lines in k€
    DSO_fixed_line = value(sum(model[:Alpha][l, k] * lines[l].length * conductors[k].cost for l in 1:L, k in 1:K))

    BASE_POWER = network.pu_basis.base_power
    SUB_COST = DSO_costs.substation
    Ns = get_nb_substations(network)
   
    # Investements to build new substations in k€
    DSO_fixed_sub = value(sum(model[:S_sub_capa][i] * BASE_POWER * SUB_COST for i in 1:Ns))

    # Operational costs of losses on the length of the horizon
    DSO_loss = value(model[:DSO_loss_costs]) * DSO_costs.amortization

    # DSO future value: DSO should get back what he must pay per year for the investment + a certain margin
    DSO_future_value = value(model[:DSO_fixed_costs]) * ((1 + DSO_costs.interest_rate)^DSO_costs.amortization) + value(model[:DSO_loss_costs]) * DSO_costs.amortization

    # DSO revenues
    DSO_revenues = value(sum(model[:grid_costs])) * DSO_costs.amortization
    # Other costs
    user_costs = value(sum(model[:user_costs]))
    PV_costs = value(sum(model[:PV_costs]))
    grid_costs = value(sum(model[:grid_costs]))
    energy_costs = value(sum(model[:energy_costs]))
    energy_revenues = value(sum(model[:energy_revenues]))

    if model[:storage]
        stor_costs = value(sum(model[:storage_costs]))
        kpis_header = (["DSO fixed line", "DSO fixed sub", "DSO loss", "DSO future value", "DSO revenues", "User", "PV", "Storage", "Grid connection", "Energy imported", "Energy exported"],
        ["kEUR", "kEUR", "kEUR", "kEUR", "kEUR", "kEUR/year", "kEUR/year", "kEUR/year", "kEUR/year", "kEUR/year", "kEUR/year"])
        kpis = sig_round([DSO_fixed_line DSO_fixed_sub DSO_loss DSO_future_value DSO_revenues user_costs PV_costs stor_costs grid_costs energy_costs energy_revenues])
    else
        kpis_header = (["DSO fixed line", "DSO fixed sub", "DSO loss", "DSO future value", "DSO revenues", "User", "PV", "Grid connection", "Energy imported", "Energy exported"],
        ["kEUR", "kEUR", "kEUR", "kEUR", "kEUR", "kEUR/year", "kEUR/year", "kEUR/year", "kEUR/year", "kEUR/year"])
        kpis = sig_round([DSO_fixed_line DSO_fixed_sub DSO_loss DSO_future_value DSO_revenues user_costs PV_costs grid_costs energy_costs energy_revenues])
    end
    
    if latex
        pretty_table(kpis, header=kpis_header, backend=Val(:latex))
    else
        pretty_table(kpis, header=kpis_header)
    end
    # returns the table as a Matrix
    return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
end


# -- PRINT NETWORK RESULTS FUNCTION --
#
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
 
    K = get_nb_conductors(network)
    cond_string = ["Nb lines cond "*string(i) for i in 1:K]

    lines_cond = [[i for (i, a) in enumerate(value.(model[:Alpha])[:, j]) if isapprox(a, 1, rtol=1e-2)] for j in 1:K]
    
    kpis_header = ([["Nb lines built", "Nb substations"];cond_string], [["-", "-"];fill("-", 1:K)])
    kpis = reshape(string.([[nb_lines_built, nb_substations];lines_cond]),1,:)
    if latex
        pretty_table(kpis, header=kpis_header, backend=Val(:latex))
    else
        pretty_table(kpis, header=kpis_header)
    end
    # returns the table as a Matrix
    return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
end


# -- PRINT GRID RESULTS FUNCTION --
#
function print_grid_results(model; latex=false)
    network = model[:network_data]
    DELTA_T = model[:delta_t]
    NB_PROFILES = model[:nb_sign_days]
    Nu = get_nb_loads(network)
    Ns = get_nb_substations(network)
    DAYS_A_YEAR = 365
    MULTIPLIER = DELTA_T/60 * DAYS_A_YEAR / NB_PROFILES
    BASE_POWER = network.pu_basis.base_power

    # Average grid capacity
    grid_capacity = (value(sum(model[:s_grid_max])) * BASE_POWER) / Nu

    # Energy exchanges
    energy_from_grid = value(sum(model[:p_imp])) * MULTIPLIER * BASE_POWER
    energy_to_grid = value(sum(model[:p_exp])) * MULTIPLIER * BASE_POWER

    # Substations
    substation_string = ["Substation "*string(i)*" capacity" for i in 1:Ns]
    buses = network.buses
    substation_capacities = [buses[i].S_rating*BASE_POWER for i in 1:Ns]

    # Computation of the losses
    lines = network.lines
    conductors = network.conductors
    T = model[:time_steps]
    L = get_nb_lines(network)
    K = get_nb_conductors(network)

    loss = value(sum(lines[l].length * conductors[k].r *
    model[:I_sqr_k][t, l, k] for t in 1:T, l in 1:L, k in 1:K)) * MULTIPLIER * BASE_POWER

    # Total grid costs
    User_costs = model[:User_costs]
    grid_cost = (value(sum(model[:p_imp])) * User_costs.EI * MULTIPLIER + value(sum(model[:p_imp])) * User_costs.DSOEI * MULTIPLIER + value(sum(model[:s_grid_max])) * User_costs.GCC) * BASE_POWER

    # Total production
    grid_production = value(sum(model[:p_imp])) * MULTIPLIER
    LCOE_grid = grid_cost / grid_production
    kpis_header = ([["Grid capacity", "Energy bought from grid", "Energy sold to grid"];substation_string; ["loss", "LCOE grid"]],
        [["average, MVA", "MWh/year", "MWh/year"];fill("MVA",1:Ns);["MWh/year", "kEUR/MWh"]])
    kpis = reshape(sig_round([[grid_capacity, energy_from_grid, energy_to_grid];substation_capacities;[loss, LCOE_grid]]),1,:)
    if latex
        pretty_table(kpis, header=kpis_header, backend=Val(:latex))
    else
        pretty_table(kpis, header=kpis_header)
    end
    # returns the table as a Matrix
    return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
end

# -- PRINT PV GLOBAL RESULTS FUNCTION --
#
function print_pv_results(model; latex=false)
    DAYS_A_YEAR = 365
    DELTA_T = model[:delta_t]
    NB_PROFILES = model[:nb_sign_days]
    MULTIPLIER = DELTA_T/60 * DAYS_A_YEAR / NB_PROFILES
    network = model[:network_data]
    BASE_POWER = network.pu_basis.base_power
    T = model[:time_steps]
    Nu = get_nb_loads(network)
    Ns = get_nb_substations(network)
    buses = network.buses 

    # ----- 1. PV installed capacity -----
    PV_capa = sum(value.(model[:p_pv_max]))

    # ----- 2. PV penetration -----
    # = amount of PV capacity installed in a grid with respect to the peak load demand
    #  It is commonly defined as either (1) the ratio of the PV installed capacity with respect to the peak load demand, or (2) the ratio of the total energy produced by PV with respect to the total energy consumed.

    P_consumed = [buses[Ns + i].load_profile.time_serie[t] * buses[Ns + i].cos_phi for t in 1:T, i in 1:Nu]
    energy_consumed = sum(P_consumed) * BASE_POWER * MULTIPLIER

    peak_demand = maximum(sum(P_consumed, dims=2)) * BASE_POWER

    # (1) First definition
    PV_penetration_power = PV_capa / peak_demand

    # (2) Second definition
    PV_penetration_energy = sum(value.(model[:p_pv])) * MULTIPLIER/energy_consumed 

    # ----- 3. Energy produced by PV -----
    PV_energy = sum(value.(model[:p_pv])) * MULTIPLIER

    # ----- 4. Global PV potential -----
   
    PV_prod = [(isnothing(buses[Ns + i].PV_installation) ? 0.0 : buses[Ns + i].PV_installation.profile.time_serie[t]) for t in 1:T, i in 1:Nu]

    PV_potential = sum(PV_prod[t, i] * value(model[:p_pv_max][i] * BASE_POWER) for t in 1:T, i in 1:Nu) * MULTIPLIER

    # ----- 5. Mean Self-Sufficiency/Day/user -----

    user_ss = ["Mean Self Sufficiency User $u" for u in 1:Nu]
    user_sc = ["Mean Self Consumption User $u" for u in 1:Nu]
    user_pr = ["Mean Production Ratio User $u" for u in 1:Nu]

    mean_user_self_sufficiency = zeros(Float64, Nu)
    mean_user_self_consumption = zeros(Float64, Nu)
    mean_user_production_ratio = zeros(Float64, Nu)
    for u in 1:Nu    
        nb_non_zeros = 0    
        for t in 1:T
            mean_user_self_sufficiency[u] += (value.(model[:p_pv])[t, u] * DELTA_T/60 - value.(model[:p_exp])[t, u] * DELTA_T/60)/(P_consumed[t, u] * DELTA_T/60) / (model[:nb_sign_days] * T)

            if PV_prod[t, u] * value(model[:p_pv_max][u]) != 0
                mean_user_self_consumption[u] += (value.(model[:p_pv])[t, u] * DELTA_T/60 - value.(model[:p_exp])[t, u] * DELTA_T/60)/(value.(model[:p_pv])[t, u] * DELTA_T/60) / (model[:nb_sign_days])

                mean_user_production_ratio[u] += (value.(model[:p_pv])[t, u] * DELTA_T/60)/(PV_prod[t, u] * value(model[:p_pv_max][u])* DELTA_T/60)/ (model[:nb_sign_days])

                nb_non_zeros += 1
            end
        end
        mean_user_self_consumption[u] /= nb_non_zeros
        mean_user_production_ratio[u] /= nb_non_zeros
    end

    # ----- 5. Global Mean Self-Sufficiency/Day -----
    mean_global_self_sufficiency = sum(mean_user_self_sufficiency ./length(mean_user_self_sufficiency))
    mean_global_self_consumption = sum(mean_user_self_consumption ./length(mean_user_self_consumption))
    mean_global_production_ratio = sum(mean_user_production_ratio ./length(mean_user_production_ratio))
    # for day in 1:model[:nb_sign_days]
    #     idx_first = Int((day - 1) * T / model[:nb_sign_days] + 1)
    #     idx_last = Int(day * T / model[:nb_sign_days])

    #     mean_global_self_suffiency += (sum(value.(model[:p_pv])[idx_first:idx_last, :]) - sum(value.(model[:p_exp])[idx_first:idx_last, :]))/(sum(P_consumed[idx_first:idx_last, :]) / model[:nb_sign_days])

    #     mean_global_self_consumption += (sum(value.(model[:p_pv])[idx_first:idx_last, :]) - sum(value.(model[:p_exp])[idx_first:idx_last, :]))/(sum(value.(model[:p_pv])[idx_first:idx_last, :]) / model[:nb_sign_days])

    #     mean_global_production_ratio += sum(value.(model[:p_pv])[idx_first:idx_last, :])/(sum(PV_prod[idx_first:idx_last, :] .* value.(model[:p_pv_max])) / model[:nb_sign_days])
    # end

    # global_self_sufficiency = (sum(value.(model[:p_pv])) - sum(value.(model[:p_exp])))/sum(P_consumed) 
    # global_self_consumption = (sum(value.(model[:p_pv])) - sum(value.(model[:p_exp])))/sum(value.(model[:p_pv]))

    # ----- 6. Production ratio -----
    #global_production_ratio = sum(value.(model[:p_pv])) / PV_potential

    # ----- 7. PV costs-----
    User_costs = model[:User_costs]
    pv_costs   = value(sum(model[:s_conv_pv])) * BASE_POWER * User_costs.PV_conv / 
            User_costs.amortization_PVC + value(sum(model[:p_pv_max])) * BASE_POWER * User_costs.PV / User_costs.amortization_PV

    # ----- 8. PV production -----
    PV_power_prod = sum(value.(model[:p_pv])) #* MULTIPLIER
    LCOE_pv = pv_costs / PV_energy

    kpis_header = (
        [
            ["PV capacity", "PV penetration (v1)", "PV penetration (v2)", "PV energy", "PV potential", "Mean Global Self sufficiency/day", "Mean Global Self Consumption/day", "Mean Global Production ratio", "LCOE pv"]
            ; user_ss; user_sc; user_pr
        ],
        [["MVA", "%/year", "%/year", "MWh/year", "MWh/year", "average, [%]", "average, [%]", "average, [%]", "kEUR/MWh"]; fill("%",3*Nu)])

    
    kpis = reshape(string.(
            sig_round(
            [
                [PV_capa, PV_penetration_power, PV_penetration_energy, PV_energy, PV_potential, mean_global_self_sufficiency, mean_global_self_consumption, mean_global_production_ratio, LCOE_pv]; 
                mean_user_self_sufficiency; 
                mean_user_self_consumption; 
                mean_user_production_ratio
            ])), 1, :)

    println(kpis)

    if latex
        pretty_table(kpis, header=kpis_header, backend=Val(:latex))
    else
        pretty_table(kpis, header=kpis_header)
    end
    # returns the table as a Matrix
    return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
end

# -- PRINT STORAGE RESULTS FUNCTION --
#
function print_storage_results(model; latex=true)
    network = model[:network_data]
    Nu      = get_nb_loads(network)

    # storage capacity
    storage_capacity = value(sum(model[:storage_capacity])) / Nu

    # Avg storage state
    avg_storage_state = mean(value.(model[:storage_state]))
    kpis_header = (["Storage capacity", "Storage state"],
        ["average, MWh", "average, [-]"])
    kpis = sig_round([storage_capacity avg_storage_state])
    if latex
        pretty_table(kpis, header=kpis_header, backend=Val(:latex))
    else
        pretty_table(kpis, header=kpis_header)
    end
    # returns the table as a Matrix
    return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
end

# -- PRINT SUMMARY RESULTS FUNCTION --
#
function print_summary(model; latex=false)
    DSO_costs = model[:DSO_costs]
    T = model[:time_steps]
    network = model[:network_data]
    Nu = get_nb_loads(network)
    Ns = get_nb_substations(network)
    DELTA_T = model[:delta_t]

    DSO_capex = value(model[:DSO_fixed_costs]) / DSO_costs.amortization / 1000
    DSO_opex = value(model[:DSO_loss_costs])
    PV_costs = value(sum(model[:PV_costs])) / 1000
    grid_costs = value(sum(model[:grid_costs])) / 1000
    net_energy_costs = (value(sum(model[:energy_costs])) - value(sum(model[:energy_revenues]))) / 1000

    buses = network.buses
    P_consumed = [buses[Ns + i].load_profile.time_serie[t] * buses[Ns + i].cos_phi for t in 1:T, i in 1:Nu]
    PV_prod = [(isnothing(buses[Ns + i].PV_installation) ? 0.0 : buses[Ns + i].PV_installation.profile.time_serie[t]) for t in 1:T, i in 1:Nu]


    mean_user_self_sufficiency = zeros(Float64, Nu)
    mean_user_self_consumption = zeros(Float64, Nu)
    for u in 1:Nu   
        nb_non_zeros = 0     
        for t in 1:T
            mean_user_self_sufficiency[u] += (value.(model[:p_pv])[t, u] * DELTA_T/60 - value.(model[:p_exp])[t, u] * DELTA_T/60)/(P_consumed[t, u] * DELTA_T/60) / (model[:nb_sign_days] * T)

            if PV_prod[t, u] * value(model[:p_pv_max][u]) != 0
                mean_user_self_consumption[u] += (value.(model[:p_pv])[t, u] * DELTA_T/60 - value.(model[:p_exp])[t, u] * DELTA_T/60)/(value.(model[:p_pv])[t, u] * DELTA_T/60) / (model[:nb_sign_days])
                nb_non_zeros += 1
            end
        end
        mean_user_self_consumption[u] /= nb_non_zeros
    end

    # ----- 5. Global Mean Self-Sufficiency/Day -----
    self_sufficiency = sum(mean_user_self_sufficiency./length(mean_user_self_sufficiency))
    self_consumption = sum(mean_user_self_consumption ./length(mean_user_self_consumption))
    kpis_header = (["CAPEX", "OPEX", "UPVC", "UGCC", "UNEEC", "USS", "USC"],
        ["M€/y", "k€/y", "M€/y", "M€/y", "M€/y", "%", "%"])
    kpis = sig_round([DSO_capex DSO_opex PV_costs grid_costs net_energy_costs self_sufficiency self_consumption])
    if latex
        pretty_table(kpis, header=kpis_header, backend=Val(:latex))
    else
        pretty_table(kpis, header=kpis_header)
    end
    # returns the table as a Matrix
    return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
end

# -- PRINTED TABLES FUNCTION --
#
function printed_tables(model; print_latex=false)
    if model[:storage]
        return hcat(
            print_case_description(model),
            print_network_results(model, latex=print_latex),
            print_cost_results(model, latex=print_latex),
            print_pv_results(model, latex=print_latex),
            print_grid_results(model, latex=print_latex),
            print_storage_results(model, latex=print_latex),
            print_summary(model, latex=print_latex),
        )
    else 
        return hcat(
            print_case_description(model),
            print_network_results(model, latex=print_latex),
            print_cost_results(model, latex=print_latex),
            print_pv_results(model, latex=print_latex),
            print_grid_results(model, latex=print_latex),
            print_summary(model, latex=print_latex),
        )
    end
end