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
                   display=true, reshape=false)

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
    style_sub = "rectangle, draw, fill=white, minimum size=1.6em, inner sep=1pt"
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
    end
end

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
                if  l.P_send[time_step] >= 0                   
                    i = l.edge.from_node.id ; j = l.edge.to_node.id
                    p = l.P_send[time_step]
                else 
                    i = l.edge.to_node.id ; j = l.edge.from_node.id
                    p = abs(l.P_rec[time_step])
                end
              
                if l.conductor.name == "Poppy"
                    write(file,
                    "    $(i) ->[\"$(round(p, digits=2))\"] $(j);\n")
                elseif l.conductor.name == "Oxlip"
                    write(file,
                    "    $(i) ->[\"$(round(p, digits=2))\", thick] $(j);\n")
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

function plot_results(model::JuMP.AbstractModel; thesis_dir=pwd(), pgfplot=true)

    colors = ["#096F7B", "#FF8D3E","#842600"]

    network = model[:network_data]
    buses = network.buses
    T = model[:time_steps]
    Ns = get_nb_substations(network)
    Nu = get_nb_loads(network)
    BASE_POWER = network.pu_basis.base_power

    P_consumed = [buses[Ns + i].load_profile.time_serie[t] * buses[Ns+i].cos_phi * BASE_POWER  for t in 1:T, i in 1:Nu] 
    PV_prod = [buses[Ns + i].PV_installation.profile.time_serie[t] * BASE_POWER  for t in 1:T, i in 1:Nu]

    for u in 1:Nu
        filename = "user@bus$(Ns+u)"
        pgfplot && Plots.pgfplotsx()
        fig = Plots.plot(title="Bus $(u)", ylabel="Power [MW]", tex_output_standalone = false, legend = :topleft)

        plot!(fig, 1:T, P_consumed[:, u], label="P_cons", color=colors[i])
        plot!(fig, 1:T, [value.(model[:p_pv_max])[u] for t in 1:T], label="P_pv_max")
        plot!(fig, 1:T, value.(model[:p_pv])[:, u], label="P_pv")
        plot!(fig, 1:T, PV_prod[:, u] * value.(model[:p_pv_max])[u], label="PV_prod max")
        plot!(fig, 1:T, value.(model[:p_exp])[:, u], label="P_exp")

        pgfplot && Plots.savefig(fig, joinpath(thesis_dir, "$filename.tikz"))
        Plots.savefig(fig, joinpath(thesis_dir, "$filename.pdf"))
    end


#     # print network
#     for t in 1:T
#         Bilevel.print_network_tikz(model[:network_data], t, 4.5, 4.5; dir=thesis_dir, filename="network@timestep_$t")
#     end
 

#     # Plot currents
#     for l in 1:L
#         filename = "current@line$(l)"
#         pgfplot && Plots.pgfplotsx()
#         fig = Plots.plot(title="Line $(l)", ylabel="Power [MW]", tex_output_standalone = false, legend = :topleft)

#         chosen_cond_index = findfirst(i -> isapprox(i, 1; rtol=1e-4), value.(model[:Alpha])[l, :])
#         if isnothing(chosen_cond_index)
#             continue
#         end
#         plot!(p.T, sqrt.(if_small(value.(model[:I_sqr])[l, :])) * p.BASE_CURRENT, label="I [kA]", title="Line $(l)")
#         plot!(plt, p.T, p.MAX_CURRENT[l, chosen_cond_index] * p.BASE_CURRENT * ones(p.T_SIZE), label="MAX_CURRENT [kA]")
#         display(plt)
#     end

#     # Display the power for each conductor
#     for l in p.L
#         plt = Plots.plot()
#         chosen_cond_index = findfirst(i -> isapprox(i, 1; rtol=1e-4), value.(model[:Alpha])[l, :])
#         title = isnothing(chosen_cond_index) ? "Line $l is not built" : "Line $l is built with conductor $chosen_cond_index"
#         for k in p.K
#             plot!(plt, p.T, sqrt.(if_small(value.(model[:I_sqr_k])[l, k, :])) * p.BASE_CURRENT, label="I_ij_$k [kA]", title=title)
#         end
#         plot!(plt, p.T, sum(sqrt.(if_small(value.(model[:I_sqr_k])[l, :, :])) * p.BASE_CURRENT, dims=1)', label="I_ij [kA]")
#         display(plt)
#     end


#     # Check voltage limits
#     for i in p.N
#         plt = Plots.plot()
#         plot!(p.T, sqrt.(value.(model[:V_sqr])[i, :]), label="V [p.u.]", title="Bus $(i)")
#         plot!(plt, p.T, p.MIN_VOLTAGE * ones(p.T_SIZE), label="MIN_VOLTAGE [p.u.]")
#         plot!(plt, p.T, p.MAX_VOLTAGE * ones(p.T_SIZE), label="MAX_VOLTAGE [p.u.]")
#         display(plt)
#     end

#     for i in p.NS
#         plt = Plots.plot()
#         plot!(p.T, value.(model[:P_sub])[i, :], label="P_sub", title="Bus $(i)")
#         display(plt)
#     end

#     # Compute total losses per time step: 
#     total_losses = [sum(p.R[:, :] .* value.(model[:I_sqr_k])[:, :, t] .* p.BASE_POWER) for t in p.T]
#     # Check P_sub, Q_sub
#     plt = Plots.plot()
#     plot!(p.T, [sum(value.(model[:P_sub])[:, t]) for t in p.T] * p.BASE_POWER, label="P_sub [MW]")
#     plot!(plt, p.T, [sum(p.P_CONSUMPTION[:, t]) * p.BASE_POWER for t in p.T], label="Total P_CONSUMPTION [MW]")
#     plot!(plt, p.T, [sum(value.(model[:p_pv])[:, t]) * p.BASE_POWER for t in p.T], label="Total p_pv [MW]")
#     plot!(plt, p.T, total_losses, label="Total Joule Losses [MW]")
#     plot!(plt, p.T, [sum(value.(model[:p_pv])[:, t]) * p.BASE_POWER for t in p.T] .+ [sum(value.(model[:P_sub])[:, t]) for t in p.T] * p.BASE_POWER .- total_losses, label="p_pv + P_sub - losses[MW]")
#     display(plt)
# end

# function plot_storage(p, model)
#     # Plot storage power transfers
#     for i in p.NU
#         plt = Plots.plot(legend=:bottomleft)
#         plot!(p.T, value.(model[:p_storage])[i, :], label="storage power [MW]", title="Bus $(i)")
#         plot!(p.T, value.(model[:storage_state])[i, :], label="storage state [MWh]")
#         hline!([value(model[:storage_capacity][i])], label="storage capacity [MWh]")
#         summer_check = round.([value(model[:storage_state][i, 1]) value(model[:storage_state][i, 24])+value(model[:p_storage][i, 1])* p.STORAGE_EFF], sigdigits=2)
#         winter_check = round.([value(model[:storage_state][i, 25]) value(model[:storage_state][i, 48])+value(model[:p_storage][i, 25])* p.STORAGE_EFF], sigdigits=2)
#         annotate!(12, 0.9 * value(model[:storage_capacity][i]), "$(summer_check[1]) ?= $(summer_check[2])")
#         annotate!(36, 0.9 * value(model[:storage_capacity][i]), "$(winter_check[1]) ?= $(winter_check[2])")
#         # ramping_max : fraction of the battery (dis)charged on 1 hour
#         ramping_max=round(maximum(abs.(diff([value.(model[:storage_state][i, 1:24]) value.(model[:storage_state][i, 25:48])], dims=1))), sigdigits=2)
#         annotate!(24, 0.7 * value(model[:storage_capacity][i]), "Ramping max $(ramping_max)/h")
#         display(plt)
#     end
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
    DELTA_T         = model[:delta_t]
    NB_PROFILES     = model[:nb_sign_days]
    T   = model[:time_steps]  
    bilevel         = model[:bilevel]
    Ns              = get_nb_substations(network)
    Nu              = get_nb_loads(network)
    N = Ns + Nu
    buses           = network.buses
    BASE_POWER      = network.pu_basis.base_power

    P_consumed = [buses[Ns + i].load_profile.time_serie[t]* buses[Ns + i].cos_phi for t in 1:T, i in 1:Nu]
    
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

# seems okay: to verify !!
function print_storage_results(model; latex=true)
    network = model[:network_data]
    Nu = get_nb_loads(network)

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



# NOT OKAY !!!!! TO WORK
# function print_pv_results(p, model; latex=true)
#     DAYS_A_YEAR = 365
#     DELTA_T = model[:delta_t]
#     NB_PROFILES = model[:nb_sign_days]
#     MULTIPLIER = DELTA_T/60 * DAYS_A_YEAR / NB_PROFILES

#     # PV penetration
#     PV_penetration = sum(value.(model[:p_pv_max]))

#     # PV energy
#     PV_energy = sum(value.(model[:p_pv])) * MULTIPLIER / PV_penetration # percentage

#     # PV potential
#     network = model[:network_data]
#     Nu = get_nb_loads(network)
#     Ns = get_nb_substations(network)
#     T  = model[:time_steps]
#     buses = network.buses 
#     PV_prod = [(isnothing(buses[Ns + i].PV_installation) ? 0.0 : buses[Ns + i].PV_installation.profile.time_serie[t]) for t in 1:T, i in 1:Nu]

#     PV_potential = sum(PV_prod[t, i] * value(model[:p_pv_max][i]) for i in 1:Nu, t in 1:T) * MULTIPLIER / PV_penetration

#     # Self-sufficiency
#     P_consumed = [buses[i].load_profile.time_serie[t] * buses[Ns + i].cos_phi for t in 1:T, i in 1:Nu]

#     # Self-suffiency and self-consumption ? => ASK BERTRAND
#     avg_self_sufficiency = 0
#     avg_self_consumption = 0

#     mean(value.(model[:p_pv]) .- value.(model[:p_exp]) ./ P_consumed)

#     for t in 1:T
#         avg_self_sufficiency += (value(sum(model[:p_pv][t, :])) - value(sum(model[:p_exp][t, :])))/(sum(P_consumed[t, :]))/T

#         if value(model[:p_pv][t, i]) != 0
#             avg_self_consumption += (value(sum(model[:p_pv][t, :])) - value(sum(model[:p_exp][t, :])))/value(model[:p_pv][t, i])/T
#         end
#     end
    
#     index_non_zeros = PV_potential .!= 0
#     production_ratio = mean(value.(model[:p_pv][p.NU, :])[index_non_zeros] ./ PV_potential[index_non_zeros])

#     potential_array = PV_prod .* value.(model[:p_pv_max])
#     index_non_zeros = potential_array .!= 0

#     production_ratio = mean(value.(model[:p_pv])[index_non_zeros] ./ potential_array[index_non_zeros])

#     # PV_costs 
#     User_costs = model[:User_costs]
#     pv_cost = value(sum(model[:s_conv_pv])) * User_costs.PV_conv/ User_costs.amortization_PVC + value(sum(model[:p_pv_max])) * User_costs.PV / User_costs.amortization_PV

#     # Pv prod
#     pv_production = MULTIPLIER * sum(value.(model[:p_pv]))
#     LCOE_pv = pv_cost / pv_production
#     kpis_header = (["PV penetration", "PV production", "PV potential", "Avg self sufficiency", "Avg self consumption", "production ratio", "LCOE pv"],
#         ["MVA", "MWh/MVA/year", "MWh/MVA/year", "average, [-]", "average, [-]", "average, [-]", "kEUR/MWh"])
#     kpis = sig_round([PV_penetration PV_energy PV_potential avg_self_sufficiency self_consumption production_ratio LCOE_pv])
#     if latex
#         pretty_table(kpis, header=kpis_header, backend=Val(:latex))
#     else
#         pretty_table(kpis, header=kpis_header)
#     end
#     # returns the table as a Matrix
#     return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
# end
