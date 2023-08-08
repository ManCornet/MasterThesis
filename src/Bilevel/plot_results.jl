"""Draws a graph using TikZ."""
# Put reshape = true for Nahman Peric graph
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
                    "    $(i) ->[\"$(round(p_s, digits=2))\"] $(j);\n")
                elseif l.conductor.name == "Oxlip"
                    write(file,
                    "    $(i) ->[\"$(round(p_s, digits=2))\", thick] $(j);\n")
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

#This one seems okay: but check with results of Geoffrey
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
    user_costs = value(sum(model[:user_costs]))
    PV_costs = value(sum(model[:PV_costs]))
    grid_costs = value(sum(model[:grid_costs]))
    energy_costs = value(sum(model[:energy_costs]))
    energy_revenues = value(sum(model[:energy_revenues]))
    kpis_header = (["DSO fixed line", "DSO fixed sub", "DSO loss", "DSO future value", "DSO revenues", "User", "PV", "Grid connection", "Energy imported", "Energy exported"],
        ["kEUR", "kEUR", "kEUR", "kEUR", "kEUR", "kEUR/year", "kEUR/year", "kEUR/year", "kEUR/year", "kEUR/year"])
    kpis = sig_round([DSO_fixed_line DSO_fixed_sub DSO_loss DSO_future_value DSO_revenues user_costs PV_costs grid_costs energy_costs energy_revenues])
    if latex
        pretty_table(kpis, header=kpis_header, backend=Val(:latex))
    else
        pretty_table(kpis, header=kpis_header)
    end
    # returns the table as a Matrix
    return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
end

# seems ok
function print_case_results(model; latex=false)
    DAYS_A_YEAR = 365
    DELTA_T = model[:delta_t]
    NB_PROFILES = model[:nb_sign_days]
   
    MULTIPLIER = DELTA_T/60 * DAYS_A_YEAR / NB_PROFILES
    # PV penetration : p_pv max = peak PV 
    PV_penetration = sum(value.(model[:p_pv_max]))
    # PV energy
    PV_energy = sum(value.(model[:p_pv])) * MULTIPLIER / PV_penetration

    # PV potential
    network = model[:network_data]
    Nu = get_nb_loads(network)
    Ns = get_nb_substations(network)
    T  = model[:time_steps]
    buses = network.buses 
    PV_prod = [(isnothing(buses[Ns + i].PV_installation) ? zeros(Float64, T) : buses[Ns + i].PV_installation.profile.time_serie) for i in 1:Nu]

    PV_potential = sum(PV_prod[i][t] * value(model[:p_pv_max][i]) for i in 1:Nu, t in 1:T) * MULTIPLIER / PV_penetration

    # Energy exchanges with grid
    energy_to_grid = value(sum(model[:p_exp])) * MULTIPLIER
    energy_from_grid = value(sum(model[:p_imp])) * MULTIPLIER

    # Substations capa
    sub_init_capa = [buses[i].S_rating for i in 1:Ns]
    substation1_capacity = sub_init_capa[1]
    substation2_capacity = sub_init_capa[2]
    kpis_header = (["PV penetration", "PV production", "PV potential", "Energy bought from grid", "Energy sold to grid", "Substation1 capacity", "Substation2 capacity"], ["MVA", "MWh/MVA/year", "MWh/MVA/year", "MWh/year", "MWh/year", "MVA", "MVA"])
    kpis = sig_round([PV_penetration PV_energy PV_potential energy_from_grid energy_to_grid substation1_capacity substation2_capacity])
    if latex
        pretty_table(kpis, header=kpis_header, backend=Val(:latex))
    else
        pretty_table(kpis, header=kpis_header)
    end
    # returns the table as a Matrix
    return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
end