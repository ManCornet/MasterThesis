# Print the main results 
using PrettyTables, Latexify

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

function print_case_description()
    energy_consumed = sum(P_CONSUMPTION) * TIME_STEP * DAYS_A_YEAR / NB_PROFILES
    peak_sub_P = value.(sum(model[:P_sub], dims=1))
    peak_sub_max = maximum(peak_sub_P)
    peak_sub_min = minimum(peak_sub_P)
    case_headers = (["Model", "Objective", "Time periods", "Energy consumed", "Peak sub max", "Peak sub min"], ["", "kEUR/year", "", "MWh/year", "MW on 1 DT", "MW on 1 DT"])
    case_params = sig_round([(SINGLE_LEVEL ? "Single Level" : "Bilevel") objective_value(model) length(T) energy_consumed peak_sub_max peak_sub_min])
    pretty_table(case_params, header=case_headers)
    # returns the table as a Matrix
    return vcat([i[j] for i in case_headers, j in 1:length(case_headers[1])], case_params)
end

function print_network_results(; latex=true)
    lines_built = sum(sum(value.(model[:Alpha]), dims=2))
    nb_substations = length(Ns_init) + sum(value.(model[:Beta])[Ns_notinit])
    lines_cond1 = [[i for (i,a) in enumerate(value.(model[:Alpha])[:,1]) if isapprox(a, 1, rtol=1e-2)]]
    lines_cond2 = [[i for (i,a) in enumerate(value.(model[:Alpha])[:,2]) if isapprox(a, 1, rtol=1e-2)]]
    lines_cond3 = [[i for (i,a) in enumerate(value.(model[:Alpha])[:,3]) if isapprox(a, 1, rtol=1e-2)]]
    lines_cond4 = [[i for (i,a) in enumerate(value.(model[:Alpha])[:,4]) if isapprox(a, 1, rtol=1e-2)]]
    kpis_header = (["Nb lines built", "Nb substations", "Nb lines cond 1", "Nb lines cond 2", "Nb lines cond 3", "Nb lines cond 4"], ["-", "-", "-", "-", "-", "-"])
    kpis = [lines_built nb_substations lines_cond1 lines_cond2 lines_cond3 lines_cond4]
    if latex
        pretty_table(kpis, header=kpis_header, backend=Val(:latex))
    else
        pretty_table(kpis, header=kpis_header)
    end
    # returns the table as a Matrix
    return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
end

function print_cost_results(; latex=true)
    DSO_fixed = value(model[:DSO_fixed_costs])
    DSO_loss = value(model[:DSO_loss_costs]) * AMORTIZATION_DSO
    DSO_future_value = (value(model[:DSO_fixed_costs]) * (1 + DSO_INTEREST_RATE)^AMORTIZATION_DSO) + value(model[:DSO_loss_costs]) * AMORTIZATION_DSO
    DSO_revenues = value(sum(model[:grid_costs])) * AMORTIZATION_DSO
    user_costs = value(sum(model[:user_costs]))
    PV_costs = value(sum(model[:PV_costs]))
    grid_costs = value(sum(model[:grid_costs]))
    energy_costs = value(sum(model[:energy_costs]))
    energy_revenues = value(sum(model[:energy_revenues]))
    kpis_header = (["DSO fixed", "DSO loss", "DSO future value", "DSO revenues", "User",      "PV",        "Grid connection", "Energy imported", "Energy exported"],
                   ["kEUR",      "kEUR",     "kEUR",              "kEUR",        "kEUR/year", "kEUR/year", "kEUR/year",       "kEUR/year",       "kEUR/year"])
    kpis = sig_round([DSO_fixed DSO_loss DSO_future_value DSO_revenues user_costs PV_costs grid_costs energy_costs energy_revenues])
    if latex
        pretty_table(kpis, header=kpis_header, backend=Val(:latex))
    else
        pretty_table(kpis, header=kpis_header)
    end
    # returns the table as a Matrix
    return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
end

function print_case_results(; latex=true)
    PV_penetration = sum(value.(model[:p_pv_max]))
    PV_energy = sum(value.(model[:p_pv])) * TIME_STEP * DAYS_A_YEAR / NB_PROFILES / PV_penetration
    PV_potential = sum(PV_PRODUCTION[i, t] * value(model[:p_pv_max][i]) for i in Nu, t in T) * TIME_STEP * DAYS_A_YEAR / NB_PROFILES / PV_penetration
    energy_to_grid = value(sum(model[:p_exp])) * TIME_STEP * DAYS_A_YEAR / NB_PROFILES
    energy_from_grid = value(sum(model[:p_imp])) * TIME_STEP * DAYS_A_YEAR / NB_PROFILES
    substation_capacity = value(sum(S_rating_init + model[:S_sub_capa]))
    kpis_header = (["PV penetration", "PV production", "PV potential", "Energy bought from grid", "Energy sold to grid", "Substation capacity"], ["MVA", "MWh/MVA/year", "MWh/MVA/year", "MWh/year", "MWh/year", "MVA"])
    kpis = sig_round([PV_penetration PV_energy PV_potential energy_from_grid energy_to_grid substation_capacity])
    if latex
        pretty_table(kpis, header=kpis_header, backend=Val(:latex))
    else
        pretty_table(kpis, header=kpis_header)
    end
    # returns the table as a Matrix
    return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
end

function print_case_processed_results(; latex=true)
    self_sufficiency = mean(min.(value.(model[:p_pv])[Nu, :], P_CONSUMPTION[Nu, :]) ./ P_CONSUMPTION[Nu, :])
    PV_potential = (PV_PRODUCTION.*value.(model[:p_pv_max]))[Nu, :]
    index_non_zeros = PV_potential .!= 0
    self_consumption = mean(min.(PV_potential, P_CONSUMPTION[Nu, :])[index_non_zeros] ./
                            PV_potential[index_non_zeros])
    loss = value(sum(R .* sum(model[:I_sqr_k], dims=3))) * DAYS_A_YEAR * TIME_STEP / NB_PROFILES
    pv_cost = value(sum(model[:s_conv_pv])) * CONVERTER_COST / AMORTIZATION_CONVERTER + value(sum(model[:p_pv_max])) * PV_COST / AMORTIZATION_PV
    pv_production = sum(value.(model[:p_pv])) * TIME_STEP * DAYS_A_YEAR / NB_PROFILES
    LCOE_pv = pv_cost / pv_production
    grid_cost = (value(sum(model[:p_imp])) * IMP_ELECTRICITY_ENRG_COST - value(sum(model[:p_exp])) * EXP_ELECTRICITY_ENRG_COST) * TIME_STEP * DAYS_A_YEAR / NB_PROFILES +
                (value(sum(model[:p_imp])) * IMP_ELECTRICITY_DSO_COST + value(sum(model[:p_exp])) * EXP_ELECTRICITY_DSO_COST) * TIME_STEP * DAYS_A_YEAR / NB_PROFILES +
                value(sum(model[:s_grid_max])) * GRID_CONNECTION_COST
    grid_production = (value(sum(model[:p_imp]) - sum(model[:p_exp]))) * TIME_STEP * DAYS_A_YEAR / NB_PROFILES
    LCOE_grid = grid_cost / grid_production
    kpis_header = (["self sufficiency", "self consumption", "loss", "LCOE pv", "LCOE grid"], ["average, [-]", "average, [-]", "MWh/year", "kEUR/MWh", "kEUR/MWh"])
    kpis = sig_round([self_sufficiency self_consumption loss LCOE_pv LCOE_grid])
    if latex
        pretty_table(kpis, header=kpis_header, backend=Val(:latex))
    else
        pretty_table(kpis, header=kpis_header)
    end
    # returns the table as a Matrix
    return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
end

print_latex = false
printed_tables = hcat(
    print_case_description(),
    print_network_results(latex=print_latex),
    print_cost_results(latex=print_latex),
    print_case_results(latex=print_latex),
    print_case_processed_results(latex=print_latex),
);
