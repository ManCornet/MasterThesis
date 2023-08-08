function print_cost_results(model; latex=true)
    DSO_fixed_line = value(sum(sum(model[:Alpha] .* p.LINE_COST)))
    DSO_fixed_sub = value(sum(model[:S_sub_capa]) * p.SUB_INSTALL_COST)
    DSO_loss = value(model[:DSO_loss_costs])* p.AMORTIZATION_DSO
    DSO_future_value = (value(model[:DSO_fixed_costs]) * (1 + p.DSO_INTEREST_RATE)^p.AMORTIZATION_DSO) + value(model[:DSO_loss_costs]) * p.AMORTIZATION_DSO
    DSO_revenues = value(sum(model[:grid_costs])) * p.AMORTIZATION_DSO
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

function print_case_results(p, model; latex=true)
    PV_penetration = sum(value.(model[:p_pv_max]))
    PV_energy = sum(value.(model[:p_pv])) * p.TIME_STEP * p.DAYS_A_YEAR / p.NB_PROFILES / PV_penetration
    PV_potential = sum(p.PV_PRODUCTION[i, t] * value(model[:p_pv_max][i]) for i in p.NU, t in p.T) * p.TIME_STEP * p.DAYS_A_YEAR / p.NB_PROFILES / PV_penetration
    energy_to_grid = value(sum(model[:p_exp])) * p.TIME_STEP * p.DAYS_A_YEAR / p.NB_PROFILES
    energy_from_grid = value(sum(model[:p_imp])) * p.TIME_STEP * p.DAYS_A_YEAR / p.NB_PROFILES
    substation1_capacity = value(p.S_RATING_INIT[1] + model[:S_sub_capa][1])
    substation2_capacity = value(p.S_RATING_INIT[2] + model[:S_sub_capa][2])
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

function print_case_processed_results(p, model; latex=true)
    self_sufficiency = mean(min.(value.(model[:p_pv])[p.NU, :], p.P_CONSUMPTION[p.NU, :]) ./ p.P_CONSUMPTION[p.NU, :])
    PV_potential = (p.PV_PRODUCTION.*value.(model[:p_pv_max]))[p.NU, :]
    index_non_zeros = PV_potential .!= 0
    self_consumption = mean(min.(PV_potential, p.P_CONSUMPTION[p.NU, :])[index_non_zeros] ./
                            PV_potential[index_non_zeros])
    production_ratio = mean(value.(model[:p_pv][p.NU, :])[index_non_zeros] ./
                            PV_potential[index_non_zeros])
    loss = value(sum(p.R .* sum(model[:I_sqr_k], dims=3))) * p.DAYS_A_YEAR * p.TIME_STEP / p.NB_PROFILES
    pv_cost = value(sum(model[:s_conv_pv])) * p.CONVERTER_COST / p.AMORTIZATION_CONVERTER + value(sum(model[:p_pv_max])) * p.PV_COST / p.AMORTIZATION_PV
    pv_production = sum(value.(model[:p_pv])) * p.TIME_STEP * p.DAYS_A_YEAR / p.NB_PROFILES
    LCOE_pv = pv_cost / pv_production
    grid_cost = (value(sum(model[:p_imp])) * p.IMP_ELECTRICITY_ENRG_COST) * p.TIME_STEP * p.DAYS_A_YEAR / p.NB_PROFILES + #=- value(sum(model[:p_exp])) * EXP_ELECTRICITY_ENRG_COST=#
                (value(sum(model[:p_imp])) * p.IMP_ELECTRICITY_DSO_COST) * p.TIME_STEP * p.DAYS_A_YEAR / p.NB_PROFILES + #=+ value(sum(model[:p_exp])) * EXP_ELECTRICITY_DSO_COST=#
                value(sum(model[:s_grid_max])) * p.GRID_CONNECTION_COST
    grid_production = (value(sum(model[:p_imp]))) * p.TIME_STEP * p.DAYS_A_YEAR / p.NB_PROFILES #=- sum(model[:p_exp])=#
    LCOE_grid = grid_cost / grid_production
    kpis_header = (["self sufficiency", "self consumption", "production ratio", "loss", "LCOE pv", "LCOE grid"], ["average, [-]", "average, [-]", "average, [-]", "MWh/year", "kEUR/MWh", "kEUR/MWh"])
    kpis = sig_round([self_sufficiency self_consumption production_ratio loss LCOE_pv LCOE_grid])
    if latex
        pretty_table(kpis, header=kpis_header, backend=Val(:latex))
    else
        pretty_table(kpis, header=kpis_header)
    end
    # returns the table as a Matrix
    return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
end

function print_summary(p, model; latex=true)
    DSO_capex = value(model[:DSO_fixed_costs]) / p.AMORTIZATION_DSO / 1000
    DSO_opex = value(model[:DSO_loss_costs]) / 1000
    PV_costs = value(sum(model[:PV_costs])) / 1000
    grid_costs = value(sum(model[:grid_costs])) / 1000
    net_energy_costs = (value(sum(model[:energy_costs])) - value(sum(model[:energy_revenues]))) / 1000
    self_sufficiency = mean(min.(value.(model[:p_pv])[p.NU, :], p.P_CONSUMPTION[p.NU, :]) ./ p.P_CONSUMPTION[p.NU, :]) *100
    PV_potential = (p.PV_PRODUCTION.*value.(model[:p_pv_max]))[p.NU, :]
    index_non_zeros = PV_potential .!= 0
    self_consumption = mean(min.(PV_potential, p.P_CONSUMPTION[p.NU, :])[index_non_zeros] ./
                            PV_potential[index_non_zeros]) *100
    kpis_header = (["CAPEX", "OPEX", "UPVC", "UGCC", "UNEEC", "USS", "USC"],
        ["M€/y", "M€/y", "M€/y", "M€/y", "M€/y", "%", "%"])
    kpis = sig_round([DSO_capex DSO_opex PV_costs grid_costs net_energy_costs self_sufficiency self_consumption])
    if latex
        pretty_table(kpis, header=kpis_header, backend=Val(:latex))
    else
        pretty_table(kpis, header=kpis_header)
    end
    # returns the table as a Matrix
    return vcat([i[j] for i in kpis_header, j in 1:length(kpis_header[1])], kpis)
end


function printed_tables(p, model; print_latex=false)
    return hcat(
        print_case_description(p, model),
        print_network_results(p, model, latex=print_latex),
        print_cost_results(p, model, latex=print_latex),
        print_case_results(p, model, latex=print_latex),
        print_case_processed_results(p, model, latex=print_latex),
        print_summary(p, model, latex=print_latex),
    )
end
