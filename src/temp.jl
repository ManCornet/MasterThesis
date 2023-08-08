# # Seems ok: weird results for grid costs etc
# function print_cost_results(model; latex=true)
#     network = model[:network_data]
#     DSO_costs = model[:DSO_costs]


#     L = get_nb_lines(network)
#     K = get_nb_conductors(network)
#     lines = network.lines 
#     conductors = network.conductors

#     # Investements to build new lines in k€
#     DSO_fixed_line = value(sum(model[:Alpha][l, k] * lines[l].length * conductors[k].cost for l in 1:L, k in 1:K))

#     BASE_POWER = network.pu_basis.base_power
#     SUB_COST = DSO_costs.substation
#     Ns = get_nb_substations(network)
   
#     # Investements to build new substations in k€
#     DSO_fixed_sub = value(sum(model[:S_sub_capa][i] * BASE_POWER * SUB_COST for i in 1:Ns))
   
#     # Operational costs of losses on the length of the horizon
#     DSO_loss = value(model[:DSO_loss_costs]) * DSO_costs.amortization

#     # DSO future value: DSO should get back what he must pay per year for the investment + a certain margin
#     DSO_future_value = value(model[:DSO_fixed_costs]) * ((1 + DSO_costs.interest_rate)^DSO_costs.amortization) + value(model[:DSO_loss_costs]) * DSO_costs.amortization

#     # DSO revenues
#     DSO_revenues = value(sum(model[:grid_costs])) * DSO_costs.amortization
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

# Seems ok
# function print_case_results(model; latex=false)
#     DAYS_A_YEAR = 365
#     DELTA_T = model[:delta_t]
#     NB_PROFILES = model[:nb_sign_days]
   
#     MULTIPLIER = DELTA_T/60 * DAYS_A_YEAR / NB_PROFILES
#     # PV penetration : p_pv max = peak PV 
#     PV_penetration = sum(value.(model[:p_pv_max]))
#     # PV energy
#     PV_energy = sum(value.(model[:p_pv])) * MULTIPLIER / PV_penetration

#     # PV potential
#     network = model[:network_data]
#     Nu = get_nb_loads(network)
#     Ns = get_nb_substations(network)
#     T  = model[:time_steps]
#     buses = network.buses 
#     PV_prod = [(isnothing(buses[Ns + i].PV_installation) ? zeros(Float64, T) : buses[Ns + i].PV_installation.profile.time_serie) for i in 1:Nu]

#     PV_potential = sum(PV_prod[i][t] * value(model[:p_pv_max][i]) for i in 1:Nu, t in 1:T) * MULTIPLIER / PV_penetration

#     # Energy exchanges with grid
#     energy_to_grid = value(sum(model[:p_exp])) * MULTIPLIER
#     energy_from_grid = value(sum(model[:p_imp])) * MULTIPLIER

#     # Substations capa
#     sub_init_capa = [buses[i].S_rating for i in 1:Ns]
#     substation1_capacity = sub_init_capa[1]
#     substation2_capacity = sub_init_capa[2]
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
