@variables(LOWER_MODEL, begin
    p_imp[N, T] >= 0, (container = Array) # active power imported at time t
    q_imp[N, T] >= 0, (container = Array)    # reactive power imported at time t
    p_exp[N, T] >= 0, (container = Array) # active power exported at time t
    q_exp[N, T] >= 0, (container = Array)   # reactive power exported at time t
    s_grid_max[N] >= 0, (container = Array) # grid imp/exp capacity
    p_pv[N, T] >= 0, (container = Array) # active power generated at time t
    #q_pv[N, T], (container = Array)   # reactive power generated at time t # DEACTIVATED in a first time to reimplement and test later.
    s_conv_pv[N] >= 0, (container = Array) # PV converter max capacity
    0 <= p_pv_max[N] <= MAX_PV_CAPACITY_PER_NODE, (container = Array) # PV capacity
    user_costs[N], (container = Array) # Total cost of each user, i.e. PV_costs + energy_costs - energy_revenues + grid_costs
    PV_costs[N] >= 0, (container = Array) # Costs related to PV investment for each user
    energy_costs[N] >= 0, (container = Array)  # Costs related to energy imported, for each user
    energy_revenues[N] >= 0, (container = Array)  # Costs related to energy imported, for each user
    grid_costs[N] >= 0, (container = Array) # Costs related to network for each user, both capacy and energy-based
end)