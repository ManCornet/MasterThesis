# Data from feeder 1, network 1 of LV network models from Manchester (https://www.enwl.co.uk/lvns)

# Topology data
line_ends = Dict(1 => (1, 2), 2 => (2, 3)) # key = line index, value = [node_from, node_to]
Omega_sending = Dict(1 => [1], 2 => [2], 3 => [])
Omega_receiving = Dict(1 => [], 2 => [1], 3 => [2])
N = 3
L = 2 # Number of lines -> added by myself.
line_length = [0.2, 0.3] # km

# Power related data
P_demand = [50, 60, 100] / BASE_POWER # Active power demand at each bus
Q_demand = [5, 6, 10] / BASE_POWER # Reactive power demand at each bus

# Conductors related data
K = 2
# TODO a datastructure to describe a conductor, a list of conductors, then map the Dictionaris below to that list
max_current = Dict(1 => [70 / BASE_CURRENT for l in 1:L],
    2 => [185 / BASE_CURRENT for l in 1:L]) # absolute, pu
conductance = Dict(1 => [2.19 / l / BASE_ADMITANCE for l in line_length],
    2 => [5.16 / l / BASE_ADMITANCE for l in line_length]) # absolute, pu
susceptance = Dict(1 => [(-0.35 / l) / BASE_ADMITANCE for l in line_length],
    2 => [(-2.11 / l) / BASE_ADMITANCE for l in line_length]) # absolute, pu
line_cost = Dict(1 => [l * 1e3 for l in line_length], 2 => [l * 2.5e3 for l in line_length]) # EUR/km

# Substation data
n_s_init = 1 # Should be <= n_s
n_s = 3 # susbstation nodes are numbered to correspond to the first n_s nodes in the network
S_rating_init = [2000, 0, 0] / BASE_POWER
S_rating_max = [2000, 2000, 2000] / BASE_POWER
# substation_cost = [1, 1, 1]*1e0 # EUR small value => shoulf build substations and no line: OK.
substation_cost = [0.01, 1, 1] * 1e8 # EUR
