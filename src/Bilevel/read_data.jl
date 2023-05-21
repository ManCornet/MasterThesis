#-----------------------------------------------------------------------------
#
#                           - TFE : Bilevel DNEP - 
#                             University of Liege
#
#-----------------------------------------------------------------------------
# Created By  : Manon Cornet
# Created Date: Saturday March 20 2023
#
# utils:
#   File containing the function required to read the data required to launch 
#   a simulation of the bilevel model
#
# =============================================================================
#                                   Imports
# =============================================================================
import DataFrames
import Random 
import XLSX

# =============================================================================
#                                  Functions
# =============================================================================
""" define_pu_basis 

    Arguments:
    ----------
        - BASE_POWER   : value of the base power in MW
        - BASE_VOLTAGE : value of the base voltage in kV
        - BASE_MONEY   : value of the base money in kA
    
    Return value:
    -------------
        - A named tuple containing the pu basis
"""
function define_pu_basis(;BASE_POWER::Float64 = 1.0, 
                        BASE_VOLTAGE::Float64 = 34.5, 
                        BASE_MONEY = 1.0
                        )
    
    BASE_CURRENT = BASE_POWER / BASE_VOLTAGE      # [kA]
    BASE_IMPEDANCE = BASE_VOLTAGE / BASE_CURRENT  # [Ohm]

    pu_basis = (base_power=BASE_POWER, 
                base_voltage=BASE_VOLTAGE,
                base_current=BASE_CURRENT, 
                base_impedance=BASE_IMPEDANCE,
                base_money=BASE_MONEY
                )
                
    return pu_basis 
end

""" get_buses_data 

    Arguments:
    ----------
        - NETWORK_PATH   : path of xlsx file contrainint the network data
        - voltage_limits : limits on the bus voltage
        - pu_basis       : named tuple containing the per-unit basis used for the simulation
    
    Return value:
    -------------
        - A structure of type DistributionNetwork that contains the network data
"""
function get_buses_data(NETWORK_PATH::String, V_limits::VLIM, pu_basis::PU_BASIS)
    # Get the sheet corresponding to network buses data
    df_bus = DataFrames.DataFrame(XLSX.readtable(NETWORK_PATH, "bus"))

    N = DataFrames.nrow(df_bus)            # Total number of buses
    Ns = sum(df_bus.type .== "substation") # Number of substation buses
    Nu = N - Ns                            # Number of user buses

    nodes = [Node(i, (x=convert(Float64, df_bus.x[i]), y=convert(Float64, df_bus.y[i]))) for i in 1:N]
    subs_buses = [Substation(nodes[i], V_limits, convert(Float64, df_bus.S_G_max_mva[i]) / pu_basis.base_power) for i in 1:Ns]
    load_buses = [User(nodes[i], V_limits) for i in Ns+1:N]

    return nodes, subs_buses, load_buses
end

""" get_lines_data 

    Arguments:
    ----------
        - NETWORK_PATH : path of xlsx file containing the network data
        - nodes        : list containing all the nodes of the network graph
    
    Return values:
    -------------
        - edges: list containing all the edges of the network graph
        - lines: list containing all the electrical lines of the network
"""
function get_lines_data(NETWORK_PATH::String, nodes::Vector{Node})
    # Get the sheet corresponding to network lines data
    df_line = DataFrames.DataFrame(XLSX.readtable(NETWORK_PATH, "line"))

    L = DataFrames.nrow(df_line) # number of lines
    # Get all the edges of the network
    edges = [Edge(l, nodes[df_line.from_bus[l]], nodes[df_line.to_bus[l]]) for l in 1:L]
    # Get all the lines
    lines = [Line(edges[l], convert(Float64, df_line.length_km[l])) for l in 1:L]
    
    return edges, lines
end

""" get_conductors_data 

    Arguments:
    ----------
        - NETWORK_PATH : path of xlsx file containing the network data
        - pu_basis     : named tuple containing the per-unit basis used for the simulation
    
    Return values:
    -------------
        - conds: list of available conductors to build a line
"""
function get_conductors_data(NETWORK_PATH::String, pu_basis::PU_BASIS)
    # Get the sheet corresponding to network buses data
    df_cond = DataFrames.DataFrame(XLSX.readtable(NETWORK_PATH, "conductor"))

    K = DataFrames.nrow(df_cond)    

    conds = [Conductor( df_cond.name[k], 
                        convert(Float64, df_cond.r_ohm_per_km[k]) / pu_basis.base_impedance, 
                        convert(Float64, df_cond.x_ohm_per_km[k]) / pu_basis.base_impedance,
                        convert(Float64, df_cond.max_i_ka[k]) / pu_basis.base_current,
                        convert(Float64, df_cond.cost_kdollars_per_km[k]) / pu_basis.base_money
                        ) 
            for k in 1:K]

    return conds
end


""" get_network_data 

    Arguments:
    ----------
        - NETWORK_PATH   : path of xlsx file containing the network data
        - voltage_limits : limits on the bus voltage in pu
        - PU_basis       : named tuple containing the per-unit basis used for the simulation
    
    Return values:
    -------------
        - A structure of type DistributionNetwork that contains the network data
        - A structure of type DistributionNetworkTopology that contains the topology of the network
"""
function get_network_data(  NETWORK_PATH::String; 
                            voltage_limits::VLIM=(V_min=0.95, V_max=1.05), 
                            pu_basis::PU_BASIS=define_pu_basis()
                            )
    
    nodes, subs_buses, load_buses = get_buses_data(NETWORK_PATH, voltage_limits, pu_basis)
    edges, lines = get_lines_data(NETWORK_PATH, nodes)
    conductors = get_conductors_data(NETWORK_PATH, pu_basis)

    return  DistributionNetwork(lines, subs_buses, load_buses, conductors, pu_basis), 
            DistributionNetworkTopology(nodes, edges)
end



