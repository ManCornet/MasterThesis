#----------------------------------------------------------------------------#
#                                                                            #
#                           - TFE : Bilevel DNEP -                           #
#                             University of Liege                            #
#                                                                            #
#----------------------------------------------------------------------------#
# Created By  : Manon Cornet
# Created Date: Wednesday May 24 2023
#
# structs:
#   file containing the data structures that contain the data required to run 
#   a simulation
#
# ============================================================================#
#                                   Imports                                   #
# ============================================================================#
using LaTeXStrings
using GraphRecipes, Graphs
using Plots
# abstract type Device end abstract type to represent PV or storage
#=
mutable struct Battery <: Device
    capa::Float64                  # in MVA
    SOC::Union{Nothing, Vector{Float64}}
    P::Union{Nothing, Vector{Float64}} # in MW, can be positive or negative
    Q::Union{Nothing, Vector{Float64}} # in MVar
end
=#
# ============================================================================#
#                       Structures for physical parameters                    #
# ============================================================================#

#------------------------------------ 1. -------------------------------------#
#      Structures representing the topology of the distribution network       #
#-----------------------------------------------------------------------------#

COORD = NamedTuple{(:x, :y), Tuple{Float64, Float64}}

""" Node: 

    Description:
    ------------
    Structure reprensenting a node in the graph of the network

    Fields:
    -------
        - id    : The node id
        - coord : The node coordinates
"""
struct Node 
    id::Int64 
    coord::COORD

    function Node(id::Int64, coord::COORD)
        if id < 0
            throw(DomainError("[Node]: The id of a node cannot be negative")) 
        end
        return new(id, coord)
    end
end

""" Edge: 

    Description:
    ------------
    Structure reprensenting an edge in the graph of the network

    Fields:
    -------
        - id        : The node i
        - from_node : The sending-end node of the edge
        - to_node   : The receiving-end node of the edge
""" 
struct Edge 
    id::Int64
    from_node::Node
    to_node::Node
end

""" NetworkTopology: 

    Description:
    ------------
    Structure reprensenting the topology of the network

    Fields:
    -------
        - nodes : The list of all nodes in the network
        - edges : The list of all the edges in the network
""" 
mutable struct NetworkTopology
    nodes::Vector{Node}
    edges::Vector{Edge} # contains the substation buses
end

#------------------------------------ 2. -------------------------------------#
#            Structures representing a time_serie (PV or load profiles)       #
#-----------------------------------------------------------------------------#
""" Profile: 

    Description:
    ------------
    Structure reprensenting a time serie

    Fields:
    -------
        - time_serie : The time serie
        - granularity : Its granularity [min.]
""" 
struct Profile 
    time_serie::Vector{Float64}    # Power profile [pu] 
    granularity::Int64             # [min] 
end

#------------------------------------ 3. -------------------------------------#
#    Structures representing devices that can be connected to a load node     #
#-----------------------------------------------------------------------------#
PQ_DIAGRAM = NamedTuple{(:max_q, :slope), Tuple{Float64, Float64}}

mutable struct PV 
    profile::Profile 
    PQ_diagram::PQ_DIAGRAM    # max Q consumed/produced wrt P_peak
    installed_capa::Union{Nothing, Vector{Float64}}
    P::Union{Nothing, Vector{Float64}} # in pu, always positive cause generation
    Q::Union{Nothing, Vector{Float64}} # in pu, can be positive or negative

    function PV(profile::Profile, PQ_diagram::PQ_DIAGRAM) 

        return new(profile, PQ_diagram, nothing, nothing, nothing) 
    end

end

#------------------------------------ 4. -------------------------------------#
#            Structures representing the electrical buses of the network      #
#-----------------------------------------------------------------------------#

VLIM = NamedTuple{(:V_min, :V_max), Tuple{Float64, Float64}}

abstract type Bus end

mutable struct User <: Bus
    node::Node
    V_limits::Union{Nothing, VLIM}        #  in KV
    max_pv_capa::Float64                  # in MVA
    V_magn::Union{Nothing, Float64}       #  in KV
    load_profile::Union{Nothing, Profile} 
    PV_installation::Union{Nothing, PV}   # 
    cos_phi::Float64                # cos(phi)

    function User(node::Node, V_limits::Union{Nothing, VLIM}, max_pv_capa::Float64,) 
        return new(node, V_limits, max_pv_capa::Float64,  nothing, nothing, nothing, 0.9) 
    end
end

mutable struct Substation <: Bus
    node::Node
    V_limits::Union{Nothing, VLIM}        # in kV
    S_rating_max::Float64                 # in MVA
    S_rating::Union{Nothing, Float64}     # in MVA
    V_magn::Union{Nothing, Float64}       # in kV 
    P_sup::Union{Nothing, Float64}        # in MVA 
    Q_sup::Union{Nothing, Float64}        # in MVA
    
    function Substation(node::Node, V_limits::Union{Nothing, VLIM}, S_rating_max::Float64) 

        return new(node, V_limits, S_rating_max, nothing, nothing, nothing, nothing)
    end
end

#------------------------------------ 5. -------------------------------------#
#            Structures representing the electrical lines of the network      #
#-----------------------------------------------------------------------------#

struct Conductor 
    name::String            
    r::Float64 # in Ohm/km
    x::Float64 # in Ohm/km
    max_i::Float64 # in kA
    cost::Float64 # in k€/km
    function Conductor( name::String, r::Float64, x::Float64, 
                        max_i::Float64, cost::Float64)
        if r < 0 || x < 0 || max_i < 0
            throw(DomainError("""[Conductor]: The resistance, reactance and max current of a 
                                              conductor cannot be negative."""))
        end

        return new(name, r, x, max_i, cost) 
    end
end

mutable struct Line
    edge::Edge
    line_length::Float64   # in km
    built::Bool
    line_cost::Union{Nothing, Float64}
    conductor::Union{Nothing, Conductor}       # the conductor that was chosen
    I_magn::Union{Nothing, Float64}            # in kA the magnitude of the current 
    P_send::Union{Nothing, Float64}            # P flowing at sending end of the line
    Q_send::Union{Nothing, Float64}            # Q flowing at sending end of the line

    function Line(edge::Edge, line_length::Float64) 
        if  line_length < 0 
            throw(DomainError("""[Line]: The length of a line cannot be negative."""))
        end
        
        return new(edge, line_length, false, nothing, nothing, nothing, nothing, nothing) 
    end
end


PU_BASIS = NamedTuple{(:base_power, :base_voltage, :base_current, :base_impedance), 
                        Tuple{Float64, Float64, Float64, Float64}} 


#------------------------------------ 5. -------------------------------------#
#              Structures representing the distribution network               #
#-----------------------------------------------------------------------------#
""" Network

    Fields:
    -------
        - nb_lines : The number of lines in the network
        - nb_buses : The number of buses in the network
        - nb_subs  : The number of substations in the network
        - lines    : The lines of the network 
        - nodes    : The nodes of the network
        - subs     : The substations of the network
"""
mutable struct Network
    lines::Vector{Line}
    sub_buses::Vector{Substation} # contains the substation buses
    load_buses::Vector{User} # contains the load buses
    conductors::Vector{Conductor}
    pu_basis::PU_BASIS
end



#------------------------------------ 6. -------------------------------------#
#              Structures representing the DSO costs                          #
#-----------------------------------------------------------------------------#

struct DSOCosts
    substation::Float64     # [kEUR/MVA]
    loss::Float64           # loss cost [€/kWh]
    amortization::Int64     # amortization period of DSO investements [years]
    interest_rate::Int64    # interest rate DSO
    money_basis::Float64    # money basis in k€
end

struct UserCosts
    PV::Float64             # PV capacity cost in [€/kWc]
    PV_conv::Float64        # PV converter cost in 
    EI::Float64             # energy imported [€/kWh]
    EE::Float64             # energy exported [€/kWh]
    DSOEI::Float64          # DSO energy imported cost [€/kWh]
    DSOEE::Float64          # DSO energy exported cost [€/kWh]
    GCC::Float64            # grid connection cost [€/kVA/year]
    amortization_PV::Int64  # amortization period of pv panels
    amortization_PVC::Int64 # amortization periof of pv converters
    money_basis::Float64    # money basis in k€
end


#------------------------------------ 6. -------------------------------------#
#                       Structures representing a simulation                  #
#-----------------------------------------------------------------------------#
struct Simulation 
    network::Network
    DSO_costs::DSOCosts 
    User_costs::UserCosts
    nb_time_steps::Int64
    delta_t::Float64

    function Simulation(network::Network, DSO_costs::DSOCosts, User_costs::UserCosts)
        nb_time_steps = get_nb_time_steps(network.load_buses[1].load_profile)
        delta_t       = network.load_buses[1].load_profile.granularity
        return new(network, DSO_costs, User_costs, nb_time_steps, delta_t)
    end
end

#------------------------------------ 7. -------------------------------------#
#                         Additionnal functions                               #
#-----------------------------------------------------------------------------#

function get_nb_load_bus(d::Network)
    return length(d.load_buses)
end

function get_load_buses(d::Network)
    return d.load_buses
end

function get_nb_sub_bus(d::Network)
    return length(d.sub_buses)
end

function get_nb_bus(d::Network)
    return length(d.load_buses) + length(d.sub_buses)
end

function get_nb_lines(d::Network)
    return length(d.lines)
end

function get_nb_conductors(d::Network)
    return length(d.conductors)
end

function get_nb_nodes(t::NetworkTopology)
    return length(t.nodes)
end

function get_nb_time_steps(p::Profile)
    return length(p.time_serie)
end

function get_granularity(p::Profile)
    return length(p.time_serie)
end

StructTypes.StructType(::Type{NetworkTopology}) = StructTypes.Mutable()
function save(R::NetworkTopology, filename::String)
    f = open(filename,"w")
    JSON3.write(f, R)
    close(f)
end



