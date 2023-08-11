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
    Structure representing an edge in the graph of the network

    Fields:
    -------
        - id        : The node i
        - from_node : The sending-end node of the edge
        - to_node   : The receiving-end node of the edge
""" 
mutable struct Edge 
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
    sending_lines::Dict{Int64, Vector{Int64}}
    receiving_lines::Dict{Int64, Vector{Int64}}
    function NetworkTopology(nodes::Vector{Node}, edges::Vector{Edge})

        N = length(nodes)
        sending_lines = Dict(n => Vector{Int64}() for n in 1:N)
        receiving_lines = Dict(n => Vector{Int64}() for n in 1:N)

        for e in edges 
            push!(sending_lines[e.from_node.id], e.id)
            push!(receiving_lines[e.to_node.id], e.id)
        end

        return new(nodes, edges, sending_lines, receiving_lines)
    end

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
    capa::Union{Nothing, Float64}
    P::Union{Nothing, Vector{Float64}} # in pu, always positive cause generation
    Q::Union{Nothing, Vector{Float64}} # in pu, can be positive or negative

    function PV(profile::Profile, PQ_diagram::PQ_DIAGRAM) 

        return new(profile, PQ_diagram, nothing, nothing, nothing) 
    end

end

mutable struct Storage 
    efficiency::Float64 
    capa::Union{Nothing, Float64}
    state::Union{Nothing, Vector{Float64}}
    P_charge::Union{Nothing, Vector{Float64}} # in pu, always positive cause generation
    P_discharge::Union{Nothing, Vector{Float64}} # in pu, always positive cause generation

    function Storage(efficiency::Float64) 
        return new(efficiency, nothing, nothing, nothing) 
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
    V_magn::Union{Nothing, Vector{Float64}}       #  in KV
    load_profile::Union{Nothing, Profile} 
    PV_installation::Union{Nothing, PV}   # 
    storage::Union{Nothing, Storage} 
    cos_phi::Float64                # cos(phi)

    function User(node::Node, V_limits::Union{Nothing, VLIM}, max_pv_capa::Float64, cos_phi::Float64) 
        return new(node, V_limits, max_pv_capa, nothing, nothing, nothing, nothing, cos_phi) 
    end
end

# Add to substation the fact that it can be built or not 
mutable struct Substation <: Bus
    node::Node
    built::Bool
    V_limits::Union{Nothing, VLIM}        # in kV
    S_rating_max::Float64                 # in MVA
    S_rating::Union{Nothing, Float64}     # in MVA
    V_magn::Union{Nothing, Vector{Float64}}       # in kV 
    P_sup::Union{Nothing, Vector{Float64}}        # in MVA 
    Q_sup::Union{Nothing, Vector{Float64}}        # in MVA
  
    
    function Substation(node::Node, V_limits::Union{Nothing, VLIM}, S_rating_max::Float64, S_rating::Float64)
        
        built = S_rating > 0 
       
        return new(node, built, V_limits, S_rating_max, S_rating, nothing, nothing, nothing)
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
    length::Float64   # in km
    built::Bool
    cost::Union{Nothing, Float64}
    conductor::Union{Nothing, Conductor}       # the conductor that was chosen
    I_magn::Union{Nothing, Vector{Float64}}            # in kA the magnitude of the current 
    P_send::Union{Nothing, Vector{Float64}}            # P flowing at sending end of the line
    Q_send::Union{Nothing, Vector{Float64}}            # Q flowing at sending end of the line
    P_rec::Union{Nothing, Vector{Float64}}            # P flowing at sending end of the line
    Q_rec::Union{Nothing, Vector{Float64}}            # Q flowing at sending end of the line

    function Line(edge::Edge, length::Float64) 
        if  length < 0 
            throw(DomainError("""[Line]: The length of a line cannot be negative."""))
        end
        
        return new(edge, length, false, nothing, nothing, nothing, nothing, nothing, nothing, nothing) 
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
    buses::Vector{Bus}
    #sub_buses::Vector{Substation} # contains the substation buses
    #load_buses::Vector{User} # contains the load buses
    conductors::Vector{Conductor}
    nb_substations::Int64 
    nb_init_subs::Int64 
    nb_loads::Int64 
    nb_lines::Int64
    nb_conductors::Int64
    pu_basis::PU_BASIS
end



#------------------------------------ 6. -------------------------------------#
#              Structures representing the DSO costs                          #
#-----------------------------------------------------------------------------#

struct DSOCosts
    substation::Float64     # [kEUR/MVA]
    loss::Float64           # loss cost [€/kWh]
    amortization::Int64     # amortization period of DSO investements [years]
    interest_rate::Float64    # interest rate DSO
    WEIGHT_I::Float64       # cost to violate 
    WEIGHT_V::Float64
    money_basis::Float64    # money basis in k€
    WEIGHT_OBJ::Float64
end

struct UserCosts
    PV::Float64             # PV capacity cost in [€/kWc]
    PV_conv::Float64        # PV converter cost in 
    storage::Float64
    EI::Float64             # energy imported [€/kWh]
    EE::Float64             # energy exported [€/kWh]
    DSOEI::Float64          # DSO energy imported cost [€/kWh]
    DSOEE::Float64          # DSO energy exported cost [€/kWh]
    GCC::Float64            # grid connection cost [€/kVA/year]
    amortization_PV::Int64  # amortization period of pv panels
    amortization_PVC::Int64 # amortization periof of pv converters
    amortization_storage::Int64
    money_basis::Float64    # money basis in k€
    WEIGHT_OBJ::Float64
end


#------------------------------------ 6. -------------------------------------#
#                       Structures representing a simulation                  #
#-----------------------------------------------------------------------------#
struct Simulation 
    network::Network
    network_topology::NetworkTopology
    DSO_costs::DSOCosts 
    User_costs::UserCosts
    nb_sign_days::Int64 # number of significative days to simulate
    nb_time_steps::Int64
    delta_t::Float64
    bilevel::Bool
    storage::Bool 
    formulation::Formulation
   

    function Simulation(network::Network, network_topology::NetworkTopology, DSO_costs::DSOCosts, User_costs::UserCosts, nb_sign_days::Int64, bilevel::Bool, storage::Bool,  
        formulation::Formulation)

        Ns            = network.nb_substations
        nb_time_steps = get_nb_time_steps(network.buses[Ns + 1].load_profile)
        delta_t       = network.buses[Ns + 1].load_profile.granularity

        println(delta_t)
        println(nb_time_steps)

        return new(network, network_topology, DSO_costs, User_costs, nb_sign_days, nb_time_steps, delta_t, bilevel, storage, formulation)
    end
end

#------------------------------------ 7. -------------------------------------#
#                         Additionnal functions                               #
#-----------------------------------------------------------------------------#

function get_nb_loads(n::Network)
    return n.nb_loads
end

function get_nb_substations(n::Network)
    return n.nb_substations
end

function get_nb_conductors(n::Network)
    return n.nb_conductors
end

function get_nb_buses(n::Network)
    return n.nb_substations + n.nb_loads
end

function get_nb_lines(n::Network)
    return n.nb_lines
end

function get_nb_nodes(t::NetworkTopology)
    return length(t.nodes)
end

function get_nb_time_steps(p::Profile)
    return length(p.time_serie)
end

StructTypes.StructType(::Type{Edge}) = StructTypes.Struct()
StructTypes.StructType(::Type{Node}) = StructTypes.Struct()
StructTypes.StructType(::Type{NetworkTopology}) = StructTypes.Mutable()
function save_struct(data::NetworkTopology, filename::String)
    f = open(filename,"w")
    JSON3.pretty(f, data)
    close(f)
end

StructTypes.StructType(::Type{Line}) = StructTypes.Mutable()
StructTypes.StructType(::Type{Network}) = StructTypes.Mutable()
function save_struct(data::Network, filename::String)
    f = open(filename,"w")
    JSON3.pretty(f, data)
    close(f)
end



