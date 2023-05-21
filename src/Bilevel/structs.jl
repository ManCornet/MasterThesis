

# abstract type Device end abstract type to represent PV or storage
#=
mutable struct Battery <: Device
    capa::Float64                  # in MVA
    SOC::Union{Nothing, Vector{Float64}}
    P::Union{Nothing, Vector{Float64}} # in MW, can be positive or negative
    Q::Union{Nothing, Vector{Float64}} # in MVar
end
=#

# All the structures required to put in argument of a model


#------------------------------------------------------------------------------------------
#                       Structure representing the topology of the network 
#------------------------------------------------------------------------------------------

COORD = NamedTuple{(:x, :y), Tuple{Float64, Float64}}

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

# has no direction
struct Edge 
    id::Int64
    from_node::Node
    to_node::Node
end

#------------------------------------------------------------------------------------------
#                               Structure for the load/PV profiles 
#------------------------------------------------------------------------------------------

struct Profile 
    nb_time_steps::Int64            # [-] Number of time steps
    granularity::Int64              # [min] 
    time_series::Vector{Float64}    # Apparent power profile [pu] 
    cos_phi::Float64                # cos(phi)
end

#------------------------------------------------------------------------------------------
#                               Structure for a PV installation 
#------------------------------------------------------------------------------------------

mutable struct PV 
    capa_max::Float64                  # in MVA
    profile::Profile 
    installed_capa::Float64
    P::Union{Nothing, Vector{Float64}} # in MW, always positive cause generation
    Q::Union{Nothing, Vector{Float64}} # in MVar, can be positive or negative

    function PV(capa_max::Float64, profile::Profile) 
        if capa_max < 0 
            throw(DomainError("[PV]: The max capa of a PV installation must be higher than 0")) 
        end
    
        return new(capa_max, profile) 
    end

end

#------------------------------------------------------------------------------------------
#                                   Bus structures 
#------------------------------------------------------------------------------------------

# We can have two types of buses: either a User or a substation
VLIM = NamedTuple{(:V_min, :V_max), Tuple{Float64, Float64}}

abstract type Bus end

mutable struct User <: Bus
    node::Node
    V_limits::Union{Nothing, VLIM}        #  in KV
    V_magn::Union{Nothing, Float64}       #  in KV
    load_profile::Union{Nothing, Profile} 
    PV_installation::Union{Nothing, PV}   # 
    
    function User(node::Node, V_limits::Union{Nothing, VLIM}) 
        return new(node, V_limits, nothing, nothing, nothing) 
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
    conductor::Union{Nothing, Conductor}    # the conductor that was chosen
    I_magn::Union{Nothing, Float64}                          # in kA the magnitude of the current 
    P_send::Union{Nothing, Float64}                        # P flowing at sending end of the line
    Q_send::Union{Nothing, Float64}                         # Q flowing at sending end of the line


    function Line(edge::Edge, line_length::Float64) 
        if  line_length < 0 
            throw(DomainError("""[Line]: The length of a line cannot be negative."""))
        end
        
        return new(edge, line_length, false, nothing, nothing, nothing, nothing, nothing) 
    end
end


PU_BASIS = NamedTuple{(:base_power, :base_voltage, :base_current, :base_impedance, :base_money), 
                        Tuple{Float64, Float64, Float64, Float64, Float64}} 

""" Distribution network structure

    Fields:
    -------
        - nb_lines : The number of lines in the network
        - nb_buses : The number of buses in the network
        - nb_subs  : The number of substations in the network
        - lines    : The lines of the network 
        - nodes    : The nodes of the network
        - subs     : The substations of the network
"""
mutable struct DistributionNetwork
    lines::Vector{Line}
    subs_buses::Vector{Substation} # contains the substation buses
    load_buses::Vector{User} # contains the load buses
    conductors::Vector{Conductor}
    pu_basis::PU_BASIS
   
end

mutable struct DistributionNetworkTopology
    nodes::Vector{Node}
    edges::Vector{Edge} # contains the substation buses
end


function get_nb_load_bus(d::DistributionNetwork)
    return length(d.load_buses)
end

function get_nb_subs_bus(d::DistributionNetwork)
    return length(d.subs_buses)
end

function get_nb_lines(d::DistributionNetwork)
    return length(d.lines)
end

function get_nb_conductors(d::DistributionNetwork)
    return length(d.conductors)
end


