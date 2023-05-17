module UpperLevel

import JuMP 
import XLSX
import DataFrames

# ------------------------------------------------------------------------------#
# ------------------------------------------------------------------------------#
#                                   1. DATA                                     #
# ------------------------------------------------------------------------------#
# ------------------------------------------------------------------------------#

# -- DEFINITION OF THE PER UNIT BASIS --
const BASE_VOLTAGE    = 34.5                             # [kV]
const BASE_POWER      = 1                                # [MVA]
const BASE_CURRENT    = BASE_POWER / BASE_VOLTAGE        # [kA]
const BASE_ADMITTANCE = BASE_CURRENT / BASE_VOLTAGE      # [S]
const BASE_IMPEDANCE  = 1/BASE_ADMITTANCE                # [Ohm]

# -- FUNCTION TO COMPUTE THE NPV --
function _PV_coeff(tau, lambda)
    return (1 - 1/(1 + tau)^lambda)/tau  
end

# -- FUNCTION TO COMPUTE THE CRF --
function _CRF(tau, n)
    # tau: interest rate 
    # n : number of annuity received
    return (tau * (1 + tau)^n)/((1 + tau)^n - 1)  
end

# -- FUNCTION TO PROCESS LINE CONDUCTORS --
function _process_conductors(df_cond::DataFrames.DataFrame, 
                            len_lines::Vector{Float64},  
                            nb_lines::Integer, 
                            )

    nb_cond = DataFrames.nrow(df_cond)
    max_i   = Array{Float64}(undef, nb_lines, nb_cond) # absolute, [pu]
    r       = Array{Float64}(undef, nb_lines, nb_cond) # absolute, [pu]
    x       = Array{Float64}(undef, nb_lines, nb_cond) # absolute, [pu]
    g       = Array{Float64}(undef, nb_lines, nb_cond) # absolute, [pu]
    b       = Array{Float64}(undef, nb_lines, nb_cond) # absolute, [pu]

    line_cost = Array{Float64}(undef, nb_lines, nb_cond) # [€/km]

    # Only take the first conductors of the list in the file
    for l in 1:nb_lines
        for k in 1:nb_cond
            max_i[l, k] = df_cond.max_i_ka[k] ./ BASE_CURRENT

            r[l, k] = len_lines[l] * df_cond.r_ohm_per_km[k] ./ BASE_IMPEDANCE
            x[l, k] = len_lines[l] * df_cond.x_ohm_per_km[k] ./ BASE_IMPEDANCE
            y = 1/(r[l, k]+im*x[l, k])

            g[l, k] = real(y) 
            b[l, k] = abs(imag(y))

            line_cost[l, k] = df_cond.cost_kdollars_per_km[k] * len_lines[l]
        end
    end
    return max_i, line_cost, r, x, g, b
end

function _read_data(filename)

    # -- FETCH THE DATA FROM THE EXCEL SHEET --
    df_bus  = DataFrames.DataFrame(XLSX.readtable(filename, "bus"))
    df_line = DataFrames.DataFrame(XLSX.readtable(filename, "line"))
    df_cond = DataFrames.DataFrame(XLSX.readtable(filename, "conductor"))

    # -- LINE PARAMETERS DEFINITION --
    L_size    = DataFrames.nrow(df_line)                                   # Number of lines in the network
    L         = 1:L_size                                                   # Line set
    line_ends = [(df_line.from_bus[l], df_line.to_bus[l]) for l in L]   # Indices of the line extremities
    len_lines = convert(Vector{Float64}, df_line.length_km)     

    # -- DEFINITION OF THE PHYSICAL QUANTITIES ASSOCIATED TO NETWORK LINES --
    K_size = DataFrames.nrow(df_cond)   # Number of conductor types
    K      = 1:K_size                   # Set of conductors

    max_i, line_cost, R, X, G, B, = _process_conductors(df_cond, len_lines, L_size)

    # -- BUS PARAMETERS DEFINITION --
    N_size = DataFrames.nrow(df_bus)            # Number of buses in the network
    N      = 1:N_size                           # Buses set

    Ns_size = sum(df_bus.type .== "substation") # Number of substation buses
    Nu_size = sum(df_bus.type .== "user")       # Number of load nodes

    Ns = 1:Ns_size                              # Set of substation buses
    Nu = (1:Nu_size) .+ Ns_size                 # Set of load buses

    # -- DEFINITION OF THE PHYSICAL QUANTITIES ASSOCIATED TO NETWORK BUSES --
    # Limits on voltage
    MIN_VOLTAGE = 0.97  # [pu]
    MAX_VOLTAGE = 1.03  # [pu]

    # Demand at buses
    # Assumption: load power factor is constant for all loads and is lagging (inductive)
    cos_phi = 0.9
    S_D     = convert(Vector{Float64}, df_bus.S_D_mva) ./ BASE_POWER
    P_D     = S_D * cos_phi
    Q_D     = S_D * sin(acos(cos_phi))

    # -- SUBSTATION PARAMETERS DEFINITION --
    S_rating_init    = convert(Vector{Float64}, df_bus.S_G_init_mva[Ns]) ./ BASE_POWER # [pu]
    S_rating_max     = convert(Vector{Float64}, df_bus.S_G_max_mva[Ns]) ./ BASE_POWER # [pu]
    sub_install_cost = 1e3      # k$
    sub_op_cost      = 0.1*1e-3 # k$/kVah^2

    # -- LINK BTW LINES AND NODES --
    Omega_sending   = Dict(n => [] for n in N)
    Omega_receiving = Dict(n => [] for n in N)
    for l in L
        push!(Omega_sending[line_ends[l][1]], l)
        push!(Omega_receiving[line_ends[l][2]], l)
    end

    # -- OBJECTIVE FUNCTION PARAMETERS --
    nb_years_planning = 1
    delta_t           = 1 # [h]

    tau = 0.1

    line_loss_factor = 0.35     # phi_l : loss factor of lines
    sub_loss_factor  = 0.35     # phi_s : cost per energy lost [€/kWh]

    K_l = _CRF(tau, 1)           # Capital recovery rate of line constructions
    K_s = _CRF(tau, 1)           # Capital recovery rate of substation construction or reinforcement

    loss_cost = 0.05*1e-3       # [k$/kWh]
    tau_l     = tau             # tau_l : interest rate of circuits
    tau_s     = tau             # tau_s : interest rate of substations

    f_l = _PV_coeff(tau_l, nb_years_planning)
    f_s = _PV_coeff(tau_s, nb_years_planning)

    # -- CREATION OF THE NETWORK DICT --

    network_dict = Dict(:bus       => (N, Omega_sending, Omega_receiving, MIN_VOLTAGE, MAX_VOLTAGE),
                        :load_bus  => (Nu, P_D, Q_D, delta_t),
                        :sub_bus   => (Ns, S_rating_init, S_rating_max),
                        :line      => (L, line_ends, max_i, R, X, G, B),
                        :conductor => (K)
                    )

    obj_fct_dict = Dict(:LF  => (line_loss_factor, sub_loss_factor),
                        :CRF => (K_l, K_s),
                        :NPV_coeff => (f_l, f_s),
                        :costs => (loss_cost, sub_install_cost, sub_op_cost, line_cost)
                        )

    return network_dict, obj_fct_dict
end



# ------------------------------------------------------------------------------#
# ------------------------------------------------------------------------------#
#                               2. CONFIGURATIONS                               #
# ------------------------------------------------------------------------------#
# ------------------------------------------------------------------------------#

# Function to add variables 
function _add_PowerFlow_variables( model::JuMP.Model, data::_KnapsackData, ::TimeIndependentConfig, ::BinaryKnapsackConfig,
    )
        return JuMP.@variable(model, x[keys(data.objects)], Bin)
    end

function _add_LoadOverSatisfaction_constraints( model::JuMP.Model, data::_KnapsackData, ::BinaryKnapsackConfig,
    )
        return JuMP.@variable(model, x[keys(data.objects)], Bin)
    end

# Idea here: functions that add constraints for the radiality when generation is added
function _add_radiality_constraints!(model::JuMP.Model, formulation=Formulation())::Nothing 

end

function _add_objective_function!()



function build_model(optimizer = nothing, formulation = Formulation())::JuMP.Model 

    if optimizer !== nothing
        set_optimizer(model, optimizer)
    end

    model = JuMP.Model(optimizer) 

    _add_knapsack_variables(model, data, config) 
    _add_PowerFlow_constraints(model, data, config) 

    _add_knapsack_objective(model, data, config) 

    JuMP.optimize!(model)

    if JuMP.termination_status(model) != JuMP.OPTIMAL
        @warn("Model not solved to optimality")
        return nothing 
    end
    return JuMP.value.(model[:x]) 
end




end