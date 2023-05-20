#-----------------------------------------------------------------------------
#
#           - Planning & Operation of electric power and energy systems - 
#                   Homework : Implementation of the DNEP 
#
#                             University of Liege
#
#-----------------------------------------------------------------------------
# Created By  : Geoffrey Bailly, Manon Cornet
# Created Date: Tuesday November 21 2023
#
# main:
#   Main file
#
# =============================================================================
#                                   Imports
# =============================================================================

import XLSX
import DataFrames
using JuMP, Gurobi
using ArgParse, Printf


#include("planning_models.jl")
include("upper_level_models.jl")
include("read_xlsx_net.jl")

# =============================================================================
#                                   Functions
# =============================================================================

# Main argument parser
function parse_commandline()

    s = ArgParseSettings()

    @add_arg_table s begin
        "--formulation"
            help = "The formulation to use"
            default = 1
            arg_type = Int
    end
    return parse_args(s)
end

# Prints a basic section title in terminal
function section(title)

    # Number of letters to determine section size
    # title_size = length(title)
    boundary  = "="
    for i in 1:(30)
        boundary *= "="
    end
    
    # Printing section
    @printf("\n \n")
    @printf("%s", boundary)
    @printf(" %s ", title)
    @printf("%s", boundary)
    @printf("\n \n")
end

function print_table(table, table_name, unity)

    end_string = 12 - length(table_name)
    @printf("%s", table_name)
    for i in 1:end_string
        @printf(" ")
    end 
    @printf(": [ ")
    for (index, value) in enumerate(table)
        if typeof(value) == ComplexF64
            @printf("%.4f + i%.4f ", real(value), imag(value)) 

        elseif typeof(value) == Float64
            if value - floor(value) == 0
                @printf("%d ", Int(round(value))) 
            else
                @printf("%.4f ", value) 
            end

        else typeof(value) == Int64
            @printf("%d", value) 
        end
        if index < length(table)
            @printf(" ")
        end
    end
    @printf("]")
    @printf(" [%s]\n", unity) 
end

function print_2Dtable(table, table_name, unity)
    
    end_string = 12 - length(table_name)
    @printf("%s", table_name)
    for i in 1:end_string
        @printf(" ")
    end 
    @printf(": [ ")

    s = size(table)
    for i in 1:s[1]
        for j in 1:s[2]
            value = table[i, j]
            if typeof(value) == ComplexF64
                @printf("%.4f + i%.4f ", real(value), imag(value)) 
            elseif typeof(value) == Float64
                if value - floor(value) < 1e-16
                    @printf("%d ", Int(round(value)))
                else
                    @printf("%.10f ", value) 
                end
            else typeof(value) == Int64
                @printf("%d ", value) 
            end
        end
        if i < s[1]
            @printf("; ")
        end
    end
    @printf("]")
    @printf(" [%s]\n", unity)  
end

function check_rotated_cones()
    println("Rotated cone constraints that are not tight:")
    for k in 1:K, l in 1:K
        x = value(cone_28[k, l])
        slack = 2 * x[1] * x[2] - (x[3]^2 + x[4]^2)
        if abs(slack) > 1e-3
            println(k, l)
        end
    end
end

function check_cones()
    for n in 1:n_s
        println("Cone constraints that are not tight:")
        x = value(substation_capacity_limit[n])
        slack = x[1]^2 - (x[2]^2 + x[3]^2)
        if abs(slack) > 1e-3
            println(n)
        end
    end
end

# =============================================================================
#                                   Main
# =============================================================================

function main()

    # ================= Parsing arguments of main command line ================
    parsed_args = parse_commandline()
    args = [arg for  (arg, val) in parsed_args]
    vals = [val for (arg, val) in parsed_args]

    formulation = vals[1]
   
    println("Execution of the formulation $formulation")

    # ======================== MINLP formulation =======================
    if formulation == 1
        section("Formulation 1: MINLP")

        I_squared, V_re, V_im, x, alpha, beta,
        P_s, Q_s, P_conductor, Q_conductor, 
        obj, time = MINLP_formulation(  N, n_s, n_s_init, L, K ,
                                        K_l, K_s, 
                                        substation_utilization,
                                        S_rating_init, S_rating_max,
                                        conductance, susceptance, max_current, 
                                        line_cost, line_ends, line_length, line_loss, 
                                        interest_rate_losses, interest_rate_substation,
                                        substation_op_cost, substation_fixed_cost,cost_unit_loss,  
                                        Omega_sending, Omega_receiving, 
                                        P_demand, Q_demand
                                    )


        @printf("-------------------------------------------------------------\n")
        @printf("-------------------------- Results --------------------------\n")
        @printf("-------------------------------------------------------------\n")

        #I = [sqrt(i) for i in I_squared]
        V = [(v_re + im*v_im) for (v_re, v_im) in zip(V_re, V_im)]

        @printf("-> Number of buses: %d \n\n", N)
        @printf("-> Number of lines: %d \n\n", L)
        @printf("-> Number of conductor types: %d \n\n", K)
        @printf("-> Binary variables: \n\n")
        print_2Dtable(alpha, "Alpha" , "-") # 2D dimensions
        print_table(x, "x", "-")
        print_table(beta, "Beta", "-")
        @printf("\n-> Continuous variables: \n\n")
        print_2Dtable(I_squared, "I" , "pu^2")
        print_table(V, "V", "pu")
        print_table(P_s .* beta, "P_s", "kW")
        print_table(Q_s .* beta, "Q_s", "kVar")
        print_2Dtable(P_conductor .* alpha, "P_conductor", "kW")
        print_2Dtable(Q_conductor .* alpha, "Q_conductor", "kVar")
        
        

    # ======================== MISOCP formulation =======================
    elseif formulation == 2
        section("Formulation 2: MISOCP")

        I_squared, V_squared, x, alpha, beta, X_ij_re, 
        X_ij_im, X_i_ij, P_s, Q_s, P_cond_forward, Q_cond_forward,
        P_cond_backward, Q_cond_backward,
        obj, time, losses1, losses2 = MISOCP_formulation( N, n_s, n_s_init, L, K ,
                                        K_l, K_s, 
                                        substation_utilization,
                                        S_rating_init, S_rating_max,
                                        conductance, susceptance, max_current, 
                                        line_cost, line_ends, line_length, line_loss, 
                                        interest_rate_losses, interest_rate_substation,
                                        substation_op_cost, substation_fixed_cost, cost_unit_loss,  
                                        Omega_sending, Omega_receiving, 
                                        P_demand, Q_demand
                                        )    

        @printf("-------------------------------------------------------------\n")
        @printf("-------------------------- Results --------------------------\n")
        @printf("-------------------------------------------------------------\n")

        s = size(X_ij_re)
        #I = [sqrt(i) for i in I_squared]
        X = Array{Complex{Float64}}(undef, s[1], s[2])
        for i = 1:s[1], j = 1:s[2]
            X[i,j] = X_ij_re[i,j] + im*X_ij_im[i,j]
        end 
        #int_alpha = convert(Array{Bool,2}, alpha)
        #int_beta  = convert(Array{Bool,1}, beta)
        #int_x     = convert(Array{Bool,1}, x)

        @printf("-> Losses 1: %g \n\n", losses1)
        @printf("-> Losses 2: %g \n\n", losses2)
        println(Omega_receiving)
        println(Omega_sending)
        @printf("-> Number of buses: %d \n\n", N)
        @printf("-> Number of lines: %d \n\n", L)
        @printf("-> Number of conductor types: %d \n\n", K)
        @printf("-> Binary variables: \n\n")
        print_2Dtable(alpha, "Alpha" , "-") # 2D dimensions
        print_table(x, "x", "-")
        print_table(beta, "Beta", "-")
        @printf("\n-> Continuous variables: \n\n")
        print_2Dtable(I_squared * BASE_CURRENT, "I_squared" , "pu^2")
        print_table(V_squared, "V_squared", "pu^2")
        print_table(abs.(P_s .* beta), "P_s", "kW")
        print_table(abs.(Q_s .* beta), "Q_s", "kVar")
        print_2Dtable(abs.(P_cond_forward .* alpha), "P_cond_forward", "kW")
        print_2Dtable(abs.(P_cond_backward .* alpha), "P_cond_backward", "kW")
        print_2Dtable(abs.(Q_cond_forward .* alpha), "Q_cond_forward", "kVar")
        print_2Dtable(abs.(Q_cond_backward .* alpha), "Q_cond_backward", "kVar")
        print_2Dtable(X, "X", "pu^2")
        println(obj)
        println(I_squared)

        check_rotated_cones()
        check_cones() 

    end
end

main()