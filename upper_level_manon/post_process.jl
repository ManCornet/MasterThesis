#-----------------------------------------------------------------------------
#
#                           - TFE : Bilevel DNEP - 
#                             University of Liege
#
#-----------------------------------------------------------------------------
# Created By  : Manon Cornet
# Created Date: Wednesday March 29 2023
#
# post_process:
#   Contains the functions used to create an XLSX file with the results
#
# =============================================================================
#                                   Imports
# =============================================================================

import XLSX
import DataFrames
import Dates

# =============================================================================
#                                  Functions
# =============================================================================
function process2D_variable(var::Array{Float64, 2})
    return vec(sum(var, dims=2))
end

function process3D_variable(var::Array{Float64, 3})
    # 1st dimension: T, 2nd dimension: L, 3rd dimension: K
    T, L, ~ = size(var)
    processed_var = Array{Float64}(undef, T, L)
    # Trick to remove the kth dimension
    for t in 1:T
        processed_var[t, :] = sum(var[t, :, :], dims=2) 
    end    
   
    return processed_var
end

function process_X_i_ij_time(var::Array{Float64, 4})
    # 1st dim: T, 2nd dim: L, 3rd dim: K, 4th dim: N
    T, L, ~ , N = size(var)
    processed_var = Array{Float64}(undef, T, L, N)

    # Here we just want to keep the case k when it is chosen 
    # but does not mean that e
    for t in 1:T, i in 1:N   
        processed_var[t, :, i] = sum(var[t, :, :, i],  dims=2)
    end

    X_i_ij = Array{Float64}(undef, T, N)
    for t in 1:T, i in 1:N
        diff0 = processed_var[t, :, i][processed_var[t, :, i] .> 1e-7]
        if length(diff0) == 0
            X_i_ij[t, i] = 0
        else
            X_i_ij[t, i] = diff0[1]
        end
    end

    return X_i_ij
end

function process_X_i_ij(var::Array{Float64, 3})
    # 1st dim: T, 2nd dim: L, 3rd dim: K, 4th dim: N
    L, ~ , N = size(var)
    processed_var = Array{Float64}(undef, L, N)

    # Here we just want to keep the case k when it is chosen 
    # but does not mean that e
    for i in 1:N   
        processed_var[:, i] = sum(var[:, :, i],  dims=2)
    end

    X_i_ij = Vector{Float64}(undef, N)
    for i in 1:N
        diff0 = processed_var[:, i][processed_var[:, i] .> 1e-7]
        if length(diff0) == 0
            X_i_ij[i] = 0
        else
            X_i_ij[i] = diff0[1]
        end
    end

    return X_i_ij
end

# --- Function that adds a variable to an XLSX file ---

function add_var_to_XLSX(XLSX_PATH::String, 
                        var, var_name::String, 
                        dimensions::Vector{String}; 
                        date_range = nothing
                        )

    # Checking if the XLSX file exists
    nb_dim = ndims(var)
    isfile(XLSX_PATH) ? mode = "rw" : mode = "w"

    # Writing the variable in the XLSX file
    XLSX.openxlsx(XLSX_PATH, mode = mode) do xf
        # Mode = "w" only when the file must be created
        mode == "w" && XLSX.rename!(xf[1], var_name)
        mode == "rw" && !XLSX.hassheet(xf, var_name) && XLSX.addsheet!(xf, var_name)

        sheet = xf[var_name]

        if nb_dim  == 0
            column_names = [var_name]
            data = reshape([var], (1,:))
            #println("$var_name: $nb_dim")

        elseif nb_dim == 1 || nb_dim == 2
            # this only works if we make the assumption that the type can be vector, matrix or array
            #println("$var_name: $nb_dim")
            #println(var)
            nb_dim == 1 ? last_col = var_name : last_col=[dimensions[2] * "_$i" for i in 1:size(var)[2]]
        
            if dimensions[1] == "time" && !isnothing(date_range)
                dates = Dates.format.(Date.(date_range), "dd-mm-yyyy")
                times = Dates.format.(Time.(date_range), "HH:MM:SS")
                column_names = ["date"; "time"; last_col]
                data = hcat(dates, times, var)
            else
                column_names = [dimensions[1]; last_col] 
                data = hcat(1:size(var)[1], var)
            end
        end

        df = DataFrames.DataFrame(data, column_names)
        XLSX.writetable!(sheet, df)

    end
end