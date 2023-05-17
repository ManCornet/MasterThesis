#-----------------------------------------------------------------------------
#
#                           - TFE : Bilevel DNEP - 
#                             University of Liege
#
#-----------------------------------------------------------------------------
# Created By  : Manon Cornet
# Created Date: Saturday March 18 2023
#
# test_results:
#   Test the results of the formulation
#
# =============================================================================
#                                   Imports
# =============================================================================


# =============================================================================
#                                  Functions
# =============================================================================
function check_rotated_cones(constraint, T, L, K)
    println("Rotated cone constraints that are not tight:")
    for t in T, l in L, k in K
        x = value(constraint[t, l, k])
        slack = 2 * x[1] * x[2] - (x[3]^2 + x[4]^2)
        if abs(slack) > 1e-3
            println(t, l, k)
        end
    end
end

function check_cones(constraint, T, Ns)
    println("Cone constraints that are not tight:")
    for t in T, n in Ns
        x = value(constraint[t, n])
        if t==1 && n==1
            println(x[1]^2)
            println(x[2]^2 + x[3]^2)
        end
        slack = x[1]^2 - (x[2]^2 + x[3]^2)
        if abs(slack) > 1e-3
            println(t)
            println(n)
        end
    end
end

