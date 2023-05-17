
function _add_powerflow_eqs!(model::JuMP.Model, network_dict::Dict
                            )::Nothing


    P_ij = model[:P_ij]
    P_ji = model[:P_ij]
    Q_ij = model[:Q_ij]
    Q_ji = model[:Q_ji]
    G    = network_dict[:]
    for l in L, k in K 
        ifrom = network_dict[:line_ends][l][1]
        ito = network_dict[:line_ends][l][2]


        # CONSTRAINT (21)
        @constraint(model, P_ij[l, k] 
                           == G[l, k] * (X_i_ij[l, k, ifrom] 
                                                   - X_ij_re[l, k])
                           + B[l, k] * X_ij_im[l, k]
        )

        
        @constraint(model, P_ji[l, k] 
                           == G[l, k] * (X_i_ij[l, k, ito] 
                                                   - X_ij_re[l, k])
                            - B[l, k] * X_ij_im[l, k]) 
                            
      
        # CONSTRAINT (22)
        @constraint(model, Q_ij[l, k] 
                            == B[l, k] * (X_i_ij[l, k, ifrom] 
                                                     - X_ij_re[l, k])
                            - G[l, k] * X_ij_im[l, k])

        
        @constraint(model, Q_ji[l, k] 
                            ==  B[l, k] * (X_i_ij[l, k, ito] 
                                                     - X_ij_re[l, k])
                             + G[l, k] * (X_ij_im[l, k])) 
    
        
        # CONSTRAINT (23)
        @constraint(model, I_sqr[l, k] 
                            == (G[l, k]^2 + B[l, k]^2) 
                                * (X_i_ij[l, k, ifrom] + X_i_ij[l, k, ito] 
                                - 2 * X_ij_re[l, k]))
        
    end

end