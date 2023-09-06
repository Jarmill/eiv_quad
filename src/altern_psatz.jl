mutable struct psatz
    mu      #equality constrained multipliers
    s
    z       #polynomials 
    Gram    #Gram matrices for the multipliers
end

mutable struct quad_mult
    s
    z
    Gram
end

mutable struct quad_mult_dense
    s
    z
    Gram
    Zeta
end

function make_mult_quad_dense(model, vars, order)
    #make_mult_quad_dense: create the multipliers associated with a single robust SOC constraint:
    #
    #[tau s\\ s I*tau] in SOS^(n+1) [A, B]. Each time index t has its own SOC constraint
    #    
    #
    #Input:
    #   model:  JuMP model which is containing all constraints
    #   vars:   the variables A and B
    #   order:  the order of polynomials to be used (degree/2)
    #
    #Output:    
    #   Gram:   Matrix that should be PSD
    #   s:      vector of polynomials
    #   z:      vector of polynomials

    #generate all monomials
    degree = 2*order;
    n = size(vars.A, 2)
    vars_flat = vec([vars.A vars.B]);
    # vars_c = col_vars(vars);
    mon = reverse(monomials(vars_flat, 0:order));

    #create the block Gram matrix
    len_gram = (n+1)*length(mon);
    nv = length(mon);
    
    #catch the polynomials and store them in arrays
    # Gram=Array{LinearAlgebra.Symmetric{VariableRef, Matrix{VariableRef}}}(undef, n, 1);
    # Zeta = Matrix{Polynomial{true, Float64}}}(undef)
    # s = Array{Polynomial}(undef, n, 1);
    # z = Array{Polynomial}(undef, n, 1);
    # zcorner = Array{Polynomial}(undef, n, 1);
    
    #the sum of the z elements is each entry of zcorner

        #declare the PSD matrix variable 
        Gram_curr = JuMP.@variable(model, [1:len_gram, 1:len_gram], PSD);
        Gram = Gram_curr;

        eye_n = 1* Matrix(I, n+1, n+1)
        mon_kron =  kron(mon, eye_n);        
        #index out the polynomial pieces
        # s[i] = mon'*Gram_s*mon;
        # z[i] = mon'*Gram_z*mon;
        Zeta = mon_kron'*Gram_curr*mon_kron;

        #fetch the multipliers
        s = Zeta[2:end, 1];
        z = Zeta[1, 1];

        nv = length(mon);

        @constraint(model, Gram_curr >= 0, PSDCone());
        for i = 1:n
            for j=i:n
                Zeta_curr = Zeta[i, j];                
                # zeta_coeff = coefficients(Zeta_curr);
                if i==j
                    #on-diagonal
                    @constraint(model, coefficients(Zeta_curr - z)==0);
                else
                    Zeta_opp = Zeta[j, i];
                    #off-diagonal
                    @constraint(model, coefficients(Zeta_opp + Zeta_curr)==0);
                end
            end
        end



        #add the constraints

    return quad_mult_dense(s, z, Gram, Zeta)


end

function make_mult_psd_dense(model, vars, order)
    #make_mult_quad_dense: create the multipliers associated with a single robust SOC constraint:
    #
    #[Zeta] in SOS^(n+1) [A, B]. Each time index t has its own SOC constraint
    #    
    #
    #Input:
    #   model:  JuMP model which is containing all constraints
    #   vars:   the variables A and B
    #   order:  the order of polynomials to be used (degree/2)
    #
    #Output:    
    #   Gram:   Matrix that should be PSD
    #   s:      vector of polynomials
    #   z:      vector of polynomials

    #generate all monomials
    degree = 2*order;
    n = size(vars.A, 2)
    vars_flat = vec([vars.A vars.B]);
    # vars_c = col_vars(vars);
    mon = reverse(monomials(vars_flat, 0:order));

    #create the block Gram matrix
    len_gram = (n+1)*length(mon);
    nv = length(mon);
    
    #catch the polynomials and store them in arrays
    # Gram=Array{LinearAlgebra.Symmetric{VariableRef, Matrix{VariableRef}}}(undef, n, 1);
    # Zeta = Matrix{Polynomial{true, Float64}}}(undef)
    # s = Array{Polynomial}(undef, n, 1);
    # z = Array{Polynomial}(undef, n, 1);
    # zcorner = Array{Polynomial}(undef, n, 1);
    
    #the sum of the z elements is each entry of zcorner

        #declare the PSD matrix variable 
        Gram_curr = JuMP.@variable(model, [1:len_gram, 1:len_gram], PSD);
        Gram = Gram_curr;

        eye_n = 1* Matrix(I, n+1, n+1)
        mon_kron =  kron(mon, eye_n);        
        #index out the polynomial pieces
        # s[i] = mon'*Gram_s*mon;
        # z[i] = mon'*Gram_z*mon;
        Zeta = mon_kron'*Gram_curr*mon_kron;

    return quad_mult_dense(Gram, Zeta)


end

# function make_mult_quad(model, vars, order)
#     #make_mult_quad: create the multipliers associated with a single robust SOC constraint:
#     #
#     #[sum_i z_it, s_it; s_it, z_it] in SOS^2[A, B]. Each time index t has its own SOC constraint
#     #    
#     #
#     #Input:
#     #   model:  JuMP model which is containing all constraints
#     #   vars:   the variables A and B
#     #   order:  the order of polynomials to be used (degree/2)
#     #
#     #Output:    
#     #   Gram:   Matrix that should be PSD
#     #   s:      vector of polynomials
#     #   z:      vector of polynomials

#     #generate all monomials
#     n = size(vars.A, 2)
#     var_flat = reshape([vars.A, vars.B], :, 1);
#     mon = reverse(monomials(vars, 0:degree));

#     #create the block Gram matrix
#     len_gram = 2*length(mon);
#     nv = length(mon);
    
#     #catch the polynomials and store them in arrays
#     Gram=Array{LinearAlgebra.Symmetric{VariableRef, Matrix{VariableRef}}}(undef, n, 1);
#     s = Array{Polynomial}(undef, n, 1);
#     z = Array{Polynomial}(undef, n, 1);
#     zcorner = Array{Polynomial}(undef, n, 1);
    
#     #the sum of the z elements is each entry of zcorner
#     zsum = 0;
#     for i = 1:n
#         #declare the PSD matrix variable 
#         Gram_curr = JuMP.@variable(model, [1:len_gram, 1:len_gram], PSD);
#         Gram[i] = Gram_curr;

#         #store the Gram matrix pieces
#         Gram_corner = Gram_curr[1:nv, 1:nv]
#         Gram_s = Gram_curr[1:nv, (nv+1):end];
#         Gram_z = Gram_curr[(nv+1):end, (nv+1):end];

#         #index out the polynomial pieces
#         s[i] = mon'*Gram_s*mon;
#         z[i] = mon'*Gram_z*mon;
#         zcorner[i] = mon'*Gram_corner'*mon;

#         zsum = zsum + z[i];
#     end

#     #now iterate back through the corners and set the corners equal to the sum
#     for i = 1:n
#         corner_diff = zsum - zcorner[i];
#         corner_coeff = coefficients(corner_diff, monomials(corner_diff));
#         @constraint!(model, corner_coeff==0);
#     end

#     return quad_mult(s, z, Gram)


# end

function make_mult(data, vars)
    
end