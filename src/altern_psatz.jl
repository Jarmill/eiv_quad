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

function make_mult_quad(model, vars, order)
    #make_mult_quad: create the multipliers associated with a single robust SOC constraint:
    #
    #[sum_i z_it, s_it; s_it, z_it] in SOS^2[A, B]. Each time index t has its own SOC constraint
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
    n = size(vars.A, 2)
    var_flat = reshape([vars.A, vars.B], :, 1);
    mon = reverse(monomials(vars, 0:degree));

    #create the block Gram matrix
    len_gram = 2*length(mon);
    nv = length(mon);
    
    #catch the polynomials and store them in arrays
    Gram=Array{LinearAlgebra.Symmetric{VariableRef, Matrix{VariableRef}}}(undef, n, 1);
    s = Array{Polynomial}(undef, n, 1);
    z = Array{Polynomial}(undef, n, 1);
    zcorner = Array{Polynomial}(undef, n, 1);
    
    #the sum of the z elements is each entry of zcorner
    zsum = 0;
    for i = 1:n
        #declare the PSD matrix variable 
        Gram_curr = @variable(model, [1:len_gram, 1:len_gram], PSD);
        Gram[i] = Gram_curr;

        #store the Gram matrix pieces
        Gram_corner = Gram_curr[1:nv, 1:nv]
        Gram_s = Gram_curr[1:nv, (nv+1):end];
        Gram_z = Gram_curr[(nv+1):end, (nv+1):end];

        #index out the polynomial pieces
        s[i] = mon'*Gram_s*mon;
        z[i] = mon'*Gram_z*mon;
        zcorner[i] = mon'*Gram_corner'*mon;

        zsum = zsum + z[i];
    end

    #now iterate back through the corners and set the corners equal to the sum
    for i = 1:n
        corner_diff = zsum - zcorner[i];
        corner_coeff = coefficients(corner_diff, monomials(corner_diff));
        @constraint!(model, corner_coeff==0);
    end

    return quad_mult(s, z, Gram)


end

function make_mult(data, vars)
    
end