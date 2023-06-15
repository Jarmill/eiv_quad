mutable struct psatz
    mu      #equality constrained multipliers
    s
    z       #polynomials 
    Gram    #Gram matrices for the multipliers
end

function make_mult_quad(vars)
    #make_mult_quad: create the multipliers associated with a single robust SOC constraint:
    #
    #[sum_i z_it, s_it; s_it, z_it] in SOS^2[A, B].
    #
    #Input:
    #   vars:   the variables A and B
    #
    #Output:    
    #   
    
end

function make_mult(data, vars)
    
end