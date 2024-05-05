function sphere_sample(N, d)
    #SPHERE_SAMPLE Randomly sample N points from a unit d-dimensional sphere
    #d=3 sphere is a 2-sphere. S^(d-1)
    #
    #Input:
    #   N:  Number of points to sample
    #   d:  Dimension of sphere
    #
    #Output:
    #   X:  Points on sphere
    
    #dropped coordinate method
    U = randn!(zeros(N, d));
    normU = sqrt.(sum(U.^2, dims=2));
    X = U ./ normU;

    return X;
end

function ball_sample(N, d)
    #BALL_SAMPLE Randomly sample N points from a unit d-dimensional ball
    #d=3 sphere is a 2-sphere. S^(d-1)
    #
    #Input:
    #   N:  Number of points to sample
    #   d:  Dimension of ball
    #
    #Output:
    #   X:  Points on ball
    
    #dropped coordinate method
    u = sphere_sample(N, d+2);
    X = u[:, 1:d];
    
    return X;

end
    

function make_poly(model, vars, degree)
    #inspired by add_poly! from https://github.com/wangjie212/SparseDynamicSystem/
    mon = reverse(monomials(vars, 0:degree))
    coeff = @variable(model, [1:length(mon)])
    return coeff'*mon    
end

function poly_recover(v)
    #POLY_RECOVER after solving an optimization problem, recover the polynomials that define the solution
    # [multi-index]->real number (moment)

    mon = monomials(v)

    coe = coefficients(v)

    coe_value = value.(coe);

    v_sub = coe_value'*mon;


    return v_sub
end