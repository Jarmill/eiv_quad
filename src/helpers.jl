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
    