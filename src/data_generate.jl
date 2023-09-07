# x Sigma^-1 x <= epsilon
mutable struct struct_data
    X           #The noisy state data
    U           #The noisy input data
    Sigma       #The covariance matrices of the noise
    epslion     #The levels of the noise    
    tol         #The tolerance of enforcing superstability (default 1e-6)
end

struct system
    A
    B    
end

function generate_data(sys, T, umax, epsilon, Sigma, rng)
    #generate a trajectory with only x-type noise
    #add u-type noise later
    n = size(sys.B, 1);
    m = size(sys.B, 2);
    Ubase = rand!(rng, zeros(m, T));
    U = umax.*(2*Ubase .-1);

    if typeof(epsilon) != Vector{Float64}
        #[eps_x, eps_u, eps_w]
        epsilon = epsilon*[1; 0; 0];
        Sigma = [Sigma, Sigma, Sigma];
    end

    #create the noise processes
    # s_Sigma = sqrt.(Sigma);

    x_noise_base = ball_sample(T, n)';
    u_noise_base = ball_sample(T, m)';
    w_noise_base = ball_sample(T, n)';

    x_noise = epsilon[1]*(sqrt(Sigma[1]) \ x_noise_base);
    u_noise = epsilon[2]*(sqrt(Sigma[2]) \ u_noise_base);
    w_noise = epsilon[3]*(sqrt(Sigma[3]) \ w_noise_base);

    X = zeros(n, T);
    X[:, 1] = rand!(rng, zeros(2, 1));
    BU = sys.B*U;
    w_noise_base = ball_sample(T, n)';
    for t in 1:T-1
        X[:, t+1] = sys.A*X[:, t] + BU[:, t] + w_noise[:, t];
    end

    
    # x_noise = epsilon[1]*(s_Sigma[1] \ x_noise_base);
    X_noise = X + x_noise;
    U_noise = U + u_noise;

    return struct_data(X_noise, U_noise, Sigma, epsilon, 1e-6);

end