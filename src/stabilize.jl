mutable struct struct_data
    X_noise     #The noisy state data
    U           #The input data (noiseless for now)
    Sigma       #The covariance of the noise
    epslion     #The level of the noise
    tol         #The tolerance of enforcing superstability (default 1e-6)
end

mutable struct sys_vars #store the polynomial variables
    A
    B
end

function make_sys_vars(data)
    #make_sys_vars: generate the system variables of the plant parameters A and B
    n = size(X_noise, 1);
    m = size(U, 1);

    @polyvar A[1:n, 1:n]
    @polyvar B[1:n, 1:m] 
    
    return sys_vars(A, B)
end