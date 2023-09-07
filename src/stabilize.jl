mutable struct sys_vars #store the polynomial variables
    A
    B
end

function make_sys_vars(data)
    #make_sys_vars: generate the system variables of the plant parameters A and B
    n = size(data.X, 1);
    m = size(data.U, 1);

    @polyvar A[1:n, 1:n]
    @polyvar B[1:n, 1:m] 
    
    return sys_vars(A, B)
end