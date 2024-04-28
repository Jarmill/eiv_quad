#try out the methods from https://arxiv.org/pdf/2402.04157.pdf
#paper by Andrea Bisoffi, Lidong Li, Claudio De Persis, Nima Monshizadeh

struct output_qmi
    status;#solver status
    K_rec; #Extracted controller
    P_rec; #Lyapunov matrix
end

function ref_theorem_1(data)
    
    n = size(data.X, 1)
    T = size(data.X, 2)
    m = size(data.U, 1)

    X0 = data.X[:, 1:(T-1)]
    U0 = data.U[:, 1:(T-1)]
    X1 = data.X[:, 2:T]

    #Theorem 1: energy bound
    scale = T*(2*data.epsilon[1] + data.epsilon[2])


    #QMI expression in (23)
    scA = [X0; U0] * [X0; U0]' - (scale*I)
    scB = -X1*[X0; U0]'
    scC = X1*X1' - (scale*I)
    #this fails assumption 1, scA is not a positive definite matrix


    status = 0;
    e_scA = eigvals(scA);

    z = zeros(n, n)


    if sum(e_scA <= 1e-5) >= 1
        #failure of Assumption 1
        return output_qmi(0, [], [])
    end

    pyt = [P Y'];
    block_control = [(-P-scC) z scB; z -P pyt; scB' py' -scA];



    @constraint(model, tau >= 0)
    @constraint(model, P - tol*I >= 0, PSDCone())
    @constraint(model, -block_control -tol*I>= 0, PSDCone())
    # @constraint(model, sum(diag(P)) == 1)

    @objective(model, Min, 0)

    #solve the SOS problem
    optimize!(model)

    objv = objective_value(model);  
    status = termination_status(model)
    if status != MOI.OPTIMAL
        println("termination status: $status")
        status = primal_status(model)
        println("solution status: $status")
        K_rec = [];
        # e_rec = [];
        P_rec = [];
    else
        P_rec = value(P)
        Y_rec = value.(Y)
        K_rec = Y_rec * inv(P_rec)

        tau_rec = value.(tau)

        # e_rec = abs.(eigvals(A+B*K_rec))
    end    
end

function ref_theorem_2(data)


    n = size(data.X, 1)
    T = size(data.X, 2)
    m = size(data.U, 1)

    X0 = data.X[:, 1:(T-1)]
    U0 = data.U[:, 1:(T-1)]
    X1 = data.X[:, 2:T]
    
    #Theorem 2: instantaneous noise bound (multi-block S-Lemma)

    #make the variables
    tol = 1e-4;
    model = Model(optimizer_with_attributes(Mosek.Optimizer));
    P = JuMP.@variable(model, [1:n, 1:n], Symmetric);
    Y = JuMP.@variable(model, [1:n, 1:n]);
    tau = JuMP.@variable(model, [1:(T-1)])
    alpha = JuMP.@variable(model, [1:1])
    #K = Y P^-1

    #make the error values
    err_term = [X1; -X0; -U0];

    # z = zeros(2, 2);
    zn = zeros(n, n);
    zm = zeros(m, m);
    znm = zeros(n, m)

    block_control = [(-P) zn znm zn; zn P Y' zn; znm' Y zm Y; zn zn Y' (-P)]

    #construct the noise bound
    theta = 2*data.epsilon[1]^2 + data.epsilon[2]^2;

    block_data = sum([tau[t]*(err_term[:, t]*(err_term[:, t]')-theta*I) for t = 1:T-1]);


    block_diag_embed = [block_data zeros(2*n+m, n); zeros(n, 2*n+m) zeros(n, n)]

    block_term = block_control - block_diag_embed

    @constraint(model, tau >= 0)
    @constraint(model, P - tol*I >= 0, PSDCone())
    @constraint(model, -block_term -tol*I>= 0, PSDCone())
    # @constraint(model, sum(diag(P)) == 1)

    @objective(model, Min, 0)

    #solve the SOS problem
    optimize!(model)

    objv = objective_value(model);  
    status = termination_status(model)
    if status != MOI.OPTIMAL
        println("termination status: $status")
        status = primal_status(model)
        println("solution status: $status")
        K_rec = [];
        # e_rec = [];
        P_rec = [];
    else
        P_rec = value(P)
        Y_rec = value.(Y)
        K_rec = Y_rec * inv(P_rec)

        tau_rec = value.(tau)

        # e_rec = abs.(eigvals(A+B*K_rec))
    end    


    return output_qmi(status, K_rec, P_rec);
end