#stabilize ()
struct output_ss
    status  #optimizing status
    K       #controller that performs superstabilization
    lambda  #superstabilizing gain bound
    blocksize #psd block sizes involved
end

function ss_clean(sys)
    #perform clean (no-noise) superstabilization 
    #
    #minimize the superstabilizing bound lambda

    #create the optimizer
    model = Model(optimizer_with_attributes(Mosek.Optimizer));

    #get the data properties and form variables
    n = size(data.X, 1);
    m = size(data.U, 1);

    #form the design variables
    lambda = JuMP.@variable(model);
    K = JuMP.@variable(model, [1:m, 1:n]);
    M = JuMP.@variable(model, [1:n, 1:n]);

    #create the constraints
    Acl = sys.A + sys.B*K;
    con_pos = vec(M - Acl);
    con_neg = vec(M + Acl);
    con_lam = vec(lambda .- sum(M, dims=2));

    con_all = [con_pos; con_neg; con_lam];

    @constraint(model, con_all >= 0);

    #impose the objective 

    @objective(model, Min, lambda)

    #solve the SOS problem
    optimize!(model)


    objv = objective_value(model);  
    status = termination_status(model)
    if status != MOI.OPTIMAL
        println("termination status: $status")
        status = primal_status(model)
        println("solution status: $status")
        lam_rec = Inf;
        K_rec = [];
    else
        K_rec = value(K);
        lam_rec = value(lambda);
    end    

    output = output_ss(lam_rec, K_rec, status)

    #recover the solution, process the output
    # output =  slice_recover(model, vars, poly, order, opts, info);

    return output
end

function ss_quad(data, order, SPARSE=false)
    #perform EIV (L2-noise) superstabilization 
    #
    #minimize the superstabilizing bound lambda

    #create the optimizer
    model = Model(optimizer_with_attributes(Mosek.Optimizer));

    #get the data properties and form variables
    n = size(data.X, 1);
    m = size(data.U, 1);
    T = size(data.U, 2);
    vs = make_sys_vars(data);
    vars_flat = vec([vs.A vs.B]);

    #form the design variables
    lambda = JuMP.@variable(model);
    K = JuMP.@variable(model, [1:m, 1:n]);

    M = Array{Polynomial}(undef, n, n);
    for i = 1:n
        for j = 1:n
            M[i, j] = make_poly(model, vars_flat, 2*order);
        end
    end

    #create the constraints

    Acl = vs.A + vs.B*K;
    con_pos = vec(M - Acl);
    con_neg = vec(M + Acl);
    con_lam = vec(lambda .- sum(M, dims=2));

    con_all = [con_pos; con_neg; con_lam];

    num_con = length(con_all);

    ps = Array{psatz}(undef, num_con, 1);
    blocksize = [];
    for k = 1:num_con
        ps[k] = quad_psatz(con_all[k], order, model, data, vs, SPARSE);
        blocksize = [blocksize; ps[k].blocksize];
    end
        

    #impose the objective 

    @objective(model, Min, lambda)

    #solve the SOS problem
    optimize!(model)


    objv = objective_value(model);  
    status = termination_status(model)
    if status != MOI.OPTIMAL
        println("termination status: $status")
        status = primal_status(model)
        println("solution status: $status")
        lam_rec = Inf;
        K_rec = [];
    else
        K_rec = value(K);
        lam_rec = value(lambda);
    end    

    output = output_ss(lam_rec, K_rec, status, blocksize);

    #recover the solution, process the output
    # output =  slice_recover(model, vars, poly, order, opts, info);

    return output
end
