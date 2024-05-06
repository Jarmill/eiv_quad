#stabilize ()
struct output_ss
    status  #optimizing status
    K       #controller that performs superstabilization
    lambda  #superstabilizing gain bound
    blocksize #psd block sizes involved
end

struct output_ess
    status  #optimizing status
    K       #controller that performs extended superstabilization
    v       #Lyapunov vector
    S       #control-adjusted matrix
    blocksize #psd block sizes involved
    vs  #monomials for system matrices
    M   #recovered M function
end


function ss_clean(sys)
    #perform clean (no-noise) superstabilization 
    #
    #minimize the superstabilizing bound lambda

    #create the optimizer
    model = Model(optimizer_with_attributes(Mosek.Optimizer));

    #get the data properties and form variables
    n = size(sys.A, 1);
    m = size(sys.B, 1);

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
        K_rec = value.(K);
        lam_rec = value(lambda);
    end    

    output = output_ss(lam_rec, K_rec, status, [])

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

    @constraint(model, lambda >= -100);

        

    #impose the objective 

    @objective(model, Min, lambda)

    #solve the SOS problem
    optimize!(model)


    objv = objective_value(model);  
    status = termination_status(model)

    status_opt = (status==MOI.OPTIMAL) || (status == SLOW_PROGRESS);
    if !status_opt
        println("termination status: $status")
        status = primal_status(model)
        println("solution status: $status")
        lam_rec = Inf;
        K_rec = [];
    else
        K_rec = value.(K);
        lam_rec = value(lambda);
    end    

    output = output_ss(status_opt, K_rec, lam_rec, blocksize);

    #recover the solution, process the output
    # output =  slice_recover(model, vars, poly, order, opts, info);

    return output
end

function ss_quad_full(data, order)
    #perform EIV (L2-noise) superstabilization on the full psatz
    #
    #minimize the superstabilizing bound lambda

    #create the optimizer
    model = Model(optimizer_with_attributes(Mosek.Optimizer));

    #get the data properties and form variables
    n = size(data.X, 1);
    m = size(data.U, 1);
    vs = make_sys_vars(data);
    vn = make_noise_vars(data);
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

    ps = Array{full_mult}(undef, num_con, 1);
    blocksize = [];
    for k = 1:num_con
        ps[k] = full_psatz(con_all[k], order, model, data, vs, vn);
        # ps[k] = quad_psatz(con_all[k], order, model, data, vs, SPARSE);
        blocksize = [blocksize; ps[k].blocksize];
    end
        

    #impose the objective 
    @constraint(model, lambda >= 0);

    @objective(model, Min, lambda)

    #solve the SOS problem
    optimize!(model)


    objv = objective_value(model);  
    status = termination_status(model)

    status_opt = (status==MOI.OPTIMAL) || (status == SLOW_PROGRESS);
    if !status_opt
        println("termination status: $status")
        status = primal_status(model)
        println("solution status: $status")
        lam_rec = Inf;
        K_rec = [];
    else
        K_rec = value.(K);
        lam_rec = value(lambda);
    end    

    output = output_ss(status_opt, K_rec, lam_rec, blocksize);

    #recover the solution, process the output
    # output =  slice_recover(model, vars, poly, order, opts, info);

    return output
end

function ess_clean(sys)
    #perform clean (no-noise) extended superstabilization 
    #
    #feasibility problem. minimize the maximum value of v subject to sum(v)==1, v>=delta.

    #create the optimizer
    model = Model(optimizer_with_attributes(Mosek.Optimizer));

    #get the data properties and form variables
    n = size(sys.A, 1);
    m = size(sys.B, 2);

    #form the design variables
    lambda = JuMP.@variable(model);
    
    v = JuMP.@variable(model, [1:n]);
    S = JuMP.@variable(model, [1:m, 1:n]);
    M = JuMP.@variable(model, [1:n, 1:n]);

    X = diagm(vec(v));

    eta = 1e-3;

    #create the constraints
    Acl = sys.A*X + sys.B*S;
    con_pos = vec(M - Acl);
    con_neg = vec(M + Acl);
    
    #TODO: CHANGE THIS BACK TO ESS!
    # con_lam = vec(lambda .- sum(M, dims=2));
    # con_v  = vec(lambda .- v);
    
    #THE ESS STUFF:
    con_lam  = vec(lambda .- v);
    # con_v = vec(v .- sum(M, dims=2) .- eta);
    con_v = vec((1-eta)*v .- sum(M, dims=2));
    con_v_eta = vec(v .- eta);
    

    con_all = [con_pos; con_neg; con_lam; con_v; con_v_eta];

    @constraint(model, con_all >= 0);
    @constraint(model, sum(v)==n);
    

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
        S_rec = [];
        v_rec = [];
    else
        v_rec = value.(v);
        S_rec = value.(S);
        Xi_rec = diagm(1 ./ v_rec);
        K_rec = S_rec*Xi_rec;
        lam_rec = value(lambda);        
    end    

    output = output_ess(status, K_rec, v_rec, S_rec, status, [], value.(M))

    #recover the solution, process the output
    # output =  slice_recover(model, vars, poly, order, opts, info);

    return output
end

function ess_quad(data, order, SPARSE=false)
    #perform EIV (L2-noise) extended superstabilization 
    #
    #minimize the superstabilizing bound lambda

    #create the optimizer
    model = Model(optimizer_with_attributes(Mosek.Optimizer));

    #get the data properties and form variables
    n = size(data.X, 1);
    m = size(data.U, 1);
    vs = make_sys_vars(data);
    vars_flat = vec([vs.A vs.B]);

    #form the design variables
    lambda = JuMP.@variable(model);
    v = JuMP.@variable(model, [1:n]);
    S = JuMP.@variable(model, [1:m, 1:n]);

    X = diagm(vec(v));

    eta = 1e-2;

    M = Array{Polynomial}(undef, n, n);
    for i = 1:n
        for j = 1:n
            M[i, j] = make_poly(model, vars_flat, 2*order);
        end
    end

    #create the constraints

    
    Acl = vs.A*X + vs.B*S;
    

    #from SS
    # K = S;
    # Acl = vs.A + vs.B*K;
    # con_pos = vec(M - Acl);
    # con_neg = vec(M + Acl);
    # con_v = vec(lambda .- sum(M, dims=2));


    con_pos = vec(M - Acl);
    con_neg = vec(M + Acl);
    con_v = vec((1-eta)*v .- sum(M, dims=2));
    
    #TODO: CONVERT BACK TO ESS
    con_lam  = vec(lambda .- v);
    # con_lam  = vec(lambda);
    # con_lam  = lambda;
    con_lam_eta = vec(v .- eta);
    
    # con_v = vec(lambda .- sum(M, dims=2));
    # @constraint(model, v.==2*ones(2, 1));
    #BACK TO NORMAL

    con_all = [con_pos; con_neg; con_v];
 
    # finite-dimensional constraints
    @constraint(model, sum(v)==n);
    # @constraint(model, sum(v)==2);
    @constraint(model, con_lam >= 0);
    @constraint(model, lambda >= 0);
    @constraint(model, con_lam_eta >= 0);


    
    #TODO: For testing only
    # @constraint(model, v.==1);
    # @constraint(model, v.==0.5);

    #SOS constraints
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

    # status_opt = (status==MOI.OPTIMAL) || (status == SLOW_PROGRESS);
    status_opt = (status==MOI.OPTIMAL);
    if !status_opt
        println("termination status: $status")
        status = primal_status(model)
        println("solution status: $status")
        lam_rec = Inf;
        K_rec = [];
        S_rec = [];
        v_rec = [];
        vs_rec = [];
        M_rec = [];
    else
        v_rec = value.(v);
        S_rec = value.(S);
        Xi_rec = diagm(1 ./ v_rec);
        K_rec = S_rec*Xi_rec;
        lam_rec = value(lambda);
        vs_rec = vs;
        M_rec = M;
        for i=1:n
            for j=1:n
                M_rec[i, j] = value.(coefficients(M[i, j]))'*monomials(M[i, j])
            end
        end
    end    

    output = output_ess(status_opt, K_rec, v_rec, S_rec, status, vs_rec, M_rec)

    #recover the solution, process the output
    # output =  slice_recover(model, vars, poly, order, opts, info);

    return output
end