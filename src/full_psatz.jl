#altern_psatz.jl
#certify that a polynomial q(A, B) is nonnegative over the L2-specified consistency set
#include all possible system and noise variables (A, B, dx, du), very expensive

function full_psatz(q, order, model, data, s_vars, n_vars)
    #certify that q>=0

    sys_vars_flat = vec([s_vars.A s_vars.B]);
    
    #handle the dx entries
    if data.epsilon[1]>0
        dx = n_vars.dx;
        dx_flat = vec(dx);
        dx_con = vec(data.epsilon[1]^2 .-sum(dx.^2, dims=1)); 
    else
        dx = n_vars.dx*0;
        dx_flat = [];
        dx_con = [];
    end

    #handle the du entries
    if data.epsilon[2]>0
        du = n_vars.du*0;
        du_flat = vec(n_vars.du);
        du_con = vec(data.epsilon[1]^2 .-sum(n_vars.dx.^2, dims=1));
    else
        du = n_vars.du*0;
        du_flat = [];
        du_con = [];
    end


    #perform the mu constraint

    n = size(data.X, 1);
    m = size(data.U, 1);
    T = size(data.U, 2);        

    vars_flat = vec([vars.A vars.B]);

    #create the residual
    Xnext = data.X[:, 2:end];
    Xprev = data.X[:, 1:end-1];
    U = data.U[:, 1:end-1];
    h = Xnext - vars.A*(Xprev+dx) - vars.B*(U+du);

    #process noise constraints
    if data.epsilon[3]>0
        w_con = vec(data.epsilon[3]^2 .-sum(h.^2, dims=1)); 
        eq_con = [];
    else
        w_con = [];
        eq_con = vec(h);
    end



    ineq_con = [dx_con; du_con; w_con];

    vars_flat_all = [sys_vars_flat; dx_flat; du_flat];

    model = add_psatz!(model, q, vars_flat_all, ineq_con, eq_con, order);







end