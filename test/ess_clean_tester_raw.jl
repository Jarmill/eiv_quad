#perform clean (no-noise) superstabilization 
#
#minimize the superstabilizing bound lambda

# 2nd order system
using eiv_quad
using JuMP
using LinearAlgebra
using Mosek

A = [1.3 -2.6 0.3 0
     -1  -1.5  2  1
     -0.5    -2.2    0.5  0.5
     1    0.5    0  1.3];

B = [1 0 1;
0 -1 1;
1 0 1;
0 1 1];  

# B = [1 0 0;
#      0 -1 0;
#      0 -1 0;
#      0 1 1];     

# A = [0.6863    0.3968
#     0.3456    1.0388];
# A = [0.6863    -0.3968
#     0.7    1.3];
# B = [0.4170
#     0.7203];

sys = system(A, B);

#create the optimizer
model = Model(optimizer_with_attributes(Mosek.Optimizer));

#get the data properties and form variables
n = size(B, 1);
m = size(B, 2);

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
con_lam  = vec(lambda .- v);
con_v = vec(v .- sum(M, dims=2) .- eta);

con_all = [con_pos; con_neg; con_lam; con_v];

@constraint(model, con_all >= 0);
@constraint(model, sum(v)==1);


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
    v_rec = value.(v);
    S_rec = value.(S);
    Xi_rec = diagm(1 ./ v_rec);
    K_rec = S_rec*Xi_rec;
    lam_rec = value(lambda);

    Acl = A + B*K_rec;
    ecl = eigvals(Acl);
end    

# output = output_ss(lam_rec, K_rec, status)

#recover the solution, process the output
# output =  slice_recover(model, vars, poly, order, opts, info);

# return output