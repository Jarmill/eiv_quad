#perform clean (no-noise) superstabilization 
#
#minimize the superstabilizing bound lambda

# 2nd order system
using eiv_quad
using JuMP
using LinearAlgebra
using Mosek
# A = [0.6863    0.3968
#     0.3456    1.0388];
A = [0.6863    0.3968
    0.3456    0.2];
B = [0.4170
    0.7203];

sys = system(A, B);

#create the optimizer
model = Model(optimizer_with_attributes(Mosek.Optimizer));

#get the data properties and form variables
n = size(B, 1);
m = size(B, 2);

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

# output = output_ss(lam_rec, K_rec, status)

#recover the solution, process the output
# output =  slice_recover(model, vars, poly, order, opts, info);

# return output