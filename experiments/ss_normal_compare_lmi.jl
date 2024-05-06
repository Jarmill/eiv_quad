#try out the methods from https://arxiv.org/pdf/2402.04157.pdf on 
#the data used in ss_normal_test.jl
using JLD
using LinearAlgebra
using eiv_quad
using JuMP
using Mosek


#get the data from the ss_normal_test.jl experiment
# stored = load("experiments/ss_normal.jld")
# data = stored["data_chi"]

A = [0.6863    0.3968
    -0.3456    1.0388];
B = [0.4170    0
    0.7203    -0.3023];

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

z = zeros(2, 2)

block_control = [(-P) z z z; z P Y' z; z Y z Y; z z Y' (-P)]

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
else
    P_rec = value(P)
    Y_rec = value.(Y)
    K_rec = Y_rec * inv(P_rec)

    tau_rec = value.(tau)

    e_rec = abs.(eigvals(A+B*K_rec))
end    
