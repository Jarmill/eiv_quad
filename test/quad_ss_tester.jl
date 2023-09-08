using eiv_quad
using Random
using Revise
using JuMP
using LinearAlgebra
rng = MersenneTwister(13);

# 2nd order system
A = [0.6863    0.3968
    0.3456    1.0388];
B = [0.4170    0.0001
    0.7203    0.3023];

#works with this (already stable) system (T=12, R=0.04)
# A = [0.6863    0.3968
#     0.3456    0.2];
# B = [0.4170
#     0.7203];

n = 2;  m = 1;

umax = 1;           # input bound
T = 12;             # Time horizon
# T = 20;
R = 0.05;           # radius for sampling
# R =
    
model = Model();
# epsilon = [R; 0; 0];
epsilon = [R; 0; 0];
sigma = [I, I, I];
sys = system(A, B);
data = generate_data(sys, T, umax, epsilon, sigma, rng);

#test out the psatz
# vs = sys_vars(A, B)
vs = make_sys_vars(data);

order = 1;

ss_out_sparse = ss_quad(data, order, true);
ss_out_dense = ss_quad(data, order, false);

#in this experiment, sparse succeeds and dense fails

Acl_sparse = sys.A + sys.B*ss_out_sparse.K;

eig_sparse = eigvals(Acl_sparse);

# q = tr(vs.A.^2) + 1;
# ps = quad_psatz(q, order, model, data, vs);
# ps_sparse = quad_psatz(q, order, model, data, vs, true);