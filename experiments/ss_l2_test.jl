#superstabilization of a system under elementwise-l2-bounded input and measurement noise
#
#with elementwise L2 bounds of R = 0.075 for du and dx
#and a time horizon of T=12 samples on an open-loop-unstable system
#
#a controller is synthesized with:
#   worst-case:     lambda = 0.9392765535
#   on true system: lambda = 0.4177393395
#

using eiv_quad
using Random
# using Revise
using JuMP
using LinearAlgebra
using StatsFuns

rng = MersenneTwister(13);

# 2nd order system
A = [0.6863    0.3968
    0.3456    1.0388];
B = [0.4170    0.0001
    0.7203    0.3023];
n = 2;  m = 2;

umax = 1;           # input bound
# T = 15;             # Time horizon
# R = 0.1;           # radius for sampling (works for R=0.5)
T = 10;
R = 0.04;
    
model = Model();
epsilon = [R; R; 0]
sigma = [I, I, I];
sys = system(A, B);
data, data_true = generate_data(sys, T, umax, epsilon, sigma, rng, false);

#test out the psatz
vs = make_sys_vars(data);

order = 1;

# ss_out_sparse = ss_quad(data, order, true);
# ss_out_dense = ss_quad(data, order, false);

# ess_out_sparse = ess_quad(data, order, true);
ess_out_dense = ess_quad(data, order, false);

#in this experiment, sparse succeeds and dense fails

Acl_dense = sys.A + sys.B*ess_out_dense.K;

eig_sparse = eigvals(Acl_dense);
lam_sparse = maximum(sum(Acl_dense, dims=2));