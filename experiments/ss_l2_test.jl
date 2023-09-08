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
using Revise
using JuMP
using LinearAlgebra
rng = MersenneTwister(13);

# 2nd order system
A = [0.6863    0.3968
    0.3456    1.0388];
B = [0.4170    0.0001
    0.7203    0.3023];
n = 2;  m = 1;

umax = 1;           # input bound
T = 12;             # Time horizon
R = 0.075;           # radius for sampling (works for R=0.5)
    
model = Model();
epsilon = [R; R; 0]
sigma = [I, I, I];
sys = system(A, B);
data = generate_data(sys, T, umax, epsilon, sigma, rng);

#test out the psatz
vs = make_sys_vars(data);

order = 1;

ss_out_sparse = ss_quad(data, order, true);
# ss_out_dense = ss_quad(data, order, false);

#in this experiment, sparse succeeds and dense fails

Acl_sparse = sys.A + sys.B*ss_out_sparse.K;

eig_sparse = eigvals(Acl_sparse);
lam_sparse = maximum(sum(Acl_sparse, dims=2));