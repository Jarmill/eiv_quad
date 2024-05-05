#superstabilization of a system under elementwise-l2-bounded input and measurement noise
#
#
#a controller is synthesized with:
#   worst-case:     lambda = 0.9392765535
#   on true system: lambda = 0.4177393395
#
#ess gives an incorrect result for this system. Why?

using eiv_quad
using Random
using Revise
using JuMP
using LinearAlgebra
rng = MersenneTwister(13);

# 2nd order system
 

A = [0.626568  1.52102  
0.143801  1.94284];
B = [0.507894  -1.02176
-0.190421   1.47051];
n = 2;  m = 2;

umax = 1;           # input bound
T = 14;             # Time horizon
# Rx = 0.5;           # radius for sampling (works for R=0.5)
# Ru = 0.25;

Rx = 0.05;
Ru = 0.05;

model = Model();
epsilon = [Rx; Ru; 0]
sigma = [I, I, I];
sys = system(A, B);
data, data_true = generate_data(sys, T, umax, epsilon, sigma, rng);

#test out the psatz
vs = make_sys_vars(data);

order = 1;

# ss_out_sparse = ss_quad(data, order, true);
# ss_out_dense = ss_quad(data, order, false);

# ss_out_clean = ss_clean(sys);
# ss_out_dense = ss_quad(data, order, false);

ess_out_clean = ess_clean(sys);
ess_out_dense = ess_quad(data, order, false);

Acl_clean = sys.A + sys.B*ess_out_clean.K;
Acl_dense = sys.A + sys.B*ess_out_dense.K;

eig_clean = abs.(eigvals(Acl_clean))
eig_dense = abs.(eigvals(Acl_dense))

# lam_dense = maximum(sum(Acl_dense, dims=2));