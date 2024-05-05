#superstabilization of a system under elementwise-l2-bounded input and measurement noise
#
#
#the noise has a normal distribution
#
#safety probability 90%
#true lambda: 0.7329679715748868 (sparse)
#worst-case lambda: 0.997191322841913 (sparse)
#
#

# using BenchmarkTools 
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
n = size(B, 1);
m = size(B, 2);

umax = 1;           # input bound
T = 14;             # Time horizon
# Rx = 0.02;           # radius for sampling (works for R=0.02)
# Ru = 0.02;           # radius for sampling (works for R=0.02)

Rx = 0.04;           # radius for sampling (works for R=0.5)
Ru = 0.03;           # radius for sampling (works for R=0.5)


model = Model();
#standard deviations of noise process
epsilon = [Rx; Ru; 0]
sigma = [I, I, I];
sys = system(A, B);
data, data_true = generate_data(sys, T, umax, epsilon, sigma, rng, true);

#quantile for chi square, (sum of squares, so need to square Rx and Ru)
# safe scheme by Lemma 5 of https://arxiv.org/pdf/2211.05639.pdf
P_safe = 0.95;
P_error_unit = P_safe^(1/(2*T-1));
chi_quantile_x = Rx*sqrt(chisqinvcdf(n, P_error_unit));
chi_quantile_u = Ru*sqrt(chisqinvcdf(m, P_error_unit));
# chi_quantile_x = chisqinvcdf(n, 1-P_error)*Rx^2;
# chi_quantile_u = chisqinvcdf(m, 1-P_error)*Ru^2;
# prob_safe = (1-P_error)^(2*T);

epsilon_chi = [chi_quantile_x; chi_quantile_u; 0];
data_chi = data;
data_chi.epsilon = epsilon_chi;


#test out the psatz
vs = make_sys_vars(data);

order = 1;

# @btime ess_out_sparse = ss_quad(data_chi, order, true);
# @btime ess_out_sparse = ss_quad(data_chi, order+1, true);
# ess_out_dense = ess_quad(data_chi, order, false);
ss_out_dense = ss_quad(data_chi, order, false);
# @btime ess_out_full = ss_quad_full(data, order);
# K_rec = ss_out_dense.K;
K_rec = ss_out_dense.K;
Acl_rec = A + B*K_rec;
e_rec = abs.(eigvals(Acl_rec))
# @btime ss_out_full = ss_quad_full(data, order);

#in this experiment, sparse succeeds and dense fails
# if ss_out_sparse.status
#     Acl_sparse = sys.A + sys.B*ss_out_sparse.K;

#     eig_sparse = eigvals(Acl_sparse);
#     lam_sparse = maximum(sum(Acl_sparse, dims=2));
# end
