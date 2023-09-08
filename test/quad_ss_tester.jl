using eiv_quad
using Random
using Revise
using JuMP
using LinearAlgebra
rng = MersenneTwister(13);

# 2nd order system
# A = [0.6863    0.3968
#     0.3456    1.0388];
A = [0.6863    0.3968
    0.3456    0.2];
B = [0.4170
    0.7203];

n = 2;  m = 1;

umax = 1;           # input bound
# T = 12;             # Time horizon
T = 20;
R = 0.01;           # radius for sampling
    
model = Model();
# epsilon = [R; 0; 0];
epsilon = [R; 0; 0];
sigma = [I, I, I];
sys = system(A, B);
data = generate_data(sys, T, umax, epsilon, sigma, rng);

#test out the psatz
# vs = sys_vars(A, B)
vs = make_sys_vars(data);

order = 2;

ss_out = ss_quad(data, order, false);
# q = tr(vs.A.^2) + 1;
# ps = quad_psatz(q, order, model, data, vs);
# ps_sparse = quad_psatz(q, order, model, data, vs, true);