using eiv_quad
using JuMP
using LinearAlgebra
using Mosek
using Random
# A = [0.6863    0.3968
#     0.3456    1.0388];

rng = MersenneTwister(13);

A = [0.6863    0.3968
    0.3456    0.2];
B = [0.4170
    0.7203];
#data generation 
umax = 1;           # input bound
T = 14;             # Time horizon
Rx = 0.5;           # radius for sampling (works for R=0.5)
Ru = 0.5;           # radius for sampling (works for R=0.5)
    

epsilon = [Rx; Ru; 0]
sigma = [I, I, I];
sys = system(A, B);
data = generate_data(sys, T, umax, epsilon, sigma, rng, false);

sys = system(A, B);

out = ess_clean(sys);