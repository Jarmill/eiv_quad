using eiv_quad
using Random
using Revise
using JuMP
using LinearAlgebra
rng = MersenneTwister(135);

# 2nd order system
A = [0.6863    0.3968
    0.3456    1.0388];
B = [0.4170
    0.7203];

n = 2;  m = 1;

umax = 1;           # input bound
sig = 0.04^2;       # covariance of delta x
T = 12;             # Time horizon
R = 0.15;           # radius for sampling
    

# epsilon = [R; 0; 0];
epsilon = [R; 0; 0];
sigma = [I, I, I];
sys = system(A, B);
data = generate_data(sys, T, umax, epsilon, sigma, rng);

