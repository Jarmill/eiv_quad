using eiv_quad
using Random
using JuMP
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
    

#perform sampling
Ubase = rand!(rng, zeros(m, T));
U = umax.*(2*Ubase .-1);

X = zeros(n, T);
X[:, 1] = rand!(rng, zeros(2, 1));
BU = B*U;
for t in 1:T-1
    X[:, t+1] = A*X[:, t] + BU[:, t];
end

noise_base = ball_sample(T, n)';
noise = R*noise_base;
X_noise = X + noise;

