using eiv_quad
using Random
# using Revise
using JuMP
using LinearAlgebra
using StatsFuns
rng = MersenneTwister(135);



n = 5;
m = 3;

M = 1;
Rx = 0.04;           # radius for sampling (works for R=0.5)
Ru = 0.03;           # radius for sampling (works for R=0.5)

umax = 1;           # input bound
T = 20;             # Time horizon

#standard deviations of noise process
epsilon = [Rx; Ru; 0]
sigma = [I, I, I];


N_experiments = 100;
out_ss_dense = Array{output_ss}(undef, N_experiments, 1);
out_qmi = Array{output_qmi}(undef, N_experiments, 1);
system_test = Array{system}(undef, N_experiments, 1);
data_test = Array{struct_data}(undef, N_experiments, 1);


# order = 1;
order = 1;

for i = 1:N_experiments

    A = randn!(rng, zeros(n, n))*M;
    B = randn!(rng, zeros(n, m))*M;

    sys = system(A, B);
    data = generate_data(sys, T, umax, epsilon, sigma, rng, false);
    data_test[i] = data;
    system_test[i] = sys;



    out_ss_dense[i] = ss_quad(data, order, false);

    # K_rec = out_ss_dense[i].K;
    # Acl_rec = A + B*K_rec;
    # e_rec = abs.(eigvals(Acl_rec))

    # quad_out_1 = ref_theorem_1(data)
    out_qmi[i] = ref_theorem_2(data);

end


report_qmi = [out_qmi[i].status==MOI.OPTIMAL for i in 1:N_experiments]'
report_ss  = [out_ss_dense[i].status==true for i in 1:N_experiments]'