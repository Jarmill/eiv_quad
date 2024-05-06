using eiv_quad
using Random
# using Revise
using JuMP
using LinearAlgebra
using StatsFuns
using JLD
rng = MersenneTwister(135);



n = 2;
m = 2;

M = 1;
Rx = 0.3;           # radius for sampling (works for R=0.5)
Ru = 0.15;           # radius for sampling (works for R=0.5)

umax = 1;           # input bound
T = 10;             # Time horizon

#standard deviations of noise process
epsilon = [Rx; Ru; 0]
sigma = [I, I, I];


N_experiments = 10;
# N_experiments = 3;
out_ss_sparse = Array{output_ss}(undef, N_experiments, 1);
out_ss_dense = Array{output_ss}(undef, N_experiments, 1);
# out_ss_full = Array{output_ss}(undef, N_experiments, 1);
out_ess_sparse = Array{output_ess}(undef, N_experiments, 1);
out_ess_dense = Array{output_ess}(undef, N_experiments, 1);
out_qmi = Array{output_qmi}(undef, N_experiments, 1);
system_test = Array{system}(undef, N_experiments, 1);
data_test = Array{struct_data}(undef, N_experiments, 1);
data_true_test = Array{struct_data}(undef, N_experiments, 1);


# order = 2;
order = 1;

for i = 1:N_experiments

    A = randn!(rng, zeros(n, n))*M;
    B = randn!(rng, zeros(n, m))*M;

    sys = system(A, B);
    data, data_true = generate_data(sys, T, umax, epsilon, sigma, rng, false);
    data_test[i] = data;
    data_true_test[i] = data_true;
    system_test[i] = sys;



    out_ss_dense[i] = ss_quad(data, order, false);
    out_ss_sparse[i] = ss_quad(data, order, true);
    out_ess_dense[i] = ess_quad(data, order, false);
    out_ess_sparse[i] = ess_quad(data, order, true);
    # K_rec = out_ss_dense[i].K;
    # Acl_rec = A + B*K_rec;

    # quad_out_1 = ref_theorem_1(data)
    out_qmi[i] = ref_theorem_2(data);

end


report_qmi = [out_qmi[i].status==MOI.OPTIMAL for i in 1:N_experiments]'
report_ss  = [(out_ss_dense[i].status==true ) && (out_ss_dense[i].lambda < 1) for i in 1:N_experiments]'
report_ss_sparse  = [(out_ss_sparse[i].status==true) && (out_ss_sparse[i].lambda < 1) for i in 1:N_experiments]'
report_ess  = [out_ess_dense[i].status==true for i in 1:N_experiments]'
report_ess_sparse  = [out_ess_sparse[i].status==true for i in 1:N_experiments]'

eig_ess_dense = [out_ess_dense[i].status ? abs.(eigvals(system_test[i].A + system_test[i].B*out_ess_dense[i].K)) : NaN  for i in 1:N_experiments]
eig_ss_dense = [out_ess_dense[i].status ? abs.(eigvals(system_test[i].A + system_test[i].B*out_ess_dense[i].K)) : NaN  for i in 1:N_experiments]
eig_ess_report = [(q!=NaN) && (sum(q.<1)==n) for q in eig_ess_dense]
eig_ss_report = [(q!=NaN) && (sum(q.<1)==n) for q in eig_ss_dense]
save("./experiments/multi_experiments/report_ess_mid_check.jld", "report_qmi", report_qmi, "report_ess", report_ess, "report_ess_sparse",
       report_ess_sparse, "report_ss_sparse", report_ss_sparse, "report_ss", report_ss,  "out_qmi", out_qmi, "out_ess_dense", out_ess_dense,
         "out_ss_dense", out_ss_dense, "out_ss_sparse", out_ss_sparse, "data_test", data_test, "system_test", system_test, "order", order, 
         "eig_ess_dense", eig_ess_dense, "eig_ess_report", eig_ess_report, "eig_ss_dense", eig_ss_dense, "eig_ss_report", eig_ss_report)