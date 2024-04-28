#try out the methods from https://arxiv.org/pdf/2402.04157.pdf on 
#the data used in ss_normal_test.jl
using JLD
using LinearAlgebra
using eiv_quad
using JuMP
using Mosek


#get the data from the ss_normal_test.jl experiment
stored = load("experiments/ss_normal.jld")
data = stored["data_chi"]

A = [0.6863    0.3968
    0.3456    1.0388];
B = [0.4170    0.0001
    0.7203    0.3023];

output = ref_theorem_2(data)