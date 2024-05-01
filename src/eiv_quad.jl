module eiv_quad

using JuMP
using DynamicPolynomials
using SparseDynamicSystem
using Random
using Revise
using LinearAlgebra
using Mosek

export ball_sample, sphere_sample, generate_data, system,
sys_vars, make_mult_quad, quad_mult, quad_psatz, make_sys_vars, output_ss, output_ess,
ss_clean, ss_quad, ess_clean, ess_quad, system, ss_quad_full, struct_data,
output_qmi, ref_theorem_1, ref_theorem_2 #comparison against reference https://arxiv.org/pdf/2402.04157.pdf

greet() = print("Hello World!")

include("make_vars.jl")
include("helpers.jl")
include("data_generate.jl")
include("stabilize.jl")
include("altern_psatz.jl")
include("full_psatz.jl")
include("reference_qmi.jl")
end # module eiv_quad
