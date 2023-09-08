module eiv_quad

using JuMP
using DynamicPolynomials
using SparseDynamicSystem
using Random
using Revise
using LinearAlgebra
using Mosek

export ball_sample, sphere_sample, generate_data, system,
sys_vars, make_mult_quad, quad_mult, quad_psatz, make_sys_vars,
ss_clean, ss_quad, system

greet() = print("Hello World!")

include("make_vars.jl")
include("helpers.jl")
include("data_generate.jl")
include("stabilize.jl")
include("altern_psatz.jl")
include("full_psatz.jl")

end # module eiv_quad
