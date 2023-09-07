module eiv_quad

using JuMP
using DynamicPolynomials
using SparseDynamicSystem
using Random
using Revise
using LinearAlgebra

export ball_sample, sphere_sample, generate_data, system,
sys_vars, make_mult_quad, quad_mult, quad_psatz, make_sys_vars

greet() = print("Hello World!")

include("helpers.jl")
include("data_generate.jl")
include("stabilize.jl")
include("altern_psatz.jl")

end # module eiv_quad
