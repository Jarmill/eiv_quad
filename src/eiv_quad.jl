module eiv_quad

using JuMP
using DynamicPolynomials
using SparseDynamicSystem
using Random
using Revise
using LinearAlgebra

export ball_sample, sphere_sample, sys_vars, make_mult_quad,  make_mult_quad_dense, quad_mult, generate_data, system

greet() = print("Hello World!")

include("helpers.jl")
include("data_generate.jl")
include("stabilize.jl")
include("altern_psatz.jl")

end # module eiv_quad
