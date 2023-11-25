using Test
using Quaternions
using Aqua
using RealDot

# Avoid tests on unbound_args.
# https://github.com/JuliaGeometry/Quaternions.jl/pull/132#discussion_r1383817950
Aqua.test_all(Quaternions, unbound_args=false)

include("helpers.jl")
include("Quaternion.jl")
