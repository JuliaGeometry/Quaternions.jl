module Quaternions

using Random
using LinearAlgebra
using RealDot

include("Quaternion.jl")

export Quaternion
export QuaternionF16, QuaternionF32, QuaternionF64
export quat
export imag_part
export slerp

end # module
