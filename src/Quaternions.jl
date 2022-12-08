module Quaternions

using Random
using LinearAlgebra

include("Quaternion.jl")

export Quaternion
export QuaternionF16, QuaternionF32, QuaternionF64
export quat
export imag_part
export quatrand
export nquatrand
export slerp

end # module
