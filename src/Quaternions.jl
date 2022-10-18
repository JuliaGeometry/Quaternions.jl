__precompile__()

module Quaternions

  using Random
  using LinearAlgebra

  include("Quaternion.jl")

  export Quaternion,
    QuaternionF16,
    QuaternionF32,
    QuaternionF64
  export quat
  export imag_part
  export angleaxis
  export angle
  export axis
  export quatrand
  export nquatrand
  export qrotation
  export rotationmatrix
  export sign_abs
  export slerp
end
