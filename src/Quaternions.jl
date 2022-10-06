__precompile__()

module Quaternions

  import Base: +, -, *, /, ^, ==
  import Base: abs, abs2, angle, conj, cos, exp, inv, isreal, isfinite, isinf, iszero, isnan, log, real, sin, sqrt
  import Base: convert, promote_rule, float
  import Base: rand, randn
  import LinearAlgebra: lyap, norm, normalize, sylvester
  using LinearAlgebra: cross, dot
  using Random


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
  export normalize
  export normalizea
  export quatrand
  export nquatrand
  export qrotation
  export rotationmatrix
  export slerp
end
