__precompile__()

module Quaternions

  import Base: +, -, *, /, ^, ==
  import Base: abs, abs2, angle, conj, cos, exp, inv, isreal, isfinite, isinf, iszero, isnan, log, real, sin, sqrt
  import Base: convert, promote_rule, float
  import Base: rand, randn
  import LinearAlgebra: lyap, norm, normalize, sylvester
  using Random

  Base.@irrational INV_SQRT_EIGHT 0.3535533905932737622004 sqrt(big(0.125))

  include("Quaternion.jl")
  include("Octonion.jl")
  include("DualQuaternion.jl")

  export Quaternion,
    QuaternionF16,
    QuaternionF32,
    QuaternionF64,
    Octonion,
    OctonionF16,
    OctonionF32,
    OctonionF64,
    DualQuaternion,
    DualQuaternionF16,
    DualQuaternionF32,
    DualQuaternionF64
  export quat
  export octo
  export dualquat
  export angleaxis
  export angle
  export axis
  export linpol
  export normalize
  export normalizea
  export dconj
  export quatrand
  export nquatrand
  export octorand
  export dualquatrand
  export ndualquatrand
  export qrotation
  export rotationmatrix
  export slerp
end
