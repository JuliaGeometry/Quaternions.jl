__precompile__()

module Quaternions

  import Base: +, -, *, /, ^, ==
  import Base: abs, abs2, angle, conj, cos, exp, inv, isreal, isfinite, isinf, iszero, isnan, log, real, sin, sqrt
  import Base: convert, promote_rule, float
  import Base: rand, randn
  import LinearAlgebra: lyap, norm, normalize, sylvester
  using Random
  using DualNumbers

  Base.@irrational INV_SQRT_EIGHT 0.3535533905932737622004 sqrt(big(0.125))

  """
      isunit(x)
  
  Return `true` if for the hypercomplex number `x`, `abs(x) == 1`.
  """
  function isunit end
  isunit(x::Real) = isone(abs(x))
  isunit(x::Complex) = isone(abs(x))
  isunit(x::DualNumbers.Dual) = isunit(DualNumbers.value(x)) & iszero(DualNumbers.epsilon(x))
  # TODO: replace this with something more efficient
  isunit(x::Number) = isone(abs(x))

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
  export isunit
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
