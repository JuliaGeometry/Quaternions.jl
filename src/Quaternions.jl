__precompile__()

module Quaternions

  using Compat

  import Base: +, -, *, /, ^
  import Base: abs, abs2, angle, conj, cos, exp, inv, isfinite, log, real, sin, sqrt
  import Base: convert, promote_rule
  import Compat.LinearAlgebra: norm, normalize

  include("Quaternion.jl")
  include("Octonion.jl")
  include("DualQuaternion.jl")

  export Quaternion
  export quat
  export Octonion
  export octo
  export DualQuaternion
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
