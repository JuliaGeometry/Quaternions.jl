module Quat
  importall Base
  
  include("Quaternion.jl")
  include("DualQuaternion.jl")

  export Quaternion
  export quat
  export DualQuaternion
  export dualquat
  export angleaxis
  export angle
  export axis
  export linpol
  export normalize
  export normalizea
  export dconj
end
