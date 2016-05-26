using Base.Test
using Quaternions

check_associative(x, y, z, ⊗) = (x⊗y)⊗z ≈ x⊗(y⊗z)
check_commutative(x, y, ⊗) = x⊗y ≈ y⊗x
check_inverse(x, eins, ⊗, inv) = (x⊗inv(x) ≈ eins) & (inv(x)⊗x ≈ eins)
check_neutral(x, eins, ⊗) = (x⊗eins ≈ x) & (eins⊗x ≈ x)
check_monoid(x, y, z, ⊗, eins) = check_associative(x, y, z, ⊗) && check_neutral(x, eins, ⊗)
check_group(x, y, z, ⊗, eins, inv) = check_monoid(x,y,z, ⊗, eins) &&check_inverse(x, eins, ⊗, inv)
check_multiplicative(x,y,⊗, f) = f(x⊗y) ≈ f(x) ⊗ f(y)

include("test_Quaternion.jl")
