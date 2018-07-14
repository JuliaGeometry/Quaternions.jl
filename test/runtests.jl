using Compat.Test
using Quaternions

test_associative(x, y, z, ⊗) = @test (x ⊗ y) ⊗ z ≈ x ⊗ (y ⊗ z)
test_commutative(x, y, ⊗) = @test x ⊗ y ≈ y ⊗ x
test_inverse(x, eins, ⊗, inv) = (@test x ⊗ inv(x) ≈ eins; @test inv(x) ⊗ x ≈ eins)
test_neutral(x, eins, ⊗) = (@test x ⊗ eins ≈ x; @test eins ⊗ x ≈ x)
test_monoid(x, y, z, ⊗, eins) = (test_associative(x, y, z, ⊗); test_neutral(x, eins, ⊗))
test_group(x, y, z, ⊗, eins, inv) = (test_monoid(x, y, z, ⊗, eins);test_inverse(x, eins, ⊗, inv))
test_multiplicative(x, y, ⊗, f) = @test f(x ⊗ y) ≈ f(x) ⊗ f(y)

include("test_Quaternion.jl")
