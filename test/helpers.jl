using Test

test_associative(x, y, z, ⊗) = @test (x ⊗ y) ⊗ z ≈ x ⊗ (y ⊗ z)

test_commutative(x, y, ⊗) = @test x ⊗ y ≈ y ⊗ x

function test_inverse(x, eins, ⊗, inv)
    @test x ⊗ inv(x) ≈ eins
    @test inv(x) ⊗ x ≈ eins
end

function test_neutral(x, eins, ⊗)
    @test x ⊗ eins ≈ x
    @test eins ⊗ x ≈ x
end

function test_monoid(x, y, z, ⊗, eins)
    test_associative(x, y, z, ⊗)
    test_neutral(x, eins, ⊗)
    return nothing
end

function test_group(x, y, z, ⊗, eins, inv)
    test_monoid(x, y, z, ⊗, eins)
    test_inverse(x, eins, ⊗, inv)
    return nothing
end

test_multiplicative(x, y, ⊗, f) = @test f(x ⊗ y) ≈ f(x) ⊗ f(y)
