using Quaternions
using LinearAlgebra
using Random
using Test

struct MyReal <: Real
    val::Real
end
Base.:(/)(a::MyReal, b::Real) = a.val / b

function _quat(c::Complex{T}) where T
    Quaternion(reim(c)...,zero(T),zero(T))
end
function _quat(a::Real)
    Quaternion(a)
end

@testset "Quaternion" begin
    @testset "type aliases" begin
        @test QuaternionF16 === Quaternion{Float16}
        @test QuaternionF32 === Quaternion{Float32}
        @test QuaternionF64 === Quaternion{Float64}
    end

    @testset "Constructors" begin
        @testset "from coefficients" begin
            cs = [(1, 2.0, 3.0f0, 4//1), (1//1, 2.0f0, 3.0f0, 4)]
            @testset for coef in cs, T in (Float32, Float64, Int)
                q = @inferred Quaternion{T}(coef...)
                @test q isa Quaternion{T}
                @test q === Quaternion{T}(convert.(T, coef)...)
                q2 = @inferred Quaternion(convert.(T, coef)...)
                @test Quaternion(convert.(T, coef)...) === q
            end
        end
        @testset "from real" begin
            @testset for x in (-1//1, 1.0, 2.0), T in (Float32, Float64, Int, Rational{Int})
                coef = T.((x, 0, 0, 0))
                @test @inferred(Quaternion{T}(x)) === Quaternion{T}(coef...)
                @test @inferred(Quaternion(T(x))) === Quaternion{T}(coef...)
            end
        end
        @testset "from quaternion" begin
            @testset for q in (Quaternion(1, 2, 3, 4), QuaternionF64(0, 1, 0, 0)),
                T in (Float32, Float64)

                coef = T.((q.s, q.v1, q.v2, q.v3))
                @test @inferred(Quaternion{T}(q)) === Quaternion{T}(coef...)
                @test @inferred(Quaternion(q)) === q
            end
        end
    end

    @testset "==" begin
        @test Quaternion(1, 2, 3, 4) == Quaternion(1.0, 2.0, 3.0, 4.0)
        @test Quaternion(1, 2, 3, 4) != Quaternion(5, 2, 3, 4)
        @test Quaternion(1, 2, 3, 4) != Quaternion(1, 5, 3, 4)
        @test Quaternion(1, 2, 3, 4) != Quaternion(1, 2, 5, 4)
        @test Quaternion(1, 2, 3, 4) != Quaternion(1, 2, 3, 5)
    end

    @testset "isequal" begin
        @test isequal(Quaternion(1, 2, 3, 4), Quaternion(1.0, 2.0, 3.0, 4.0))
        @test !isequal(Quaternion(1, 2, 3, 4), Quaternion(5, 2, 3, 4))
        @test isequal(Quaternion(NaN, -0.0, Inf, -Inf), Quaternion(NaN, -0.0, Inf, -Inf))
        @test !isequal(Quaternion(NaN, 0.0, Inf, -Inf), Quaternion(NaN, -0.0, Inf, -Inf))
    end

    @testset "convert" begin
        @test convert(Quaternion{Float64}, 1) === Quaternion(1.0)
        @test convert(Quaternion{Float64}, Quaternion(1, 2, 3, 4)) ===
            Quaternion(1.0, 2.0, 3.0, 4.0)
        @test convert(Quaternion{Float64}, Quaternion(1.0, 2.0, 3.0, 4.0)) ===
            Quaternion(1.0, 2.0, 3.0, 4.0)
        @test convert(Quaternion{Float64}, Quaternion(0, 1, 0, 0)) ===
            Quaternion(0.0, 1.0, 0.0, 0.0)
    end

    @testset "promote" begin
        @test promote(Quaternion(1.0, 2, 3, 4), 1.0) ===
            (Quaternion(1.0, 2, 3, 4), Quaternion(1.0))
        @test promote(Quaternion(1.0f0, 2, 3, 4), 2.0) ===
            (Quaternion(1.0, 2, 3, 4), Quaternion(2.0))
        @test promote(Quaternion(1.0f0), Quaternion(2.0)) ===
            (Quaternion(1.0), Quaternion(2.0))

        @test Quaternion(1) == 1.0
    end

    @testset "shorthands" begin
        @test quat(1) === Quaternion(1)
        @test quat(1, 2, 3, 4) === Quaternion(1, 2, 3, 4)
        @test quat(Quaternion(1, 2, 3, 4)) === Quaternion(1, 2, 3, 4)
        @test quat([2, 3, 4]) == Quaternion{Int}[2, 3, 4]
        @test_throws ErrorException quat(Real[1,2,3])
        @test quat(Quaternion[1,2]) == Quaternion[1,2]

        @test quat(Int) === Quaternion{Int}
        @test quat(Float32) === Quaternion{Float32}
        @test quat(Quaternion{Int}) === Quaternion{Int}
        @test quat(Quaternion{Float32}) === Quaternion{Float32}

        # Note that `quat(1,missing,0,0)` throws an error.
        # This is the same behavior as `complex(1,missing)`.
        @test quat(missing) === missing
        @test quat(Missing) === Missing
        @test quat(Union{Missing, Int}) === Union{Missing, Quaternion{Int}}
    end

    @testset "random generation" begin
        @testset "quatrand" begin
            @test_deprecated quatrand()
            rng = Random.MersenneTwister(42)
            q1 = quatrand(rng)
            @test q1 isa Quaternion

            q2 = quatrand()
            @test q2 isa Quaternion
        end

        @testset "nquatrand" begin
            @test_deprecated nquatrand()
            rng = Random.MersenneTwister(42)
            q1 = nquatrand(rng)
            @test q1 isa Quaternion

            q2 = nquatrand()
            @test q2 isa Quaternion
        end

        @testset "rand($H)" for H in (QuaternionF32, QuaternionF64)
            rng = Random.MersenneTwister(42)
            q = rand(rng, H)
            @test q isa H

            qs = rand(rng, H, 1000)
            @test eltype(qs) === H
            @test length(qs) == 1000
            xs = map(qs) do q
                return [real(q); imag_part(q)...]
            end
            xs_mean = sum(xs) / length(xs)
            xs_var = sum(x -> abs2.(x .- xs_mean), xs) / (length(xs) - 1)
            @test all(isapprox.(xs_mean, 0.5; atol=0.1))
            @test all(isapprox.(xs_var, 1 / 12; atol=0.01))
        end

        @testset "randn($H)" for H in (QuaternionF32, QuaternionF64)
            rng = Random.MersenneTwister(42)
            q = randn(rng, H)
            @test q isa H

            qs = randn(rng, H, 10000)
            @test eltype(qs) === H
            @test length(qs) == 10000
            xs = map(qs) do q
                return [real(q); imag_part(q)...]
            end
            xs_mean = sum(xs) / length(xs)
            xs_var = sum(x -> abs2.(x .- xs_mean), xs) / (length(xs) - 1)
            @test all(isapprox.(xs_mean, 0; atol=0.1))
            @test all(isapprox.(xs_var, 1 / 4; atol=0.1))
        end
    end

    @testset "basic" begin
        q = randn(QuaternionF64)
        qnorm = sign(q)
        @test real(q) === q.s
        @test_throws MethodError imag(q)
        @test imag_part(q) === (q.v1, q.v2, q.v3)
        @test conj(q) === Quaternion(q.s, -q.v1, -q.v2, -q.v3)
        @test conj(qnorm) === Quaternion(qnorm.s, -qnorm.v1, -qnorm.v2, -qnorm.v3)
        @test conj(conj(q)) === q
        @test conj(conj(qnorm)) === qnorm
        @test float(Quaternion(1, 2, 3, 4)) === float(Quaternion(1.0, 2.0, 3.0, 4.0))
        @test Quaternions.abs_imag(q) ≈ abs(Quaternion(0, q.v1, q.v2, q.v3))
    end

    @testset "abs/abs_imag don't over/underflow" begin
        for x in [1e-300, 1e300, -1e-300, -1e300]
            @test abs(quat(x, 0, 0, 0)) == abs(x)
            @test abs(quat(0, x, 0, 0)) == abs(x)
            @test abs(quat(0, 0, x, 0)) == abs(x)
            @test abs(quat(0, 0, 0, x)) == abs(x)
            @test Quaternions.abs_imag(quat(0, x, 0, 0)) == abs(x)
            @test Quaternions.abs_imag(quat(0, 0, x, 0)) == abs(x)
            @test Quaternions.abs_imag(quat(0, 0, 0, x)) == abs(x)
        end
        @test isnan(abs(quat(NaN, NaN, NaN, NaN)))
        @test abs(quat(NaN, Inf, NaN, NaN)) == Inf
        @test abs(quat(-Inf, NaN, NaN, NaN)) == Inf
        @test abs(quat(0.0)) == 0.0
        @test abs(quat(Inf)) == Inf
        @test abs(quat(1, -Inf, 2, 3)) == Inf
        @test isnan(Quaternions.abs_imag(quat(0, NaN, NaN, NaN)))
        @test Quaternions.abs_imag(quat(0, Inf, NaN, NaN)) == Inf
        @test Quaternions.abs_imag(quat(0, NaN, -Inf, NaN)) == Inf
        @test Quaternions.abs_imag(quat(0.0)) == 0.0
        @test Quaternions.abs_imag(quat(0.0, 0.0, Inf, 0.0)) == Inf
        @test Quaternions.abs_imag(quat(0, 1, -Inf, 2)) == Inf
    end

    @testset "algebraic properties" begin
        for _ in 1:100, T in (Float32, Float64, Int32, Int64)
            if T <: Integer
                q, q1, q2, q3 = [Quaternion(rand((-T(100)):T(100), 4)...) for _ in 1:4]
            else
                q, q1, q2, q3 = randn(Quaternion{T}, 4)
            end

            # skewfield
            test_group(q1, q2, q3, +, zero(q), -)
            test_group(q1, q2, q3, *, one(q), inv)
            test_multiplicative(q1, q2, *, norm)
        end
    end

    @testset "inv does not under/overflow" begin
        x = 1e-300
        y = inv(x)
        @test isequal(inv(quat(x, 0.0, 0.0, 0.0)), quat(y, -0.0, -0.0, -0.0))
        @test isequal(inv(quat(0.0, x, 0.0, 0.0)), quat(0.0, -y, -0.0, -0.0))
        @test isequal(inv(quat(0.0, 0.0, x, 0.0)), quat(0.0, -0.0, -y, -0.0))
        @test isequal(inv(quat(0.0, 0.0, 0.0, x)), quat(0.0, -0.0, -0.0, -y))
        @test isequal(inv(quat(y, 0.0, 0.0, 0.0)), quat(x, -0.0, -0.0, -0.0))
        @test isequal(inv(quat(0.0, y, 0.0, 0.0)), quat(0.0, -x, -0.0, -0.0))
        @test isequal(inv(quat(0.0, 0.0, y, 0.0)), quat(0.0, -0.0, -x, -0.0))
        @test isequal(inv(quat(0.0, 0.0, 0.0, y)), quat(0.0, -0.0, -0.0, -x))
        @test isequal(inv(quat(-Inf, 1, -2, 3)), quat(-0.0, -0.0, 0.0, -0.0))
        @test isequal(inv(quat(1, -2, Inf, 3)), quat(0.0, 0.0, -0.0, -0.0))
    end

    @testset "isreal" begin
        @test isreal(Quaternion(1, 0, 0, 0))
        @test !isreal(Quaternion(2, 1, 0, 0))
        @test !isreal(Quaternion(2, 0, 1, 0))
        @test !isreal(Quaternion(2, 0, 0, 1))
    end

    @testset "iszero" begin
        @test iszero(Quaternion(0.0, 0.0, 0.0, 0.0))
        @test !iszero(Quaternion(1.0, 0.0, 0.0, 0.0))
        @test !iszero(Quaternion(0.0, 1.0, 0.0, 0.0))
        @test !iszero(Quaternion(0.0, 0.0, 1.0, 0.0))
        @test !iszero(Quaternion(0.0, 0.0, 0.0, 1.0))
    end

    @testset "isone" begin
        @test isone(Quaternion(1))
        @test !isone(Quaternion(-1))
        @test !isone(Quaternion(0, 1, 0, 0))
        @test !isone(Quaternion(1, 1, 0, 0))
        @test !isone(Quaternion(1, 0, 1, 0))
        @test !isone(Quaternion(1, 0, 0, 1))
    end

    @testset "isfinite" begin
        @test isfinite(Quaternion(1.0, 2.0, 3.0, 4.0))
        for value in (Inf, -Inf, NaN)
            @test !isfinite(Quaternion(value, 0.0, 0.0, 0.0))
            @test !isfinite(Quaternion(0.0, value, 0.0, 0.0))
            @test !isfinite(Quaternion(0.0, 0.0, value, 0.0))
            @test !isfinite(Quaternion(0.0, 0.0, 0.0, value))
            @test !isfinite(Quaternion(fill(value, 4)...))
        end
    end

    @testset "isinf" begin
        @test !isinf(Quaternion(1.0, 2.0, 3.0, 4.0))
        @test !isinf(Quaternion(1.0, 2.0, 3.0, NaN))
        for inf in (Inf, -Inf)
            @test isinf(Quaternion(inf, 0.0, 0.0, 0.0))
            @test isinf(Quaternion(0.0, inf, 0.0, 0.0))
            @test isinf(Quaternion(0.0, 0.0, inf, 0.0))
            @test isinf(Quaternion(0.0, 0.0, 0.0, inf))
        end
    end

    @testset "isnan" begin
        @test !isnan(Quaternion(1, 2, 3, 4))
        @test !isnan(Quaternion(1, 2, 3, Inf))
        @test !isnan(Quaternion(1, 2, 3, -Inf))
        @test isnan(Quaternion(NaN, 2, 3, 4))
        @test isnan(Quaternion(1, NaN, 3, 4))
        @test isnan(Quaternion(1, 2, NaN, 4))
        @test isnan(Quaternion(1, 2, 3, NaN))
    end

    @testset "isinteger" begin
        @test isinteger(quat(3))
        @test isinteger(quat(4.0))
        @test !isinteger(quat(4.1))
        @test !isinteger(quat(3, 1, 2, 3))
        @test !isinteger(quat(4, 0, 1, 0))
    end

    @testset "*" begin
        # verify basic correctness
        q1 = Quaternion(1,0,0,0)
        qi = Quaternion(0,1,0,0)
        qj = Quaternion(0,0,1,0)
        qk = Quaternion(0,0,0,1)
        @test q1 * q1 == q1
        @test q1 * qi == qi
        @test q1 * qj == qj
        @test q1 * qk == qk
        @test qi * q1 == qi
        @test qi * qi == -q1
        @test qi * qj == qk
        @test qi * qk == -qj
        @test qj * q1 == qj
        @test qj * qi == -qk
        @test qj * qj == -q1
        @test qj * qk == qi
        @test qk * q1 == qk
        @test qk * qi == qj
        @test qk * qj == -qi
        @test qk * qk == -q1
    end

    @testset "abs2 with $(T)" for T in (Float16, Float32, Float64)
        for _ in 1:100
            q = rand(Quaternion{T})
            @test abs2(q) == q'*q
        end
    end

    @testset "/" begin
        for _ in 1:100
            q, q2 = randn(QuaternionF64, 2)
            x = randn()
            @test q / q ≈ q \ q ≈ one(q)
            @test q / q2 ≈ q * inv(q2)
            @test q2 \ q ≈ inv(q2) * q
            @test q / x ≈ x \ q ≈ inv(x) * q
        end
        @testset "no overflow/underflow" begin
            @testset for x in [1e-300, 1e300, -1e-300, -1e300]
                @test quat(x) / quat(x) == quat(1)
                @test quat(x) / quat(0, x, 0, 0) == quat(0, -1, 0, 0)
                @test quat(x) / quat(0, 0, x, 0) == quat(0, 0, -1, 0)
                @test quat(x) / quat(0, 0, 0, x) == quat(0, 0, 0, -1)
                @test quat(0, x, 0, 0) / quat(x, 0, 0, 0) == quat(0, 1, 0, 0)
                @test quat(0, x, 0, 0) / quat(0, x, 0, 0) == quat(1, 0, 0, 0)
                @test quat(0, x, 0, 0) / quat(0, 0, x, 0) == quat(0, 0, 0, -1)
                @test quat(0, x, 0, 0) / quat(0, 0, 0, x) == quat(0, 0, 1, 0)
            end
            @testset for T in [Float32, Float64]
                o = one(T)
                z = zero(T)
                inf = T(Inf)
                nan = T(NaN)
                @testset for s in [1, -1], t in [1, -1]
                    @test isequal(quat(o) / quat(s*inf), quat(s*z, -z, -z, -z))
                    @test isequal(quat(o) / quat(s*inf, t*o, z, t*z), quat(s*z, -t*z, -z, -t*z))
                    @test isequal(quat(o) / quat(s*inf, t*nan, t*z, z), quat(s*z, nan, -t*z, -z))
                    @test isequal(quat(o) / quat(s*inf, t*inf, t*z, z), quat(s*z, -t*z, -t*z, -z))
                end
                @test isequal(quat(inf) / quat(inf, 1, 2, 3), quat(nan, nan, nan, nan))
                @test isequal(quat(inf) / quat(inf, 1, 2, -inf), quat(nan, nan, nan, nan))
            end
        end
    end

    @testset "^" begin
        @testset "^(::Quaternion, ::Real)" begin
            for _ in 1:100
                q = randn(QuaternionF64)
                @test q^2.0 ≈ q * q
                @test q^1.0 ≈ q
                @test q^-1.0 ≈ inv(q)
                @test q^1.3 ≈ exp(1.3 * log(q))
                @test q^7.8 ≈ exp(7.8 * log(q))
                @test q^1.3f0 ≈ exp(1.3f0 * log(q))
                @test q^7.8f0 ≈ exp(7.8f0 * log(q))
            end
        end
        @testset "^(::Quaternion, ::Quaternion)" begin
            @test Quaternion(ℯ, 0, 0, 0)^Quaternion(0, 0, π / 2, 0) ≈ Quaternion(0, 0, 1, 0)
            @test Quaternion(3.5, 0, 0, 2.3)^Quaternion(0.2, 0, 0, 1.7) ≈ Quaternion(
                real((3.5 + 2.3im)^(0.2 + 1.7im)), 0, 0, imag((3.5 + 2.3im)^(0.2 + 1.7im))
            )
            for _ in 1:100
                q, p = randn(QuaternionF64, 2)
                @test q^p ≈ exp(p * log(q))
            end
        end
    end

    @testset "non-analytic functions" begin
        q, q2 = randn(Quaternion{Float64}, 2)
        unary_funs = [conj, abs, abs2, norm, sign]
        # since every quaternion is conjugate to a complex number,
        # one can establish correctness as follows:
        @testset for fun in unary_funs
            for _ in 1:100
                c = randn(ComplexF64)
                q = _quat(c)
                @test @inferred(fun(q)) ≈ _quat(fun(c))
                @test q2 * fun(q) * inv(q2) ≈ fun(q2 * q * inv(q2))
            end
        end
    end

    @testset "extended complex analytic functions" begin
        # all complex analytic functions can be extended to the quaternions
        #! format: off
        unary_funs = [
            sqrt, inv, exp, exp2, exp10, expm1, log, log2, log10, log1p,
            sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh,
            csc, sec, cot, acsc, asec, acot, csch, sech, coth, acsch, asech, acoth,
            sinpi, cospi,
        ]
        #! format: on
        # since every quaternion is conjugate to a complex number,
        # one can establish correctness as follows:
        @testset for fun in unary_funs
            q, q2 = randn(QuaternionF64, 2)
            for _ in 1:100
                c = randn(ComplexF64)
                q = _quat(c)
                @test @inferred(fun(q)) ≈ _quat(fun(c))
                @test q2 * fun(q) * inv(q2) ≈ fun(q2 * q * inv(q2))
            end
        end

        @testset "identities" begin
            for _ in 1:100
                q = randn(QuaternionF64)
                @test inv(q) * q ≈ q * inv(q) ≈ one(q)
                @test sqrt(q) * sqrt(q) ≈ q
                @test exp(log(q)) ≈ q
                @test exp(zero(q)) === one(q)
                @test log(one(q)) === zero(q)
                @test exp2(log2(q)) ≈ q
                @test exp10(log10(q)) ≈ q
                @test expm1(log1p(q)) ≈ q
                @test sinpi(q) ≈ sin(π * q)
                @test cospi(q) ≈ cos(π * q)
                @test all(sincos(q) .≈ (sin(q), cos(q)))
                @test all(sincos(zero(q)) .≈ (sin(zero(q)), cos(zero(q))))
                if VERSION ≥ v"1.6"
                    @test all(sincospi(q) .≈ (sinpi(q), cospi(q)))
                    @test all(sincospi(zero(q)) .≈ (sinpi(zero(q)), cospi(zero(q))))
                end
                @test tan(q) ≈ cos(q) \ sin(q) ≈ sin(q) / cos(q)
                @test tanh(q) ≈ cosh(q) \ sinh(q) ≈ sinh(q) / cosh(q)
                @testset for (f, finv) in [
                    (sin, csc),
                    (cos, sec),
                    (tan, cot),
                    (sinh, csch),
                    (cosh, sech),
                    (tanh, coth),
                ]
                    @test f(q) ≈ inv(finv(q))
                end
                @testset for (f, finv) in [
                    (asin, acsc),
                    (acos, asec),
                    (atan, acot),
                    (asinh, acsch),
                    (acosh, asech),
                    (atanh, acoth),
                ]
                    @test f(q) ≈ finv(inv(q))
                end
            end
        end

        @testset "additional properties" begin
            @testset "log" begin
                @test log(zero(QuaternionF64)) === Quaternion(-Inf, 0, 0, 0)
                @test log(one(QuaternionF64)) === Quaternion(0.0, 0, 0, 0)
                @test log(-one(QuaternionF64)) ≈ Quaternion(0, π, 0, 0)
                x = rand()
                @test log(quat(x)) ≈ quat(log(x))
                @test log(quat(-x)) ≈ Quaternion(reim(log(complex(-x)))..., 0, 0)
            end

            @testset "exp" begin
                @test exp(Quaternion(0, 0, 0, 0)) === Quaternion(1.0, 0.0, 0.0, 0.0)
                @test exp(Quaternion(2, 0, 0, 0)) === Quaternion(exp(2), 0, 0, 0)
                @test exp(Quaternion(0, 2, 0, 0)) === Quaternion(cos(2), sin(2), 0, 0)
                @test exp(Quaternion(0, 0, 2, 0)) === Quaternion(cos(2), 0, sin(2), 0)
                @test exp(Quaternion(0, 0, 0, 2)) === Quaternion(cos(2), 0, 0, sin(2))

                @test norm(exp(Quaternion(0, 0, 0, 0))) ≈ 1
                @test norm(exp(Quaternion(2, 0, 0, 0))) ≠ 1
                @test norm(exp(Quaternion(0, 2, 0, 0))) ≈ 1
                @test norm(exp(Quaternion(0, 0, 2, 0))) ≈ 1
                @test norm(exp(Quaternion(0, 0, 0, 2))) ≈ 1

                @test exp(Quaternion(0.0, 0.0, 0.0, 0.0)) ===
                    Quaternion(1.0, 0.0, 0.0, 0.0)
                @test exp(Quaternion(2.0, 0.0, 0.0, 0.0)) ===
                    Quaternion(exp(2), 0, 0, 0)
                @test exp(Quaternion(0.0, 2.0, 0.0, 0.0)) ===
                    Quaternion(cos(2), sin(2), 0, 0)
                @test exp(Quaternion(0.0, 0.0, 2.0, 0.0)) ===
                    Quaternion(cos(2), 0, sin(2), 0)
                @test exp(Quaternion(0.0, 0.0, 0.0, 2.0)) ===
                    Quaternion(cos(2), 0, 0, sin(2))

                @test norm(exp(Quaternion(0.0, 0.0, 0.0, 0.0))) ≈ 1
                @test norm(exp(Quaternion(2.0, 0.0, 0.0, 0.0))) ≠ 1
                @test norm(exp(Quaternion(0.0, 2.0, 0.0, 0.0))) ≈ 1
                @test norm(exp(Quaternion(0.0, 0.0, 2.0, 0.0))) ≈ 1
                @test norm(exp(Quaternion(0.0, 0.0, 0.0, 2.0))) ≈ 1

                @test exp(Quaternion(0, 0, 0, 0)) isa Quaternion{Float64}
                @test exp(Quaternion(0.0, 0, 0, 0)) isa Quaternion{Float64}
                @test exp(Quaternion(0//1, 0, 0, 0)) isa Quaternion{Float64}
                @test exp(Quaternion(BigFloat(0), 0, 0, 0)) isa Quaternion{BigFloat}

                # https://github.com/JuliaGeometry/Quaternions.jl/issues/39
                @testset "exp(::Quaternion{Int})" begin
                    @test exp(Quaternion(1, 1, 1, 1)) ≈ exp(Quaternion(1.0, 1.0, 1.0, 1.0))
                end
            end
        end
    end

    @testset "sign" begin
        for _ in 1:100
            q = quatrand()
            qnorm = @inferred sign(q)
            @test abs(qnorm) ≈ 1
            @test q ≈ abs(q) * qnorm
            @test sign(qnorm) ≈ qnorm
        end
        @inferred(sign(Quaternion(1, 2, 3, 4)))
    end

    @testset "slerp" begin
        function qrotation(axis, theta)
            s, c = sincos(theta / 2)
            axis = normalize(axis)
            return Quaternion(c, s*axis[1], s*axis[2], s*axis[3])
        end
        @testset "q1=1" begin
            a = quat(1, 0, 0, 0.0)
            b = quat(0, 0, 0, 1.0)
            @test slerp(a, b, 0.0) ≈ a
            @test slerp(a, b, 1.0) ≈ b
            @test slerp(a, b, 0.5) ≈ qrotation([0, 0, 1], deg2rad(90))
            @test abs(slerp(a, b, 0.0)) ≈ 1
            @test abs(slerp(a, b, 1.0)) ≈ 1
            @test abs(slerp(a, b, 0.5)) ≈ 1
            @testset "scale $scale" for scale in (1, 1e-5, 1e-10)
                for _ in 1:100
                    q1 = quat(1, 0, 0, 0.0)
                    θ = rand() * π * scale
                    ax = randn(3)
                    q2 = qrotation(ax, θ)
                    qsmall = qrotation(ax, cbrt(eps()))
                    t = rand()
                    slerp(q1, q2, 0.0) ≈ q1
                    @test slerp(q1, q2, 0.0) ≈ q1
                    @test slerp(q1, q2, 1.0) ≈ q2
                    @test slerp(q1, q2, t) ≈ qrotation(ax, t * θ)
                    @test norm(slerp(q1, q2, t)) ≈ 1
                    @test slerp(q1, q2, 0.5) ≈ qrotation(ax, 0.5 * θ)
                    @test slerp(q1, q1, 0.5) ≈ q1
                    @test slerp(q1, qsmall, 0.5) ≈ sign((q1 + qsmall) / 2)
                end
            end
        end

        @testset "conjugate invariance" begin
            for _ in 1:100
                q, q1, q2 = randn(QuaternionF64, 3)
                ⊗(s, t) = s * t * inv(s)
                t = rand()
                @test q ⊗ slerp(q1, q2, t) ≈ slerp(q ⊗ q1, q ⊗ q2, t)
            end
        end

        @testset "type promotion" begin
            @test slerp(quat(1),quat(1),1) isa Quaternion{Float64}
            @test slerp(quat(1),quat(1),big(1)) isa Quaternion{BigFloat}
            @test slerp(quat(1),quat(1),Float32(1)) isa Quaternion{Float32}
            @test slerp(quat(1),quat(Float32(1)),Float32(1)) isa Quaternion{Float32}
            @test slerp(quat(Float64(1)),quat(Float32(1)),Float32(1)) isa Quaternion{Float64}
        end

        @testset "DomainError" begin
            @test_throws DomainError slerp(quat(1),quat(0),1)
            @test_throws DomainError slerp(quat(0),quat(1),0)
        end

        @testset "Normalizing input quaternions" begin
            for _ in 1:100
                q1 = randn(QuaternionF64)
                q2 = randn(QuaternionF64)
                t = rand()
                @test slerp(sign(q1),sign(q2),t) ≈ slerp(q1,q2,t)
            end
        end
    end

    @testset "sylvester/lyap" begin
        Ts = (Float64, QuaternionF64)
        Ttrips = [(Ta, Tb, Tc) for Ta in Ts for Tb in Ts for Tc in Ts]
        Ttrips = filter(x -> any(y -> y <: Quaternion, x), Ttrips)
        @testset "($Ta, $Tb, $Tc)" for (Ta, Tb, Tc) in Ttrips
            for _ in 1:100
                a = randn(Ta)
                b = randn(Tb)
                c = randn(Tc)
                x = @inferred sylvester(a, b, c)
                @test a * x + x * b ≈ -c
                x = @inferred sylvester(b, a, c)
                @test b * x + x * a ≈ -c
                @test iszero(sylvester(a, b, zero(c)))
                @test sylvester(a, zero(b), c) ≈ a \ -c
                @test sylvester(zero(a), b, c) ≈ -c / b
                @test iszero(sylvester(zero(a), b, zero(c)))
                @test iszero(sylvester(a, zero(b), zero(c)))
                @test iszero(sylvester(a, b, zero(c)))
                # @test isnan(sylvester(zero(a), zero(b), c))

                @test @inferred(lyap(a, c)) ≈ sylvester(a, a', c)
                @test @inferred(lyap(b, c)) ≈ sylvester(b, b', c)
                @test iszero(lyap(a, zero(c)))
            end
            @testset "nan/inf return same as for complex" begin
                Tza, Tzb, Tzc = map(
                    T -> T <: Quaternion ? complex(real(T)) : T, (Ta, Tb, Tc)
                )
                a, b = zero(Ta), zero(Tb)
                za, zb = zero(Tza), zero(Tzb)
                @testset for f in (one, zero, randn)
                    x = sylvester(a, b, f(Tc))
                    zx = sylvester(za, zb, f(Tzc))
                    if isinf(zx)
                        @test isinf(x)
                    elseif isnan(zx)
                        @test isnan(x)
                    end
                    if VERSION ≥ v"1.7"
                        x = lyap(a, f(Tc))
                        zx = lyap(za, f(Tzc))
                        if isinf(zx)
                            @test isinf(x)
                        elseif isnan(zx)
                            @test isnan(x)
                        end
                    end
                end
            end
        end
        @testset "rational" begin
            a = Quaternion(1, 2, 3, 4)
            b = Quaternion(1//2, 2//2, 3//2, 4//2)
            c = Quaternion(-1//2, 2//2, -4//2, -3//2)
            @test @inferred(sylvester(a, b, c)) isa Quaternion{Rational{Int}}
            @test @inferred(lyap(a, c)) isa Quaternion{Rational{Int}}
            a = randn(QuaternionF32)
            @test @inferred(sylvester(a, b, c)) isa QuaternionF32
            @test @inferred(lyap(a, c)) isa QuaternionF32
            null = zero(Quaternion{Rational{Int}})
            @test_throws DivideError sylvester(null, null, null)
            @test_throws DivideError lyap(null, null)
        end
    end

    @testset "RealDot with $T" for T in (Float32, Float64)
        for _ in 1:10
            q1 = randn(Quaternion{T})
            q2 = randn(Quaternion{T})
            # Check real∘dot is equal to realdot.
            @test real(dot(q1,q2)) == @inferred(realdot(q1,q2))
            # Check realdot is commutative.
            @test realdot(q1,q2) == realdot(q2,q1)
            # Check real∘dot is also commutative just in case.
            @test real(dot(q1,q2)) == real(dot(q2,q1))
            # Check the return type of realdot is correct.
            @test realdot(q1,q2) isa T
        end
    end

    @testset "widen" begin
        @test widen(Quaternion{Int}) === Quaternion{Int128}
        @test widen(QuaternionF32) === QuaternionF64
        @test widen(QuaternionF64) === Quaternion{BigFloat}
        @test widen(quat(1, 2, 3, 4)) === Quaternion{Int128}(1, 2, 3, 4)
        q = rand(QuaternionF32)
        @test widen(q) == convert(QuaternionF64, q)
        q = rand(QuaternionF64)
        @test widen(q) == convert(Quaternion{BigFloat}, q)
    end

    @testset "flipsign" begin
        q = rand(QuaternionF64)
        @test flipsign(q, 2) == q
        @test flipsign(q, -3) == -q
    end

    @testset "read/write" begin
        @testset "$T" for T in (Int16, Float32, Float64)
            io = IOBuffer(; read=true, write=true)
            q = rand(Quaternion{T})
            write(io, q)
            seek(io, 0)
            q2 = read(io, Quaternion{T})
            @test q == q2
        end
    end

    @testset "big" begin
        @test big(Quaternion{Int}) === Quaternion{BigInt}
        @test big(QuaternionF64) === Quaternion{BigFloat}
        @test big(quat(1, 2, 3, 4)) == Quaternion{BigInt}(1, 2, 3, 4)
        q = rand(QuaternionF64)
        @test big(q) == convert(Quaternion{BigFloat}, q)
    end

    @testset "round" begin
        q = quat(1.1, 2.5, -3.5, 2.3)
        @test round(q) == quat(1.0, 2.0, -4.0, 2.0)
        @test round(q; digits=1) == q
        @test round(q, RoundUp) == quat(2.0, 3.0, -3.0, 3.0)
        @test round(q, RoundUp; digits=1) == q
        @test round(q, RoundUp, RoundToZero) == quat(2.0, 2.0, -3.0, 2.0)
        @test round(q, RoundUp, RoundToZero; digits=1) == q
        rmodes = (RoundUp, RoundDown, RoundNearestTiesAway, RoundToZero)
        @test round(q, rmodes...) == quat(2.0, 2.0, -4.0, 2.0)
        @test round(q, rmodes...; digits=1) == q
    end
end
