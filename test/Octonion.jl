using Quaternions
using LinearAlgebra
using Random
using Test

@testset "Octonion" begin
    @testset "type aliases" begin
        @test OctonionF16 === Octonion{Float16}
        @test OctonionF32 === Octonion{Float32}
        @test OctonionF64 === Octonion{Float64}
    end

    @testset "Constructors" begin
        @testset "from coefficients" begin
            cs = [(1, 2.0, 3.0f0, 4//1, 5, 6, 7, 8), (1//1, 2.0f0, 3.0f0, 4, 5, 6, 7, 8)]
            @testset for coef in cs, T in (Float32, Float64, Int), norm in (true, false)
                q = @inferred Octonion{T}(coef..., norm)
                @test q isa Octonion{T}
                @test q.norm === norm
                @test q === Octonion{T}(convert.(T, coef)..., norm)
                q2 = @inferred Octonion(convert.(T, coef)..., norm)
                @test Octonion(convert.(T, coef)..., norm) === q
                if !norm
                    @test Octonion(convert.(T, coef)...) === q
                end
            end
        end
        @testset "from real" begin
            @testset for x in (-1//1, 1.0, 2.0), T in (Float32, Float64, Int, Rational{Int})
                coef = T.((x, zeros(7)...))
                @test @inferred(Octonion{T}(x)) === Octonion{T}(coef..., isone(abs(x)))
                @test @inferred(Octonion(T(x))) === Octonion{T}(coef..., isone(abs(x)))
            end
        end
        @testset "from complex" begin
            @testset for z in (1 + 0im, -im, 1 + 2im),
                T in (Float32, Float64, Int, Rational{Int})

                coef = T.((reim(z)..., zeros(6)...))
                z2 = Complex{T}(z)
                norm = isone(abs(z))
                @test Octonion{T}(z) === Octonion{T}(coef..., norm)
                @test @inferred(Octonion(z2)) === Octonion{T}(coef..., norm)
            end
        end
        @testset "from quaternion" begin
            qs = (Quaternion(1, 2, 3, 4), QuaternionF64(0, 1, 0, 0, true))
            @testset for q in qs, T in (Float32, Float64)
                coef = T.((q.s, q.v1, q.v2, q.v3, zeros(4)...))
                q2 = Quaternion{T}(q)
                @test @inferred(Octonion{T}(q)) === Octonion{T}(coef..., q.norm)
                @test @inferred(Octonion(q2)) === Octonion{T}(coef..., q.norm)
            end
        end
        @testset "from octonion" begin
            os = (
                Octonion(1, 2, 3, 4, 5, 6, 7, 8), OctonionF64(0, 1, 0, 0, 0, 0, 0, 0, true)
            )
            @testset for o in os, T in (Float32, Float64)
                coef = T.((o.s, o.v1, o.v2, o.v3, o.v4, o.v5, o.v6, o.v7))
                @test @inferred(Octonion{T}(o)) === Octonion{T}(coef..., o.norm)
                @test @inferred(Octonion(o)) === o
            end
        end
        @testset "from vector" begin
            s = randn()
            v = randn(7)
            @test @inferred(Octonion(s, v)) === Octonion(s, v...)
            @test @inferred(Octonion(v)) === Octonion(0, v)
        end
    end

    @testset "==" begin
        @test Octonion(1.0, 2, 3, 4, 5, 6, 7, 8) == Octonion(1, 2, 3, 4, 5, 6, 7, 8)
        @test Octonion(1.0, 2, 3, 4, 5, 6, 7, 8) != Octonion(1, 2, 3, 4, 1, 2, 3, 4)
        @test Octonion(1, 0, 0, 0, 0, 0, 0, 0, false) ==
            Octonion(1, 0, 0, 0, 0, 0, 0, 0, true) # test that .norm field does not affect equality
    end

    @testset "convert" begin
        @test convert(Octonion{Float64}, 1) === Octonion(1.0)
        @test convert(Octonion{Float64}, Complex(1, 2)) ===
            Octonion(1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        @test convert(Octonion{Float64}, Quaternion(1, 2, 3, 4)) ===
            Octonion(1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0)
        @test convert(Octonion{Float64}, Quaternion(0, 1, 0, 0, true)) ===
            Octonion(0.0, 1.0, 0.0, 0.0, 0, 0, 0, 0, true)
        @test convert(Octonion{Float64}, Octonion(1, 2, 3, 4, 5, 6, 7, 8)) ===
            Octonion(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0)
        @test convert(Octonion{Float64}, Octonion(0, 0, 0, 0, 1, 0, 0, 0, true)) ===
            Octonion(0.0, 0.0, 0.0, 0.0, 1, 0, 0, 0, true)
    end

    @testset "promote" begin
        @test promote(Octonion(1.0, 2:8...), 1.0) === (Octonion(1.0, 2:8...), Octonion(1.0))
        @test promote(Octonion(1.0f0, 2:8...), 2.0) ===
            (Octonion(1.0, 2:8...), Octonion(2.0))
        @test promote(Octonion(1.0f0), 2 + 3im) ===
            (Octonion(1.0f0), Octonion(2.0f0 + 3.0f0im))
        @test promote(Octonion(1.0f0), Quaternion(1, 2, 3, 4)) ===
            (Octonion(1.0f0), Octonion(1.0f0:4.0f0..., fill(0, 4)...))
        @test promote(Octonion(1.0f0), Octonion(2.0)) === (Octonion(1.0), Octonion(2.0))

        @test Octonion(1) == 1.0
        @test Octonion(1, 2, fill(0, 6)...) == Complex(1.0, 2.0)
        @test Octonion(1) == 1.0
        @test Octonion(1:4..., fill(0, 4)...) == Quaternion(1.0:4.0...)
    end

    @testset "shorthands" begin
        @test octo(1) === Octonion(1) # checking the .norm field in particular
        @test octo(1, 0, 0, 0, 0, 0, 0, 0) === Octonion(1, 0, 0, 0, 0, 0, 0, 0) # checking the .norm field in particular
        @test octo(1, 2, 3, 4, 5, 6, 7, 8) === Octonion(1, 2, 3, 4, 5, 6, 7, 8)
        @test octo(Octonion(1, 0, 0, 0, 0, 0, 0, 0)) === Octonion(1, 0, 0, 0, 0, 0, 0, 0) # checking the .norm field in particular
        @test octo(Octonion(1, 2, 3, 4, 5, 6, 7, 8)) === Octonion(1, 2, 3, 4, 5, 6, 7, 8)
        @test octo(1, 0, 0, 0, 0, 0, 0, 0, false).norm == false # respect the .norm input (even if wrong)
        @test octo(1, 2, 3, 4, 5, 6, 7, 8, true).norm == true # respect the .norm input (even if wrong)
        @test octo(1, collect(2:8)) === Octonion(1:8...)
        @test octo(collect(2:8)) === Octonion(0, 2:8...)
    end

    @testset "random generation" begin
        @testset "octorand" begin
            o = octorand()
            @test o isa Octonion
            @test !o.norm
        end

        @testset "rand($H)" for H in (OctonionF32, OctonionF64)
            rng = Random.MersenneTwister(42)
            o = rand(rng, H)
            @test o isa H
            @test !o.norm

            os = rand(rng, H, 1000)
            @test eltype(os) === H
            @test length(os) == 1000
            xs = map(os) do o
                return [real(o); imag_part(o)...]
            end
            xs_mean = sum(xs) / length(xs)
            xs_var = sum(x -> abs2.(x .- xs_mean), xs) / (length(xs) - 1)
            @test all(isapprox.(xs_mean, 0.5; atol=0.1))
            @test all(isapprox.(xs_var, 1 / 12; atol=0.01))
        end

        @testset "randn($H)" for H in (OctonionF32, OctonionF64)
            rng = Random.MersenneTwister(42)
            o = randn(rng, H)
            @test o isa H
            @test !o.norm

            os = randn(rng, H, 10000)
            @test eltype(os) === H
            @test length(os) == 10000
            xs = map(os) do o
                return [real(o); imag_part(o)...]
            end
            xs_mean = sum(xs) / length(xs)
            xs_var = sum(x -> abs2.(x .- xs_mean), xs) / (length(xs) - 1)
            @test all(isapprox.(xs_mean, 0; atol=0.1))
            @test all(isapprox.(xs_var, 1 / 8; atol=0.1))
        end
    end

    @testset "basic" begin
        q = randn(OctonionF64)
        qnorm = normalize(q)
        @test real(q) === q.s
        @test_throws MethodError imag(q)
        @test @test_deprecated(Quaternions.imag(q)) == [q.v1, q.v2, q.v3, q.v4, q.v5, q.v6, q.v7]
        @test imag_part(q) === (q.v1, q.v2, q.v3, q.v4, q.v5, q.v6, q.v7)
        @test conj(q) ===
            Octonion(q.s, -q.v1, -q.v2, -q.v3, -q.v4, -q.v5, -q.v6, -q.v7, q.norm)
        @test conj(qnorm) === Octonion(
            qnorm.s,
            -qnorm.v1,
            -qnorm.v2,
            -qnorm.v3,
            -qnorm.v4,
            -qnorm.v5,
            -qnorm.v6,
            -qnorm.v7,
            qnorm.norm,
        )
        @test conj(conj(q)) === q
        @test conj(conj(qnorm)) === qnorm
        @test float(Octonion(1:8...)) === Octonion(1.0:8.0...)
        @test Quaternions.abs_imag(q) ==
            abs(Octonion(0, q.v1, q.v2, q.v3, q.v4, q.v5, q.v6, q.v7))
    end

    @testset "algebraic properties" begin
        for _ in 1:100, T in (Float32, Float64, Int32, Int64)
            if T <: Integer
                q, q1, q2, q3 = [octo(rand((-T(100)):T(100), 8)...) for _ in 1:4]
                c1, c2 = [complex(rand((-T(100)):T(100), 2)...) for _ in 1:2]
            else
                q, q1, q2, q3 = randn(Octonion{T}, 4)
                c1, c2 = randn(Complex{T}, 2)
            end

            # skewfield
            test_group(q1, q2, q3, +, zero(q), -)
            # test all group properties but associativity
            test_neutral(q1, one(q1), *)
            test_inverse(q1, one(q1), *, inv)
            test_multiplicative(q1, q2, *, norm)

            # complex embedding
            test_multiplicative(c1, c2, *, Octonion)
            test_multiplicative(c1, c2, +, Octonion)
        end
    end

    @testset "isreal" begin
        @test isreal(octo(1))
        @test !isreal(octo(1, 1, 0, 0, 0, 0, 0, 0))
        @test !isreal(octo(1, 0, 1, 0, 0, 0, 0, 0))
        @test !isreal(octo(1, 0, 0, 1, 0, 0, 0, 0))
        @test !isreal(octo(1, 0, 0, 0, 1, 0, 0, 0))
        @test !isreal(octo(1, 0, 0, 0, 0, 1, 0, 0))
        @test !isreal(octo(1, 0, 0, 0, 0, 0, 1, 0))
        @test !isreal(octo(1, 0, 0, 0, 0, 0, 0, 1))
    end

    @testset "iszero" begin
        @test iszero(octo(0))
        @test !iszero(octo(1))
        @test !iszero(octo(0, 1, 0, 0, 0, 0, 0, 0))
        @test !iszero(octo(0, 0, 1, 0, 0, 0, 0, 0))
        @test !iszero(octo(0, 0, 0, 1, 0, 0, 0, 0))
        @test !iszero(octo(0, 0, 0, 0, 1, 0, 0, 0))
        @test !iszero(octo(0, 0, 0, 0, 0, 1, 0, 0))
        @test !iszero(octo(0, 0, 0, 0, 0, 0, 1, 0))
        @test !iszero(octo(0, 0, 0, 0, 0, 0, 0, 1))
        @test !iszero(octo(0, 0, 0, 0, 0, 0, 0, 0, true))
    end

    @testset "isone" begin
        @test isone(octo(1))
        @test !isone(octo(-1))
        @test !isone(octo(0, 1, 0, 0, 0, 0, 0, 0))
        @test !isone(octo(1, 1, 0, 0, 0, 0, 0, 0))
        @test !isone(octo(1, 0, 1, 0, 0, 0, 0, 0))
        @test !isone(octo(1, 0, 0, 1, 0, 0, 0, 0))
        @test !isone(octo(1, 0, 0, 0, 1, 0, 0, 0))
        @test !isone(octo(1, 0, 0, 0, 0, 1, 0, 0))
        @test !isone(octo(1, 0, 0, 0, 0, 0, 1, 0))
        @test !isone(octo(1, 0, 0, 0, 0, 0, 0, 1))
    end

    @testset "isfinite" begin
        @test isfinite(octo(1:8...))
        for val in (Inf, -Inf, NaN)
            @test !isfinite(octo(val, 0, 0, 0, 0, 0, 0, 0))
            @test !isfinite(octo(0, val, 0, 0, 0, 0, 0, 0))
            @test !isfinite(octo(0, 0, val, 0, 0, 0, 0, 0))
            @test !isfinite(octo(0, 0, 0, val, 0, 0, 0, 0))
            @test !isfinite(octo(0, 0, 0, 0, val, 0, 0, 0))
            @test !isfinite(octo(0, 0, 0, 0, 0, val, 0, 0))
            @test !isfinite(octo(0, 0, 0, 0, 0, 0, val, 0))
            @test !isfinite(octo(0, 0, 0, 0, 0, 0, 0, val))
            @test isfinite(octo(fill(val, 8)..., true))
        end
    end

    @testset "isinf" begin
        @test !isinf(octo(1, 2, 3, 4, 5, 6, 7, 8))
        @test !isinf(octo(1, 2, 3, 4, 5, 6, 7, NaN))
        for inf in (Inf, -Inf)
            @test isinf(octo(inf, 0, 0, 0, 0, 0, 0, 0))
            @test isinf(octo(0, inf, 0, 0, 0, 0, 0, 0))
            @test isinf(octo(0, 0, inf, 0, 0, 0, 0, 0))
            @test isinf(octo(0, 0, 0, inf, 0, 0, 0, 0))
            @test isinf(octo(0, 0, 0, 0, inf, 0, 0, 0))
            @test isinf(octo(0, 0, 0, 0, 0, inf, 0, 0))
            @test isinf(octo(0, 0, 0, 0, 0, 0, inf, 0))
            @test isinf(octo(0, 0, 0, 0, 0, 0, 0, inf))
            @test !isinf(octo(fill(inf, 8)..., true))
        end
    end

    @testset "isnan" begin
        @test !isnan(octo(1, 2, 3, 4, 5, 6, 7, 8))
        @test !isnan(octo(1, 2, 3, 4, 5, 6, 7, Inf))
        @test !isnan(octo(1, 2, 3, 4, 5, 6, 7, -Inf))
        @test isnan(octo(NaN, 2, 3, 4, 5, 6, 7, 8))
        @test isnan(octo(1, NaN, 3, 4, 5, 6, 7, 8))
        @test isnan(octo(1, 2, NaN, 4, 5, 6, 7, 8))
        @test isnan(octo(1, 2, 3, NaN, 5, 6, 7, 8))
        @test isnan(octo(1, 2, 3, 4, NaN, 6, 7, 8))
        @test isnan(octo(1, 2, 3, 4, 5, NaN, 7, 8))
        @test isnan(octo(1, 2, 3, 4, 5, 6, NaN, 8))
        @test isnan(octo(1, 2, 3, 4, 5, 6, 7, NaN))
    end

    @testset "*" begin
        # verify basic correctness
        q0 = Octonion(1,0,0,0,0,0,0,0)
        q1 = Octonion(0,1,0,0,0,0,0,0)
        q2 = Octonion(0,0,1,0,0,0,0,0)
        q3 = Octonion(0,0,0,1,0,0,0,0)
        q4 = Octonion(0,0,0,0,1,0,0,0)
        q5 = Octonion(0,0,0,0,0,1,0,0)
        q6 = Octonion(0,0,0,0,0,0,1,0)
        q7 = Octonion(0,0,0,0,0,0,0,1)
        @test q0 * q0 == q0
        @test q0 * q1 == q1
        @test q0 * q2 == q2
        @test q0 * q3 == q3
        @test q0 * q4 == q4
        @test q0 * q5 == q5
        @test q0 * q6 == q6
        @test q0 * q7 == q7
        @test q1 * q0 == q1
        @test q1 * q1 == -q0
        @test q1 * q2 == q3
        @test q1 * q3 == -q2
        @test q1 * q4 == -q7
        @test q1 * q5 == -q6
        @test q1 * q6 == q5
        @test q1 * q7 == q4
        @test q2 * q0 == q2
        @test q2 * q1 == -q3
        @test q2 * q2 == -q0
        @test q2 * q3 == q1
        @test q2 * q4 == q6
        @test q2 * q5 == -q7
        @test q2 * q6 == -q4
        @test q2 * q7 == q5
        @test q3 * q0 == q3
        @test q3 * q1 == q2
        @test q3 * q2 == -q1
        @test q3 * q3 == -q0
        @test q3 * q4 == -q5
        @test q3 * q5 == q4
        @test q3 * q6 == -q7
        @test q3 * q7 == q6
        @test q4 * q0 == q4
        @test q4 * q1 == q7
        @test q4 * q2 == -q6
        @test q4 * q3 == q5
        @test q4 * q4 == -q0
        @test q4 * q5 == -q3
        @test q4 * q6 == q2
        @test q4 * q7 == -q1
        @test q5 * q0 == q5
        @test q5 * q1 == q6
        @test q5 * q2 == q7
        @test q5 * q3 == -q4
        @test q5 * q4 == q3
        @test q5 * q5 == -q0
        @test q5 * q6 == -q1
        @test q5 * q7 == -q2
        @test q6 * q0 == q6
        @test q6 * q1 == -q5
        @test q6 * q2 == q4
        @test q6 * q3 == q7
        @test q6 * q4 == -q2
        @test q6 * q5 == q1
        @test q6 * q6 == -q0
        @test q6 * q7 == -q3
        @test q7 * q0 == q7
        @test q7 * q1 == -q4
        @test q7 * q2 == -q5
        @test q7 * q3 == -q6
        @test q7 * q4 == q1
        @test q7 * q5 == q2
        @test q7 * q6 == q3
        @test q7 * q7 == -q0

        @testset "* same between Quaternions and Octonion" begin
            # make sure this tracks the `*` tests for Quaternions
            q1 = Octonion(1,0,0,0,0,0,0,0)
            qi = Octonion(0,1,0,0,0,0,0,0)
            qj = Octonion(0,0,1,0,0,0,0,0)
            qk = Octonion(0,0,0,1,0,0,0,0)
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
    end

    @testset "abs2" for _ in 1:100, T in (Float16, Float32, Float64)
        o = rand(Octonion{T})
        @test abs2(o) == o'*o
    end

    @testset "/" begin
        for _ in 1:100
            o, o2 = randn(OctonionF64, 2)
            x = randn()
            @test o / o ≈ o \ o ≈ one(o)
            @test o / o2 ≈ o * inv(o2)
            @test o2 \ o ≈ inv(o2) * o
            @test o / x ≈ x \ o ≈ inv(x) * o
        end
    end

    @testset "^" begin
        @testset "^(::Octonion, ::Real)" begin
            for _ in 1:100
                o = randn(OctonionF64)
                @test @inferred(o^2.0) ≈ o * o
                @test o^1.0 ≈ o
                @test o^-1.0 ≈ inv(o)
                @test o^1.3 ≈ exp(1.3 * log(o))
                @test o^7.8 ≈ exp(7.8 * log(o))
                @test o^1.3f0 ≈ exp(1.3f0 * log(o))
                @test o^7.8f0 ≈ exp(7.8f0 * log(o))
            end
        end
        @testset "^(::Octonion, ::Octonion)" begin
            @test octo(Float64(ℯ))^octo(0, 0, 0, 0, 0, 0, π / 2, 0) ≈
                octo(0, 0, 0, 0, 0, 0, 1, 0)
            z = (3.5 + 2.3im)^(0.2 + 1.7im)
            @test octo(3.5, 0, 0, 0, 0, 0, 2.3, 0)^octo(0.2, 0, 0, 0, 0, 0, 1.7, 0) ≈
                octo(real(z), 0, 0, 0, 0, 0, imag(z), 0)
            for _ in 1:100
                q, p = randn(OctonionF64, 2)
                @test @inferred(q^p) ≈ exp(p * log(q))
            end
        end
    end

    @testset "non-analytic functions" begin
        unary_funs = [conj, abs, abs2, norm, sign]
        # since every octonion is conjugate to a quaternion,
        # one can establish correctness as follows:
        @testset for fun in unary_funs
            for _ in 1:100
                o1, o2 = randn(OctonionF64, 2)
                q = randn(QuaternionF64)
                o = octo(q)
                @test @inferred(fun(o)) ≈ fun(q)
                @test o2 * fun(o1) * inv(o2) ≈ fun(o2 * o1 * inv(o2))
            end
        end
    end

    @testset "analytic functions" begin
        unary_funs = [sqrt, inv, exp, log]
        # since every octonion is conjugate to a quaternion,
        # one can establish correctness as follows:
        @testset for fun in unary_funs
            for _ in 1:100
                o1, o2 = randn(OctonionF64, 2)
                q = randn(QuaternionF64)
                o = octo(q)
                @test @inferred(fun(o)) ≈ fun(q)
                @test o2 * fun(o1) * inv(o2) ≈ fun(o2 * o1 * inv(o2))
            end
        end

        @testset "identities" begin
            for _ in 1:100
                o = randn(OctonionF64)
                @test inv(o) * o ≈ o * inv(o) ≈ one(o)
                @test sqrt(o) * sqrt(o) ≈ o
                @test exp(log(o)) ≈ o
                @test exp(zero(o)) === one(o)
                @test log(one(o)) === zero(o)
            end
            @test log(zero(OctonionF64)) === octo(-Inf)
        end
    end

    @testset "normalize" begin
        for _ in 1:100
            q = randn(OctonionF64)
            qnorm = @inferred normalize(q)
            @test abs(qnorm) ≈ 1
            @test qnorm.norm
            @test q ≈ abs(q) * qnorm
            @test normalize(qnorm) === qnorm
        end
        @test_broken @inferred(normalize(octo(1:8...)))
    end

    @testset "normalizea" begin
        for _ in 1:100
            q = randn(OctonionF64)
            qnorm, a = @inferred normalizea(q)
            @test abs(qnorm) ≈ 1
            @test qnorm.norm
            @test a isa Real
            @test a ≈ abs(q)
            @test q ≈ a * qnorm
            @test normalizea(qnorm) === (qnorm, one(real(q)))
        end
        @test_broken @inferred(normalizea(octo(1:8...)))
    end
end
