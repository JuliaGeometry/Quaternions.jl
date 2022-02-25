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
            cs = [
                (1, 2.0, 3f0, 4//1, 5, 6, 7, 8),
                (1//1, 2f0, 3f0, 4, 5, 6, 7, 8),
            ]
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
            @testset for z in (1+0im, -im, 1+2im), T in (Float32, Float64, Int, Rational{Int})
                coef = T.((reim(z)..., zeros(6)...))
                z2 = Complex{T}(z)
                norm = isone(abs(z))
                @test Octonion{T}(z) === Octonion{T}(coef..., norm)
                @test @inferred(Octonion(z2)) === Octonion{T}(coef..., norm)
            end
        end
        @testset "from quaternion" begin
            qs = (Quaternion(1,2,3,4), QuaternionF64(0,1,0,0,true))
            @testset for q in qs, T in (Float32,Float64)
                coef = T.((q.s, q.v1, q.v2, q.v3, zeros(4)...))
                q2 = Quaternion{T}(q)
                @test @inferred(Octonion{T}(q)) === Octonion{T}(coef..., q.norm)
                @test @inferred(Octonion(q2)) === Octonion{T}(coef..., q.norm)
            end
        end
        @testset "from octonion" begin
            os = (
                Octonion(1,2,3,4,5,6,7,8),
                OctonionF64(0,1,0,0,0,0,0,0,true)
            )
            @testset for o in os, T in (Float32,Float64)
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
        @test promote(Octonion(1.0, 2:8...), 1.0) ===
            (Octonion(1.0, 2:8...), Octonion(1.0))
        @test promote(Octonion(1f0, 2:8...), 2.0) === (Octonion(1.0, 2:8...), Octonion(2.0))
        @test promote(Octonion(1f0), 2+3im) === (Octonion(1f0), Octonion(2f0+3f0im))
        @test promote(Octonion(1f0), Quaternion(1,2,3,4)) ===
            (Octonion(1f0), Octonion(1f0:4f0..., fill(0, 4)...))
        @test promote(Octonion(1f0), Octonion(2.0)) === (Octonion(1.0), Octonion(2.0))

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
    end

    @testset "random generation" begin
        @testset "rand($H)" for H in (OctonionF32, OctonionF64)
            rng = Random.MersenneTwister(42)
            o = rand(rng, H)
            @test o isa H
            @test !o.norm

            os = rand(rng, H, 1000)
            @test eltype(os) === H
            @test length(os) == 1000
            xs = map(os) do o
                return [real(o); Quaternions.imag(o)]
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
                return [real(o); Quaternions.imag(o)]
            end
            xs_mean = sum(xs) / length(xs)
            xs_var = sum(x -> abs2.(x .- xs_mean), xs) / (length(xs) - 1)
            @test all(isapprox.(xs_mean, 0; atol=0.1))
            @test all(isapprox.(xs_var, 1 / 8; atol=0.1))
        end
    end

    @testset "basic" begin
        # TODO: test real, imag, conj, float, and Quaternions.abs_imag
    end

    @testset "algebraic properties" begin
        for _ in 1:10, T in (Float32, Float64, Int32, Int64)
            if T <: Integer
                q, q1, q2, q3 = [Quaternion(rand((-T(100)):T(100), 4)...) for _ in 1:4]
                c1, c2 = [complex(rand((-T(100)):T(100), 2)...) for _ in 1:2]
            else
                q, q1, q2, q3 = randn(Quaternion{T}, 4)
                c1, c2 = randn(Complex{T}, 2)
            end

            # skewfield
            test_group(q1, q2, q3, +, zero(q), -)
            test_group(q1, q2, q3, *, one(q), inv)
            test_multiplicative(q1, q2, *, norm)

            # complex embedding
            test_multiplicative(c1, c2, *, Quaternion)
            test_multiplicative(c1, c2, +, Quaternion)
        end
    end

    @testset "isreal" begin end

    @testset "iszero" begin end

    @testset "isone" begin end

    @testset "isfinite" begin end

    @testset "isinf" begin end

    @testset "isnan" begin end

    @testset "+" begin end

    @testset "-" begin end

    @testset "*" begin end

    @testset "/" begin end

    @testset "^" begin end

    @testset "non-analytic functions" begin end

    @testset "analytic functions" begin end

    @testset "normalize" begin end

    @testset "normalizea" begin end
end
