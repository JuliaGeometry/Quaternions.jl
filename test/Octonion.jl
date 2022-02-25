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

    @testset "Constructors" begin end

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
        @test convert(Octonion{Float64}, Octonion(1, 2, 3, 4, 5, 6, 7, 8)) ===
            Octonion(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0)
    end

    @testset "promote" begin
        @test Quaternion(1, 2, 3, 4) == Octonion(1, 2, 3, 4, 0, 0, 0, 0)
        @test Quaternion(1, 2, 3, 4) != Octonion(1, 2, 3, 4, 5, 6, 7, 8)
        @test Octonion(1) == 1.0
        @test Octonion(Complex(1, 2)) == Complex(1, 2)
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
