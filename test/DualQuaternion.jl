using Quaternions
using DualNumbers
using LinearAlgebra
using Random
using Test

@testset "DualQuaternion" begin
    @testset "type aliases" begin
        @test DualQuaternionF16 === DualQuaternion{Float16}
        @test DualQuaternionF32 === DualQuaternion{Float32}
        @test DualQuaternionF64 === DualQuaternion{Float64}
    end

    @testset "Constructors" begin end

    @testset "==" begin
        @test DualQuaternion(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) ==
            DualQuaternion(Quaternion(1.0, 2, 3, 4), Quaternion(5, 6, 7, 8))
        @test DualQuaternion(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) !=
            DualQuaternion(Quaternion(1.0, 2, 3, 4), Quaternion(1, 2, 3, 4))
    end

    @testset "convert" begin
        @test convert(DualQuaternion{Float64}, 1) === DualQuaternion(1.0)
        @test convert(DualQuaternion{Float64}, DualNumbers.Dual(1, 2)) ===
            DualQuaternion(Quaternion(1.0), Quaternion(2.0))
        @test convert(DualQuaternion{Float64}, Quaternion(1, 2, 3, 4)) ===
            DualQuaternion(Quaternion(1.0, 2.0, 3.0, 4.0))
        @test convert(
            DualQuaternion{Float64},
            DualQuaternion(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)),
        ) === DualQuaternion(Quaternion(1.0, 2.0, 3.0, 4.0), Quaternion(5.0, 6.0, 7.0, 8.0))
    end

    @testset "promote" begin
        @test Quaternion(1, 2, 3, 4) == DualQuaternion(Quaternion(1, 2, 3, 4))
        @test Quaternion(1, 2, 3, 4) !=
            DualQuaternion(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8))
        @test DualQuaternion(1) == 1.0
    end

    @testset "shorthands" begin
        @test dualquat(Quaternion(1, 0, 0, 0)) == Quaternion(1, 0, 0, 0)
        @test dualquat(Quaternion(1, 2, 3, 4)) == Quaternion(1, 2, 3, 4)
        @test dualquat(Quaternion(1, 0, 0, 0)) === DualQuaternion(Quaternion(1, 0, 0, 0)) # checking the .norm field in particular
        @test dualquat(Quaternion(1, 2, 3, 4)) === DualQuaternion(Quaternion(1, 2, 3, 4))
        @test dualquat(1) === DualQuaternion(1)
        @test dualquat(Dual(1, 2)) === DualQuaternion(Dual(1, 2))
        @test dualquat(Dual(1, 2), Dual(0), Dual(0), Dual(0)) ===
            DualQuaternion(Dual(1, 2), Dual(0), Dual(0), Dual(0))
        @test dualquat(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) ==
            DualQuaternion(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8))
        @test dualquat(Quaternion(1, 0, 0, 0), Quaternion(0)).norm == false
        @test dualquat(Quaternion(1, 0, 0, 0), Quaternion(0), false).norm == false # respect the .norm input (even if wrong)
        @test dualquat(Quaternion(1, 2, 3, 4), Quaternion(0), true).norm == true # respect the .norm input (even if wrong)
        @test dualquat(Dual(2, 0), Dual(0), Dual(0), Dual(0), true).norm == true # respect the .norm input (even if wrong)
        @test dualquat(Dual(1, 0), Dual(0), Dual(0), Dual(0), false).norm == false # respect the .norm input (even if wrong)
    end

    @testset "random generation" begin
        @testset "rand($H)" for H in (DualQuaternionF32, DualQuaternionF64)
            rng = Random.MersenneTwister(42)
            dq = rand(rng, H)
            @test dq isa H
            @test !dq.norm

            dqs = rand(rng, H, 1000)
            @test eltype(dqs) === H
            @test length(dqs) == 1000
            xs = map(dqs) do dq
                return [
                    real(dq.q0)
                    Quaternions.imag(dq.q0)
                    real(dq.qe)
                    Quaternions.imag(dq.qe)
                ]
            end
            xs_mean = sum(xs) / length(xs)
            xs_var = sum(x -> abs2.(x .- xs_mean), xs) / (length(xs) - 1)
            @test all(isapprox.(xs_mean, 0.5; atol=0.1))
            @test all(isapprox.(xs_var, 1 / 12; atol=0.01))
        end
    end

    @testset "basic" begin
        # TODO: test real, imag, conj, dconj, float, and Quaternions.abs_imag
    end

    @testset "algebraic properties" begin end

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
