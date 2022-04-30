using Quaternions
using DualNumbers
using ForwardDiff
using LinearAlgebra
using Random
using Test

function fdquat_to_dualquat(q::Quaternion{<:ForwardDiff.Dual})
    q0 = quat(
        ForwardDiff.value(q.s),
        ForwardDiff.value(q.v1),
        ForwardDiff.value(q.v2),
        ForwardDiff.value(q.v3),
    )
    qe = quat(
        ForwardDiff.partials(q.s, 1),
        ForwardDiff.partials(q.v1, 1),
        ForwardDiff.partials(q.v2, 1),
        ForwardDiff.partials(q.v3, 1),
    )
    return dualquat(q0, qe)
end

function dualquat_to_fdquat(dq::DualQuaternion)
    return quat(
        ForwardDiff.Dual(dq.q0.s, dq.qe.s),
        ForwardDiff.Dual(dq.q0.v1, dq.qe.v1),
        ForwardDiff.Dual(dq.q0.v2, dq.qe.v2),
        ForwardDiff.Dual(dq.q0.v3, dq.qe.v3),
    )
end

@testset "DualQuaternion" begin
    @testset "type aliases" begin
        @test DualQuaternionF16 === DualQuaternion{Float16}
        @test DualQuaternionF32 === DualQuaternion{Float32}
        @test DualQuaternionF64 === DualQuaternion{Float64}
    end

    @testset "Constructors" begin
        @testset "from quaternions" begin
            q0 = Quaternion{Int}(1, 2, 3, 4, false)
            qe = QuaternionF32(5, 6, 7, 8, false)
            @test @inferred(DualQuaternion(q0, qe)) isa DualQuaternionF32
            @test DualQuaternion(q0, qe) === DualQuaternionF32(
                QuaternionF32(1, 2, 3, 4, false), QuaternionF32(5, 6, 7, 8, false), false
            )
            @test @inferred(DualQuaternionF64(q0, qe, false)) === DualQuaternionF64(
                QuaternionF64(1, 2, 3, 4, false), QuaternionF64(5, 6, 7, 8, false), false
            )
            @test DualQuaternionF64(q0, qe, true) === DualQuaternionF64(
                QuaternionF64(1, 2, 3, 4, false), QuaternionF64(5, 6, 7, 8, false), true
            )

            qnorm = Quaternion(0//1, 1//1, 0//1, 0//1, true)
            @test @inferred(DualQuaternion(q0)) isa DualQuaternion{Int}
            @test DualQuaternion(qnorm) === DualQuaternion(qnorm, zero(qnorm), true)
            @test @inferred(DualQuaternionF64(q0)) ===
                DualQuaternionF64(QuaternionF64(q0), zero(QuaternionF64), false)
        end
        @testset "from dual" begin
            coef = (dual(1, 2), dual(3.0, 4.0), dual(5.0f0, 6.0f0), dual(7//1, 8//1))
            @testset for T in (Float32, Float64, Int), norm in (true, false)
                dq = @inferred DualQuaternion{T}(coef..., norm)
                @test dq isa DualQuaternion{T}
                @test dq.norm === norm
                @test dq === DualQuaternion{T}(convert.(Dual{T}, coef)..., norm)
                @test dq == DualQuaternion(
                    Quaternion(DualNumbers.value.(coef)...),
                    Quaternion(DualNumbers.epsilon.(coef)...),
                )
                dq2 = @inferred DualQuaternion(convert.(Dual{T}, coef)..., norm)
                @test DualQuaternion(convert.(Dual{T}, coef)..., norm) === dq
                @test DualQuaternion(coef[1]) ==
                    DualQuaternion(coef[1], fill(zero(coef[1]), 3)...)
                @test DualQuaternion{T}(coef[1]) ==
                    DualQuaternion(convert(Dual{T}, coef[1]), fill(zero(Dual{T}), 3)...)
                if !norm
                    @test DualQuaternion(convert.(Dual{T}, coef)...) === dq
                end
            end
        end
        @testset "from real" begin
            @testset for x in (-1//1, 1.0, 2.0), T in (Float32, Float64, Int, Rational{Int})
                coef = (Quaternion{T}(x, 0, 0, 0, isone(abs(x))), zero(Quaternion{T}))
                @test @inferred(DualQuaternion{T}(x)) ===
                    DualQuaternion{T}(coef..., isone(abs(x)))
                @test @inferred(DualQuaternion(T(x))) ===
                    DualQuaternion{T}(coef..., isone(abs(x)))
            end
        end
        @testset "from dual quaternion" begin
            dq = DualQuaternion(
                QuaternionF32(1, 2, 3, 4, false), QuaternionF32(4, 5, 6, 7, false), false
            )
            @test @inferred(DualQuaternion(dq)) === dq
            @test @inferred(DualQuaternionF64(dq)) ===
                DualQuaternionF64(dq.q0, dq.qe, false)
        end
        @testset "from vector" begin
            v = randn(3)
            @test @inferred(DualQuaternion(v)) ===
                DualQuaternion(Quaternion(zero(v)), Quaternion(v))
        end
    end

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
        @test promote(DualQuaternion(1.0), 1.0) ===
            (DualQuaternion(1.0), DualQuaternion(1.0))
        @test promote(DualQuaternion(1.0f0), 2.0) ===
            (DualQuaternion(1.0), DualQuaternion(2.0))
        @test promote(DualQuaternion(1.0f0), dual(1, 2)) ===
            (DualQuaternion(1.0f0), DualQuaternion(dual(1.0f0, 2.0f0)))
        @test promote(DualQuaternion(1.0f0), Quaternion(3//1)) ===
            (DualQuaternion(1.0f0), DualQuaternion(3.0f0))
        @test promote(DualQuaternion(1.0f0), DualQuaternion(2.0)) ===
            (DualQuaternion(1.0), DualQuaternion(2.0))

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
        @testset "dualquatrand" begin
            dq = dualquatrand()
            @test dq isa DualQuaternionF64
            @test !dq.norm
        end

        @testset "ndualquatrand" begin
            dq = ndualquatrand()
            @test dq isa DualQuaternionF64
            @test dq.norm
        end

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
                    imag_part(dq.q0)...
                    real(dq.qe)
                    imag_part(dq.qe)...
                ]
            end
            xs_mean = sum(xs) / length(xs)
            xs_var = sum(x -> abs2.(x .- xs_mean), xs) / (length(xs) - 1)
            @test all(isapprox.(xs_mean, 0.5; atol=0.1))
            @test all(isapprox.(xs_var, 1 / 12; atol=0.01))
        end
    end

    @testset "basic" begin
        q = rand(DualQuaternionF64)
        qnorm = normalize(q)
        @test_throws MethodError imag(q)
        @test_throws MethodError Quaternions.imag(q)
        @test_throws MethodError imag_part(q)
        @test_throws MethodError reim(q)
        @test conj(q) === dualquat(conj(q.q0), conj(q.qe), q.norm)
        @test conj(qnorm) === dualquat(conj(qnorm.q0), conj(qnorm.qe), qnorm.norm)
        @test conj(conj(q)) === q
        @test conj(conj(qnorm)) === qnorm
        @test float(dualquat(dual.(1:4, 5:8)...)) === dualquat(dual.(1.0:4.0, 5.0:8.0)...)
        @test dconj(q) === dualquat(q.q0, -q.qe, q.norm)
        @test dconj(qnorm) === dualquat(qnorm.q0, -qnorm.qe, qnorm.norm)
        @test Quaternions.Q0(q) === q.q0
        @test Quaternions.Qe(q) === q.qe
    end

    @testset "algebraic properties" begin
        @testset "addition/subtraction" begin
            for _ in 1:100
                dq1, dq2 = rand(DualQuaternionF64, 2)
                @test (dq1 + dq2).q0 ≈ dq1.q0 + dq2.q0
                @test (dq1 + dq2).qe ≈ dq1.qe + dq2.qe
                @test (dq1 - dq2).q0 ≈ dq1.q0 - dq2.q0
                @test (dq1 - dq2).qe ≈ dq1.qe - dq2.qe
                @test (+dq1).q0 === dq1.q0
                @test (+dq1).qe === dq1.qe
                @test (-dq1).q0 === -dq1.q0
                @test (-dq1).qe === -dq1.qe
            end
        end
        @testset "division" begin
            for _ in 1:100
                dq = rand(DualQuaternionF64)
                dq2 = dq / dq
                @test dq2.q0 ≈ 1
                @test dq2.qe ≈ zero(dq2.qe) atol = 1e-6
                dq2 = dq \ dq
                @test dq2.q0 ≈ 1
                @test dq2.qe ≈ zero(dq2.qe) atol = 1e-6
                s = randn()
                @test (dq / s).q0 ≈ dq.q0 / s
                @test (dq / s).qe ≈ dq.qe / s
                d = dual(randn(2)...)
                @test (dq / d).q0 ≈ (dq / dualquat(d)).q0
                @test (dq / d).qe ≈ (dq / dualquat(d)).qe
            end
        end
        @testset "multiplication is associative" begin
            for _ in 1:100
                dq1, dq2, dq3 = rand(DualQuaternionF64, 3)
                dq4 = (dq1 * dq2) * dq3
                dq5 = dq1 * (dq2 * dq3)
                @test dq4.q0 ≈ dq5.q0
                @test dq4.qe ≈ dq5.qe
            end
        end
        @testset "multiplication by dual" begin
            for _ in 1:100
                dq = rand(DualQuaternionF64)
                d = dual(randn(2)...)
                @test (dq * d).q0 ≈ (dq * dualquat(d)).q0
                @test (dq * d).qe ≈ (dq * dualquat(d)).qe
                @test (d * dq).q0 ≈ (dualquat(d) * dq).q0
                @test (d * dq).qe ≈ (dualquat(d) * dq).qe
            end
        end
    end

    @testset "iszero" begin
        @test iszero(dualquat(0))
        @test !iszero(dualquat(1))
        @test !iszero(dualquat(quat(0, 1, 0, 0)))
        @test !iszero(dualquat(quat(0, 0, 1, 0)))
        @test !iszero(dualquat(quat(0, 0, 0, 1)))
        @test !iszero(dualquat(quat(0), quat(1)))
        @test !iszero(dualquat(quat(0), quat(0, 1, 0, 0)))
        @test !iszero(dualquat(quat(0), quat(0, 0, 1, 0)))
        @test !iszero(dualquat(quat(0), quat(0, 0, 0, 1)))
    end

    @testset "isone" begin
        @test isone(dualquat(1))
        @test !isone(dualquat(-1))
        @test !isone(dualquat(quat(0, 1, 0, 0)))
        @test !isone(dualquat(quat(1, 1, 0, 0)))
        @test !isone(dualquat(quat(1, 0, 1, 0)))
        @test !isone(dualquat(quat(1, 0, 0, 1)))
        @test !isone(dualquat(quat(0), quat(1)))
        @test !isone(dualquat(quat(0), quat(0, 1, 0, 0)))
        @test !isone(dualquat(quat(0), quat(0, 0, 1, 0)))
        @test !isone(dualquat(quat(0), quat(0, 0, 0, 1)))
    end

    @testset "^" begin
        @testset "^(::DualQuaternion, ::Real)" begin
            for _ in 1:100
                dq = rand(DualQuaternionF64)
                @test dq^2 === dq * dq
                @test dq^1 === dq
                @test_broken (dq^2.0).q0 ≈ (dq * dq).q0
                @test_broken (dq^2.0).qe ≈ (dq * dq).qe
                @test_broken (dq^1.0).q0 ≈ dq.q0
                @test_broken (dq^1.0).qe ≈ dq.qe
                @test_broken (dq^-1.0).q0 ≈ inv(dq).q0
                @test_broken (dq^-1.0).qe ≈ inv(dq).qe
                for p in (1.3, 7.8, 1.3f0, 7.8f0)
                    @test (dq^p).q0 ≈ exp(p * log(dq)).q0
                    @test (dq^p).qe ≈ exp(p * log(dq)).qe
                end
            end
        end
        @testset "^(::DualQuaternion, ::DualQuaternion)" begin
            for _ in 1:100
                q, p = randn(QuaternionF64, 2)
                @test q^p ≈ exp(p * log(q))
            end
        end
    end

    @testset "non-analytic functions" begin
        unary_funs = [conj, abs, abs2, norm, sign]
        # since every dual quaternion is conjugate to a dual complex number,
        # one can establish correctness as follows:
        @testset for fun in unary_funs
            for _ in 1:100
                dq1, dq2 = rand(DualQuaternionF64, 2)
                c = dual(randn(ComplexF64), randn(ComplexF64))
                dq = dualquat(reim(c)..., reim(zero(c))...)
                p = dualquat(@inferred(fun(dq)))
                @test p.q0 ≈ DualNumbers.value(fun(c))
                @test p.qe ≈ DualNumbers.epsilon(fun(c))
                p2 = dq2 * fun(dq1) * inv(dq2)
                p3 = dualquat(fun(dq2 * dq1 * inv(dq2)))
                @test p2.q0 ≈ p3.q0
                @test p2.qe ≈ p3.qe
            end
        end
    end

    @testset "analytic functions" begin
        unary_funs = [sqrt, inv, exp, log]
        # since every dual quaternion is conjugate to a dual complex number,
        # one can establish correctness as follows:
        @testset for fun in unary_funs
            (fun === log || fun === sqrt) && continue  # log is currently broken
            for _ in 1:100
                dq1, dq2 = rand(DualQuaternionF64, 2)
                c = dual(randn(ComplexF64), randn(ComplexF64))
                dq = dualquat(reim(c)..., reim(zero(c))...)
                p = dualquat(@inferred(fun(dq)))
                @test p.q0 ≈ DualNumbers.value(fun(c))
                @test p.qe ≈ DualNumbers.epsilon(fun(c))
                p2 = dq2 * fun(dq1) * inv(dq2)
                p3 = dualquat(fun(dq2 * dq1 * inv(dq2)))
                @test p2.q0 ≈ p3.q0
                @test p2.qe ≈ p3.qe
            end
        end

        @testset "identities" begin
            for _ in 1:100
                dq = rand(DualQuaternionF64)
                @test (inv(dq) * dq).q0 ≈ one(dq.q0)
                @test (inv(dq) * dq).qe ≈ zero(dq.qe) atol = 1e-6
                @test (dq * inv(dq)).q0 ≈ one(dq.q0)
                @test (dq * inv(dq)).qe ≈ zero(dq.qe) atol = 1e-6
                @test_broken (sqrt(dq) * sqrt(dq)).q0 ≈ dq.q0
                @test_broken (sqrt(dq) * sqrt(dq)).qe ≈ dq.qe
                @test_broken exp(log(dq)).q0 ≈ dq.q0
                @test_broken exp(log(dq)).qe ≈ dq.qe
                @test exp(zero(dq)).q0 ≈ one(dq.q0)
                @test exp(zero(dq)).qe ≈ zero(dq.qe) atol = 1e-6
                @test log(one(dq)).q0 ≈ zero(dq.q0) atol = 1e-6
                @test log(one(dq)).qe ≈ zero(dq.qe) atol = 1e-6
            end
            @test_broken log(zero(DualQuaternionF64)) === dualquat(-Inf)
        end
    end

    @testset "normalize" begin
        for _ in 1:100
            dq = rand(DualQuaternionF64)
            q = dualquat_to_fdquat(dq)
            dqnorm = normalize(dq)
            qnorm = normalize(q)
            dqnorm2 = fdquat_to_dualquat(qnorm)
            @test dqnorm.q0 ≈ dqnorm2.q0
            @test dqnorm.qe ≈ dqnorm2.qe
            @test normalize(dqnorm) === dqnorm
            @test iszero(normalize(zero(dq)))
        end
    end

    @testset "normalizea" begin
        for _ in 1:100
            dq = rand(DualQuaternionF64)
            q = dualquat_to_fdquat(dq)
            dqnorm, dqa = normalizea(dq)
            qnorm, qa = normalizea(q)
            dqnorm2 = fdquat_to_dualquat(qnorm)
            @test dqnorm.q0 ≈ dqnorm2.q0
            @test dqnorm.qe ≈ dqnorm2.qe
            @test DualNumbers.value(dqa) ≈ ForwardDiff.value(qa)
            @test DualNumbers.epsilon(dqa) ≈ ForwardDiff.partials(qa, 1)
        end
    end

    @testset "angle/axis/angleaxis" begin
        dq1 = rand(DualQuaternionF64)
        dq2 = dq1 + (1 - dq1.q0.s)
        for dq in (dq1, dq2)
            t1, p1 = angleaxis(dq)
            t2, v = angleaxis(dualquat_to_fdquat(dq))
            p2 = fdquat_to_dualquat(quat(v))
            @test DualNumbers.value(t1) ≈ ForwardDiff.value(t2)
            @test_broken DualNumbers.epsilon(t1) ≈ ForwardDiff.partials(t2, 1)
            @test p2.q0 ≈ p1.q0
            @test_broken p2.qe ≈ p1.qe
            t3 = angle(dq)
            @test DualNumbers.value(t3) ≈ DualNumbers.value(t1)
            @test DualNumbers.epsilon(t3) ≈ DualNumbers.epsilon(t1)
            p3 = axis(dq)
            @test p3.q0 ≈ p1.q0
            @test p3.qe ≈ p1.qe
        end
    end
end
