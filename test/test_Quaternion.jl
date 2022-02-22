
using Quaternions: argq
using DualNumbers
using LinearAlgebra
using Random

@testset "type aliases" begin
    @test QuaternionF16 === Quaternion{Float16}
    @test QuaternionF32 === Quaternion{Float32}
    @test QuaternionF64 === Quaternion{Float64}
    @test OctonionF16 === Octonion{Float16}
    @test OctonionF32 === Octonion{Float32}
    @test OctonionF64 === Octonion{Float64}
    @test DualQuaternionF16 === DualQuaternion{Float16}
    @test DualQuaternionF32 === DualQuaternion{Float32}
    @test DualQuaternionF64 === DualQuaternion{Float64}
end

# creating random examples
sample(QT::Type{Quaternion{T}}) where {T <: Integer} = QT(rand(-100:100, 4)..., false)
sample(QT::Type{Quaternion{T}}) where {T <: AbstractFloat} = QT(rand(Bool) ? quatrand() : nquatrand())
sample(CT::Type{Complex{T}}) where {T <: Integer} = CT(rand(-100:100, 2)...)
sample(CT::Type{Complex{T}}) where {T <: AbstractFloat} = CT(randn(2)...)
sample(T, n) = T[sample(T) for _ in 1:n]

# test algebraic properties of quaternions
for _ in 1:10, T in (Float32, Float64, Int32, Int64)
    q, q1, q2, q3 = sample(Quaternion{T}, 4)

    # skewfield
    test_group(q1, q2, q3, +, zero(q), -)
    test_group(q1, q2, q3, *, one(q), inv)
    test_multiplicative(q1, q2, *, norm)

    # complex embedding
    c1, c2 = sample(Complex{T}, 2)
    test_multiplicative(c1, c2, *, Quaternion)
    test_multiplicative(c1, c2, +, Quaternion)
end

@testset "promotions and equalities" begin
    @test Quaternion(1,0,0,0,false) == Quaternion(1,0,0,0,true) # test that .norm field does not affect equality
    @test Quaternion(1) == 1.0 # test promotion
    @test Quaternion(1,2,0,0) == Complex(1.0,2.0) # test promotion
    @test Quaternion{Float64}(1) === Quaternion(1.0) # explicit type construction
    @test quat(1) === Quaternion(1) # checking the .norm field in particular
    @test quat(1,0,0,0) === Quaternion(1,0,0,0) # checking the .norm field in particular
    @test quat(1,2,3,4) === Quaternion(1,2,3,4)
    @test quat(Quaternion(1,0,0,0)) === Quaternion(1,0,0,0) # checking the .norm field in particular
    @test quat(Quaternion(1,2,3,4)) === Quaternion(1,2,3,4)
    @test quat(1,0,0,0,false).norm == false # respect the .norm input (even if wrong)
    @test quat(1,2,3,4,true).norm == true # respect the .norm input (even if wrong)

    @test Quaternion(1,2,3,4) == DualQuaternion(Quaternion(1,2,3,4))
    @test Quaternion(1,2,3,4) != DualQuaternion(Quaternion(1,2,3,4),Quaternion(5,6,7,8))
    @test DualQuaternion(1) == 1.0 # test promotion
    @test DualQuaternion(Quaternion(1,2,3,4),Quaternion(5,6,7,8)) == DualQuaternion(Quaternion(1.0,2,3,4),Quaternion(5,6,7,8))
    @test DualQuaternion(Quaternion(1,2,3,4),Quaternion(5,6,7,8)) != DualQuaternion(Quaternion(1.0,2,3,4),Quaternion(1,2,3,4))
    @test dualquat(Quaternion(1,0,0,0)) == Quaternion(1,0,0,0)
    @test dualquat(Quaternion(1,2,3,4)) == Quaternion(1,2,3,4)
    @test dualquat(Quaternion(1,0,0,0)) === DualQuaternion(Quaternion(1,0,0,0)) # checking the .norm field in particular
    @test dualquat(Quaternion(1,2,3,4)) === DualQuaternion(Quaternion(1,2,3,4))
    @test dualquat(1) === DualQuaternion(1)
    @test dualquat(Dual(1,2)) === DualQuaternion(Dual(1,2))
    @test dualquat(Dual(1,2),Dual(0),Dual(0),Dual(0)) === DualQuaternion(Dual(1,2),Dual(0),Dual(0),Dual(0))
    @test dualquat(Quaternion(1,2,3,4),Quaternion(5,6,7,8)) == DualQuaternion(Quaternion(1,2,3,4),Quaternion(5,6,7,8))
    @test dualquat(Quaternion(1,0,0,0),Quaternion(0)).norm == false
    @test dualquat(Quaternion(1,0,0,0),Quaternion(0),false).norm == false # respect the .norm input (even if wrong)
    @test dualquat(Quaternion(1,2,3,4),Quaternion(0),true).norm == true # respect the .norm input (even if wrong)
    @test dualquat(Dual(2,0),Dual(0),Dual(0),Dual(0),true).norm == true # respect the .norm input (even if wrong)
    @test dualquat(Dual(1,0),Dual(0),Dual(0),Dual(0),false).norm == false # respect the .norm input (even if wrong)

    @test Quaternion(1,2,3,4) == Octonion(1,2,3,4,0,0,0,0)
    @test Quaternion(1,2,3,4) != Octonion(1,2,3,4,5,6,7,8)
    @test Octonion(1,0,0,0,0,0,0,0,false) == Octonion(1,0,0,0,0,0,0,0,true) # test that .norm field does not affect equality
    @test Octonion(1) == 1.0 # test promotion
    @test Octonion(Complex(1,2)) == Complex(1,2)
    @test Octonion(1.0,2,3,4,5,6,7,8) == Octonion(1,2,3,4,5,6,7,8)
    @test Octonion(1.0,2,3,4,5,6,7,8) != Octonion(1,2,3,4,1,2,3,4)
    @test octo(1) === Octonion(1) # checking the .norm field in particular
    @test octo(1,0,0,0,0,0,0,0) === Octonion(1,0,0,0,0,0,0,0) # checking the .norm field in particular
    @test octo(1,2,3,4,5,6,7,8) === Octonion(1,2,3,4,5,6,7,8)
    @test octo(Octonion(1,0,0,0,0,0,0,0)) === Octonion(1,0,0,0,0,0,0,0) # checking the .norm field in particular
    @test octo(Octonion(1,2,3,4,5,6,7,8)) === Octonion(1,2,3,4,5,6,7,8)
    @test octo(1,0,0,0,0,0,0,0,false).norm == false # respect the .norm input (even if wrong)
    @test octo(1,2,3,4,5,6,7,8,true).norm == true # respect the .norm input (even if wrong)
end

@testset "conversions" begin
    @test convert(Quaternion{Float64},1) === Quaternion(1.0)
    @test convert(Quaternion{Float64},Complex(1,2)) === Quaternion(1.0,2.0,0.0,0.0)
    @test convert(Quaternion{Float64},Quaternion(1,2,3,4)) === Quaternion(1.0,2.0,3.0,4.0)

    @test convert(DualQuaternion{Float64},1) === DualQuaternion(1.0)
    @test convert(DualQuaternion{Float64},DualNumbers.Dual(1,2)) === DualQuaternion(Quaternion(1.0),Quaternion(2.0))
    @test convert(DualQuaternion{Float64},Quaternion(1,2,3,4)) === DualQuaternion(Quaternion(1.0,2.0,3.0,4.0))
    @test convert(DualQuaternion{Float64},DualQuaternion(Quaternion(1,2,3,4),Quaternion(5,6,7,8))) === DualQuaternion(Quaternion(1.0,2.0,3.0,4.0),Quaternion(5.0,6.0,7.0,8.0))

    @test convert(Octonion{Float64},1) === Octonion(1.0)
    @test convert(Octonion{Float64},Complex(1,2)) === Octonion(1.0,2.0,0.0,0.0,0.0,0.0,0.0,0.0)
    @test convert(Octonion{Float64},Quaternion(1,2,3,4)) === Octonion(1.0,2.0,3.0,4.0,0.0,0.0,0.0,0.0)
    @test convert(Octonion{Float64},Octonion(1,2,3,4,5,6,7,8)) === Octonion(1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0)
end

@testset "rotations" begin # test rotations
    @test qrotation([0, 0, 0], 1.0) == Quaternion(1.0) # a zero axis should act like zero rotation
    @test qrotation([1, 0, 0], 0.0) == Quaternion(1.0)
    @test qrotation([0, 0, 0]) == Quaternion(1.0)
    qx = qrotation([1, 0, 0], pi / 4)
    @test qx * qx ≈ qrotation([1, 0, 0], pi / 2)
    @test qx^2 ≈ qrotation([1, 0, 0], pi / 2)
    theta = pi / 8
    qx = qrotation([1, 0, 0], theta)
    c = cos(theta); s = sin(theta)
    Rx = [1 0 0; 0 c -s; 0 s c]
    @test rotationmatrix(qx) ≈ Rx
    theta = pi / 6
    qy = qrotation([0, 1, 0], theta)
    c = cos(theta); s = sin(theta)
    Ry = [c 0 s; 0 1 0; -s 0 c]
    @test rotationmatrix(qy) ≈ Ry
    theta = 4pi / 3
    qz = qrotation([0, 0, 1], theta)
    c = cos(theta); s = sin(theta)
    Rz = [c -s 0; s c 0; 0 0 1]
    @test rotationmatrix(qz) ≈ Rz

    @test rotationmatrix(qx * qy * qz) ≈ Rx * Ry * Rz
    @test rotationmatrix(qy * qx * qz) ≈ Ry * Rx * Rz
    @test rotationmatrix(qz * qx * qy) ≈ Rz * Rx * Ry

    a, b = qrotation([0, 0, 1], deg2rad(0)), qrotation([0, 0, 1], deg2rad(180))
    @test slerp(a, b, 0.0) ≈ a
    @test slerp(a, b, 1.0) ≈ b
    @test slerp(a, b, 0.5) ≈ qrotation([0, 0, 1], deg2rad(90))

    @test angle(qrotation([1, 0, 0], 0)) ≈ 0
    @test angle(qrotation([0, 1, 0], pi / 4)) ≈ pi / 4
    @test angle(qrotation([0, 0, 1], pi / 2)) ≈ pi / 2

    let # test numerical stability of angle
        ax = randn(3)
        for θ in [1e-9, π - 1e-9]
            q = qrotation(ax, θ)
            @test angle(q) ≈ θ
        end
    end
end

# Regression test for
# https://github.com/JuliaGeometry/Quaternions.jl/issues/8#issuecomment-610640094
struct MyReal <: Real
    val::Real
end
Base.:(/)(a::MyReal, b::Real) = a.val / b
# this used to throw an error
@test qrotation([1, 0, 0], MyReal(1.5)) == qrotation([1, 0, 0], 1.5)

@testset "non-analytic functions" begin
    q, q2 = randn(Quaternion{Float64}, 2)
    unary_funs = [conj, abs, abs2, norm, sign]
    # since every quaternion is conjugate to a complex number,
    # one can establish correctness as follows:
    @testset for fun in unary_funs
        for i in 1:100
            c = randn(ComplexF64)
            @test fun(Quaternion(c)) ≈ fun(c)
            @test q2 * fun(q) * inv(q2) ≈ fun(q2 * q * inv(q2))
        end
    end
end

@testset "extended complex analytic functions" begin
    # all complex analytic functions can be extended to the quaternions
    unary_funs = [
        sqrt, inv, exp, exp2, exp10, expm1, log, log2, log10, log1p, cis,
        sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh,
        csc, sec, cot, acsc, asec, acot, csch, sech, coth, acsch, asech, acoth,
        sinpi, cospi,
    ]
    # since every quaternion is conjugate to a complex number,
    # one can establish correctness as follows:
    @testset for fun in unary_funs
        q, q2 = randn(QuaternionF64, 2)
        for i in 1:100
            c = randn(ComplexF64)
            fun !== cis && @test fun(Quaternion(c)) ≈ fun(c)
            @test q2 * fun(q) * inv(q2) ≈ fun(q2 * q * inv(q2))
        end
    end

    @testset "identities" begin
        for _ in 1:100
            q = randn(QuaternionF64)
            @test inv(q) * q ≈ q * inv(q) ≈ one(q)
            @test sqrt(q) * sqrt(q) ≈ q
            @test exp(log(q)) ≈ q
            @test exp(zero(q)) ≈ one(q)
            @test log(one(q)) ≈ zero(q)
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
            @testset for (f, finv) in [(sin, csc), (cos, sec), (tan, cot), (sinh, csch), (cosh, sech), (tanh, coth)]
                @test f(q) ≈ inv(finv(q))
            end
            @testset for (f, finv) in [(asin, acsc), (acos, asec), (atan, acot), (asinh, acsch), (acosh, asech), (atanh, acoth)]
                @test f(q) ≈ finv(inv(q))
            end
            @test cis(q) ≈ exp(normalize(q - real(q)) * q)
            VERSION ≥ v"1.6" && @test cispi(q) ≈ cis(π * q)
        end
    end

    @testset "additional properties" begin
        @testset "log" begin
            @test log(zero(QuaternionF64)) === Quaternion(-Inf, 0, 0, 0)
        end

        @testset "exp" begin
            @test exp(Quaternion(0, 0, 0, 0)) == Quaternion(1, 0, 0, 0, true)
            @test exp(Quaternion(2, 0, 0, 0)) == Quaternion(exp(2), 0, 0, 0, false)
            @test exp(Quaternion(0, 2, 0, 0)) == Quaternion(cos(2), sin(2), 0, 0, true)
            @test exp(Quaternion(0, 0, 2, 0)) == Quaternion(cos(2), 0, sin(2), 0, true)
            @test exp(Quaternion(0, 0, 0, 2)) == Quaternion(cos(2), 0, 0, sin(2), true)

            @test norm(exp(Quaternion(0, 0, 0, 0))) ≈ 1
            @test norm(exp(Quaternion(2, 0, 0, 0))) ≠ 1
            @test norm(exp(Quaternion(0, 2, 0, 0))) ≈ 1
            @test norm(exp(Quaternion(0, 0, 2, 0))) ≈ 1
            @test norm(exp(Quaternion(0, 0, 0, 2))) ≈ 1

            @test exp(Quaternion(0., 0., 0., 0.)) == Quaternion(1, 0, 0, 0, true)
            @test exp(Quaternion(2., 0., 0., 0.)) == Quaternion(exp(2), 0, 0, 0, false)
            @test exp(Quaternion(0., 2., 0., 0.)) == Quaternion(cos(2), sin(2), 0, 0, true)
            @test exp(Quaternion(0., 0., 2., 0.)) == Quaternion(cos(2), 0, sin(2), 0, true)
            @test exp(Quaternion(0., 0., 0., 2.)) == Quaternion(cos(2), 0, 0, sin(2), true)

            @test norm(exp(Quaternion(0., 0., 0., 0.))) ≈ 1
            @test norm(exp(Quaternion(2., 0., 0., 0.))) ≠ 1
            @test norm(exp(Quaternion(0., 2., 0., 0.))) ≈ 1
            @test norm(exp(Quaternion(0., 0., 2., 0.))) ≈ 1
            @test norm(exp(Quaternion(0., 0., 0., 2.))) ≈ 1

            @test exp(Quaternion(0,0,0,0)) isa Quaternion{Float64}
            @test exp(Quaternion(0.,0,0,0)) isa Quaternion{Float64}
            @test exp(Quaternion(0//1,0,0,0)) isa Quaternion{Float64}
            @test exp(Quaternion(BigFloat(0),0,0,0)) isa Quaternion{BigFloat}

            # https://github.com/JuliaGeometry/Quaternions.jl/issues/39
            @testset "exp(::Quaternion{Int})" begin
                @test exp(Quaternion(1,1,1,1)) ≈ exp(Quaternion(1.0,1.0,1.0,1.0))
            end
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
        @test Quaternion(ℯ,0,0,0)^Quaternion(0,0,π/2,0) ≈ Quaternion(0,0,1,0)
        @test Quaternion(3.5,0,0,2.3)^Quaternion(0.2,0,0,1.7) ≈
            Quaternion(real((3.5+2.3im)^(0.2+1.7im)),0,0,imag((3.5+2.3im)^(0.2+1.7im)))
        for _ in 1:100
            q, p = randn(QuaternionF64, 2)
            @test q^p ≈ exp(p * log(q))
        end
    end
end

for _ in 1:100
    let # test qrotation and angleaxis inverse
        ax = randn(3); ax = ax / norm(ax)
        Θ = π * rand()
        q = qrotation(ax, Θ)
        @test angle(q) ≈ Θ
        @test axis(q) ≈ ax
        @test angleaxis(q)[1] ≈ Θ
        @test angleaxis(q)[2] ≈ ax
    end

    let # test argq
        q, q2 = sample(Quaternion{Float64}, 2)
        @test q2 * argq(q) * inv(q2) ≈ argq(q2 * q * inv(q2))
        v = Quaternion(0, randn(3)...)
        @test argq(v) * norm(v) ≈ v
    end

    let # test normalize
        q = quatrand()
        @test norm(normalize(q)) ≈ 1
        @test normalize(q).norm
        @test q ≈ norm(q) * normalize(q)
        qn = nquatrand()
        @test qn.norm
        @test normalize(qn) === qn
    end

    let # test rotation <=> rotationmatrix
        q1 = nquatrand()
        q2 = qrotation(rotationmatrix(q1), q1)
        @test q1 ≈ q2
    end

    let # test slerp and linpol if q1 = 1
        q1 = quat(1, 0, 0, 0.)
        # there are numerical stability issues with slerp atm
        θ = clamp(rand() * 3.5, deg2rad(5e-1), π)
        ax = randn(3)
        q2 = qrotation(ax, θ)
        t = rand()
        slerp(q1, q2, 0.) ≈ q1
        @test slerp(q1, q2, 0.) ≈ q1
        @test slerp(q1, q2, 1.) ≈ q2
        @test slerp(q1, q2, t) ≈ qrotation(ax, t * θ)
        @test norm(slerp(q1, q2, t)) ≈ 1
        @test slerp(q1, q2, 0.5) ≈ qrotation(ax, 0.5 * θ)
        @test linpol(q1, q2, 0.5) ≈ qrotation(ax, 0.5 * θ)

    end
    let # test conjugation invariance
        q, q1, q2 = sample(Quaternion{Float64}, 3)
        ⊗(s, t) = s * t * inv(s)
        t = rand()
        @test q ⊗ slerp(q1, q2, t) ≈ slerp(q ⊗ q1, q ⊗ q2, t)
        @test q ⊗ linpol(q1, q2, t) ≈ linpol(q ⊗ q1, q ⊗ q2, t)
    end
end

@testset "random quaternions" begin
    @testset "quatrand" begin
        rng = Random.MersenneTwister(42)
        q1 = quatrand(rng)
        @test q1 isa Quaternion
        @test !q1.norm

        q2 = quatrand()
        @test q2 isa Quaternion
        @test !q2.norm
    end

    @testset "nquatrand" begin
        rng = Random.MersenneTwister(42)
        q1 = nquatrand(rng)
        @test q1 isa Quaternion
        @test q1.norm

        q2 = nquatrand()
        @test q2 isa Quaternion
        @test q2.norm
    end

    @testset "rand($H)" for H in (Quaternion, DualQuaternion, Octonion)
        rng = Random.MersenneTwister(42)
        q1 = rand(rng, H{Float64})
        @test q1 isa H{Float64}
        @test !q1.norm

        q2 = rand(rng, H{Float32})
        @test q2 isa H{Float32}
        @test !q2.norm

        qs = rand(rng, H{Float64}, 1000)
        @test eltype(qs) === H{Float64}
        @test length(qs) == 1000
        xs = map(qs) do q
            if q isa DualQuaternion
                return [real(q.q0); Quaternions.imag(q.q0); real(q.qe); Quaternions.imag(q.qe)]
            else
                return [real(q); Quaternions.imag(q)]
            end
        end
        xs_mean = sum(xs) / length(xs)
        xs_var = sum(x -> abs2.(x .- xs_mean), xs) / (length(xs) - 1)
        @test all(isapprox.(xs_mean, 0.5; atol=0.1))
        @test all(isapprox.(xs_var, 1/12; atol=0.01))
    end

    @testset "randn($H)" for H in (Quaternion, Octonion)
        rng = Random.MersenneTwister(42)
        q1 = randn(rng, H{Float64})
        @test q1 isa H{Float64}
        @test !q1.norm

        q2 = randn(rng, H{Float32})
        @test q2 isa H{Float32}
        @test !q2.norm

        qs = randn(rng, H{Float64}, 10000)
        @test eltype(qs) === H{Float64}
        @test length(qs) == 10000
        xs = map(qs) do q
            if q isa DualQuaternion
                return [real(q.q0); Quaternions.imag(q.q0); real(q.qe); Quaternions.imag(q.qe)]
            else
                return [real(q); Quaternions.imag(q)]
            end
        end
        xs_mean = sum(xs) / length(xs)
        xs_var = sum(x -> abs2.(x .- xs_mean), xs) / (length(xs) - 1)
        @test all(isapprox.(xs_mean, 0; atol=0.1))
        if H === Quaternion
            @test all(isapprox.(xs_var, 1/4; atol=0.1))
        else
            @test all(isapprox.(xs_var, 1/8; atol=0.1))
        end
    end
end
