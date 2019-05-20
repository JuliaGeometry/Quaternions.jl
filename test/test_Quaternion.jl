
using Quaternions: argq
using LinearAlgebra

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

let # test rotations
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

for _ in 1:100
    let # test specialfunctions
        c = Complex(randn(2)...)
        q, q2 = sample(Quaternion{Float64}, 4)
        unary_funs = [exp, log, sin, cos, sqrt, inv, conj, abs2, norm]
        # since every quaternion is conjugate to a complex number,
        # one can establish correctness as follows:
        for fun in unary_funs
            @test fun(Quaternion(c)) ≈ Quaternion(fun(c))
            @test q2 * fun(q) * inv(q2) ≈ fun(q2 * q * inv(q2))
        end

        @test exp(log(q)) ≈ q
        @test exp(zero(q)) ≈ one(q)
    end

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
