using Base.Test
using Quaternions

let # test rotations
    qx = qrotation([1,0,0], pi/4)
    @test_approx_eq qx*qx qrotation([1,0,0], pi/2)
    @test_approx_eq qx^2 qrotation([1,0,0], pi/2)
    theta = pi/8
    qx = qrotation([1,0,0], theta)
    c = cos(theta); s = sin(theta)
    Rx = [1 0 0; 0 c -s; 0 s c]
    @test_approx_eq rotationmatrix(qx) Rx
    theta = pi/6
    qy = qrotation([0,1,0], theta)
    c = cos(theta); s = sin(theta)
    Ry = [c 0 s; 0 1 0; -s 0 c]
    @test_approx_eq rotationmatrix(qy) Ry
    theta = 4pi/3
    qz = qrotation([0,0,1], theta)
    c = cos(theta); s = sin(theta)
    Rz = [c -s 0; s c 0; 0 0 1]
    @test_approx_eq rotationmatrix(qz) Rz

    @test_approx_eq rotationmatrix(qx*qy*qz) Rx*Ry*Rz
    @test_approx_eq rotationmatrix(qy*qx*qz) Ry*Rx*Rz
    @test_approx_eq rotationmatrix(qz*qx*qy) Rz*Rx*Ry

    a, b = qrotation([0,0,1], deg2rad(0)), qrotation([0,0,1], deg2rad(180))
    @test_approx_eq slerp(a,b,0.0) a
    @test_approx_eq slerp(a,b,1.0) b
    @test_approx_eq slerp(a,b,0.5) qrotation([0,0,1], deg2rad(90))

    @test_approx_eq angle(qrotation([1,0,0], 0)) 0
    @test_approx_eq angle(qrotation([0,1,0], pi/4)) pi/4
    @test_approx_eq angle(qrotation([0,0,1], pi/2)) pi/2

    # test qrotation and angleaxis inverse
    for _ in 1:100
        ax = randn(3); ax = ax/norm(ax)
        Θ = π * rand()
        q = qrotation(ax, Θ)
        @test_approx_eq angle(q) Θ
        @test_approx_eq axis(q) ax
    end

    let # test numerical stability of angle
        ax = randn(3)
        for θ in [1e-9, π - 1e-9]
            q = qrotation(ax, θ)
            @test_approx_eq angle(q) θ
        end
    end
end

for _ in 1:100 # test slerp
    # correctness can be proven if we show
    # 1) slerp is correct for qa = 1.
    # 2) slerp is conjugation invariant.

    let # test slerp if q1 = 1
        q1 = Quaternion(1.,0,0,0)
        θ_max = rand() * π
        ax = randn(3)
        q2 = qrotation(ax, θ_max)
        t = rand()
        @test slerp(q1, q2, 0.) ≈ q1
        @test slerp(q1, q2, 1.) ≈ q2
        @test norm(slerp(q1, q2, t)) ≈ 1
        @test slerp(q1, q2, t) ≈ qrotation(ax, t*θ_max)
    end

    let # test slerp invariant under conjugation action
        q, q1, q2 = [qrotation(randn(3), rand() * π) for _ in 1:3]
        ⊗(s, t) = s*t*inv(s)
        t = rand()
        @test q ⊗ slerp(q1, q2, t) ≈ slerp(q ⊗ q1, q ⊗ q2, t)
    end
end
