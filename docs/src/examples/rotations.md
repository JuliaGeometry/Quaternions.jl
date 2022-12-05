# Rotations with quaternions

One of the most useful application of quaternions is representation of 3D-rotations.
See also [Rotations.jl documentation](https://juliageometry.github.io/Rotations.jl/stable/3d_quaternion/)

```@example rotation
using Quaternions
using LinearAlgebra
```

## Basics
A 3D rotation can be represented by a [unit quaternion (versor)](https://en.wikipedia.org/wiki/Versor).
For example, a 90° rotation around the ``y``-axis is ``q = \frac{1}{\sqrt{2}} + 0i + \frac{1}{\sqrt{2}} j + 0k``.
Rotations with quaternions have the following properties:

* A unit quaternion (4 real numbers) is more efficient for representing a rotation than a rotation matrix (9 real numbers).
    * This results in higher computational performance in terms of time, memory usage, and accuracy.
* The negative of a unit quaternion represents the same rotation.
* The conjugate of a unit quaternion represents the inverse rotation.
    * The quaternion has unit length, so conjugate and multiplicative inverse is the same.
* The set of unit quaternion ``\left\{w + ix + jy + kz \in \mathbb{H} \ | \ x, y, z \in \mathbb{R} \right\} = U(1,\mathbb{H}) \simeq S^3`` is isomorphic to ``SU(2)``.
    * These groups are homomorphic to ``SO(3)``.
    * These groups are double covering of ``SO(3)``.

## Rotation around a vector
A ``\theta`` rotation around a unit vector ``v = (v_x, v_y, v_z)`` can be obtained as
```math
q = \cos(\theta/2) + \sin(\theta/2)(iv_x + jv_y + kv_z).
```

```@example rotation
function quat_from_axisangle(axis, theta)
    if length(axis) != 3
        error("Must be a 3-vector")
    end
    s, c = sincos(theta / 2)
    axis = normalize(axis)
    return Quaternion(c, s*axis[1], s*axis[2], s*axis[3])
end
nothing  # hide
```

```@repl rotation
q1 = quat_from_axisangle([0,1,0], deg2rad(90))  # 90° rotation around y-axis
q2 = quat_from_axisangle([1,1,1], deg2rad(120))
q3 = -q2  # additive inverse quaternion represents the same rotation
```

## Rotate a vector with a quaternion
A vector ``v = (v_x, v_y, v_z)`` can be rotated by a unit quaternion ``q``.
The rotated vector ``v' = (v_x', v_y', v_z')`` can be obtained as
```math
\begin{aligned}
q_v &= iv_x + jv_y + kv_z \\
q_v' &= q q_v \bar{q} = 0 + iv_x + jv_y + kv_z \\
v' &= (v_x', v_y', v_z').
\end{aligned}
```

```@example rotation
function rotate_vector(q::Quaternion, vector)
    if length(vector) != 3
        error("Must be a 3-vector")
    end
    q_v = Quaternion(0, vector[1], vector[2], vector[3])
    q_v′ = q*q_v*conj(q)
    return [imag_part(q_v′)...]
end
nothing  # hide
```

```@repl rotation
rotate_vector(q1, [1,2,3])
rotate_vector(q2, [1,2,3])
rotate_vector(q3, [1,2,3])  # Same as q2
```

## Convert a quaternion to a rotation matrix
A unit quaternion can be converted to a rotation matrix.

```@example rotation
function rotmatrix_from_quat(q::Quaternion)
    sx, sy, sz = 2q.s * q.v1, 2q.s * q.v2, 2q.s * q.v3
    xx, xy, xz = 2q.v1^2, 2q.v1 * q.v2, 2q.v1 * q.v3
    yy, yz, zz = 2q.v2^2, 2q.v2 * q.v3, 2q.v3^2
    r = [1 - (yy + zz)     xy - sz     xz + sy;
            xy + sz   1 - (xx + zz)    yz - sx;
            xz - sy      yz + sx  1 - (xx + yy)]
    return r
end
nothing  # hide
```

```@repl rotation
m1 = rotmatrix_from_quat(q1)
m2 = rotmatrix_from_quat(q2)
m3 = rotmatrix_from_quat(q3)  # Same as q2
```

This function does not return [`StaticMatrix`](https://juliaarrays.github.io/StaticArrays.jl/dev/pages/api/#StaticArraysCore.StaticArray), so the implementation is not much effective.
If you need more performance, please consider using [Rotations.jl](https://github.com/JuliaGeometry/Rotations.jl).

## Convert a rotation matrix to a quaternion
A rotation matrix can be converted to a unit quaternion.
The following implementation is based on [https://arxiv.org/pdf/math/0701759.pdf](https://arxiv.org/pdf/math/0701759.pdf).
Note that the following mapping ``SO(3) \to SU(2)`` is not surjective.

```@example rotation
function quat_from_rotmatrix(dcm::AbstractMatrix{T}) where {T<:Real}
    a2 = 1 + dcm[1,1] + dcm[2,2] + dcm[3,3]
    a = sqrt(a2)/2
    b,c,d = (dcm[3,2]-dcm[2,3])/4a, (dcm[1,3]-dcm[3,1])/4a, (dcm[2,1]-dcm[1,2])/4a
    return Quaternion(a,b,c,d)
end
nothing  # hide
```

```@repl rotation
quat_from_rotmatrix(m1)
quat_from_rotmatrix(m2)
quat_from_rotmatrix(m3)
quat_from_rotmatrix(m1) ≈ q1
quat_from_rotmatrix(m2) ≈ q2
quat_from_rotmatrix(m3) ≈ q3  # q2 == -q3
```

## Interpolate two rotations (slerp)
Slerp (spherical linear interpolation) is a method to interpolate between two unit quaternions.
This function [`slerp`](@ref) equates antipodal points, and interpolates the shortest path.
Therefore, the output `slerp(q1, q2, 1)` may be different from `q2`. (`slerp(q1, q2, 0)` is always equal to `q1`.)

```@repl rotation
slerp(q1, q2, 0) ≈ q1
slerp(q1, q2, 1) ≈ q2
slerp(q1, q3, 1) ≈ q3
slerp(q1, q3, 1) ≈ -q3
r = slerp(q1, q2, 1/2)
abs(q1-r) ≈ abs(q2-r)  # Same distance
abs(r)  # Interpolates on the unit sphere S³
```
