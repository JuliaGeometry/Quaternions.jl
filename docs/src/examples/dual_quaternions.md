# Dual quaternions

## Introduction

The [dual quaternions](https://en.wikipedia.org/wiki/Dual_quaternion) are an example of "biquaternions."
They can be represented equivalently either as a [dual number](https://en.wikipedia.org/wiki/Dual_number) where both both the "primal" and "tangent" part are quaternions

```math
d = q_0 + q_e \epsilon = (s_0 + a_0 i + b_0 j + c_0 k) + (s_e + a_e i + b_e j + c_e k) \epsilon
```

or as a quaternion where the scalar part and three imaginary parts are all dual numbers

```math
d = s + ai + bj + ck = (s_0 + s_e \epsilon) + (a_0 + a_e \epsilon) i + (b_0 + b_e \epsilon) j + (c_0 + c_e \epsilon) k.
```

Like unit quaternions can compactly representation rotations in 3D space, dual quaternions can compactly represent rigid transformations (rotation with translation).

Without any special glue code, we can construct a dual quaternion by composing `ForwardDiff.Dual` and [`Quaternion`](@ref); this uses the second representation described above:

!!! note
    Previously this package contained a specialized `DualQuaternion` type.
    This was removed in v0.6.0 because offered nothing extra over composing [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) and Quaternions.

## Utility functions

First let's load the packages:

```@example dualquat
using Quaternions, ForwardDiff, Random
```

Then we'll create some utility types/functions:

```@example dualquat
const DualQuaternion{T} = Quaternion{ForwardDiff.Dual{Nothing,T,1}}

purequat(p::AbstractVector) = quat(false, @views(p[begin:begin+2])...)

dual(x::Real, v::Real) = ForwardDiff.Dual(x, v)

function dualquat(_q0::Union{Real,Quaternion}, _qe::Union{Real,Quaternion})
    q0 = quat(_q0)
    qe = quat(_qe)
    Quaternion(
        dual(real(q0), real(qe)),
        dual.(imag_part(q0), imag_part(qe))...,
    )
end

function primal(d::DualQuaternion)
    return Quaternion(
        ForwardDiff.value(real(d)),
        ForwardDiff.value.(imag_part(d))...,
    )
end

function tangent(d::DualQuaternion)
    return Quaternion(
        ForwardDiff.partials(real(d), 1),
        ForwardDiff.partials.(imag_part(d), 1)...,
    )
end

function dualconj(d::DualQuaternion)
    de = tangent(d)
    return dualquat(conj(primal(d)), quat(-real(de), imag_part(de)...))
end

rotation_part(d::DualQuaternion) = primal(d)

translation_part(d::DualQuaternion) = dualquat(true, conj(rotation_part(d)) * tangent(d))

# first=true returns the translation performed before the rotation: R(p+t)
# first=false returns the translation performed after the rotation: R(p)+t
function translation(d::DualQuaternion; first::Bool=true)
    v = first ? primal(d)' * tangent(d) : tangent(d) * primal(d)'
    return collect(2 .* imag_part(v))
end

function transform(d::DualQuaternion, p::AbstractVector)
    dp = dualquat(true, purequat(p))
    dpnew = d * dp * dualconj(d)
    pnew_parts = imag_part(tangent(dpnew))
    pnew = similar(p, eltype(pnew_parts))
    pnew .= pnew_parts
    return pnew
end

function transformationmatrix(d::DualQuaternion)
    R = rotationmatrix(rotation_part(d))
    t = translation(d; first=false)
    T = similar(R, 4, 4)
    T[1:3, 1:3] .= R
    T[1:3, 4] .= t
    T[4, 1:3] .= 0
    T[4, 4] = 1
    return T
end

randdualquat(rng::AbstractRNG,T=Float64) = dualquat(rand(rng, Quaternion{T}), rand(rng, Quaternion{T}))
randdualquat(T=Float64) = randdualquat(Random.GLOBAL_RNG,T)
nothing  # hide
```

## Example: transforming a point

Now we'll create a dual quaternion.
```@repl dualquat
x = sign(randdualquat())
```

`sign(q) == q / abs(q)` both normalizes the primal part of the dual quaternion and makes the tangent part perpendicular to it.

```@repl dualquat
abs(primal(x)) ≈ 1
isapprox(real(primal(x)' * tangent(x)), 0; atol=1e-10)
```

Here's how we use dual quaternions to transform a point:

```@repl dualquat
p = randn(3)
```

```@repl dualquat
transform(x, p)
```

We can check that this gives the same result as transforming with an affine transformation matrix:

```@example dualquat
transform(x, p) ≈ (transformationmatrix(x) * vcat(p, 1))[1:3]
```

## Example: homomorphism from dual quaternions to matrices

```@repl dualquat
y = sign(randdualquat())
```

```@repl dualquat
X = transformationmatrix(x)
Y = transformationmatrix(y)
XY = transformationmatrix(x*y)
X*Y ≈ XY
```
## Example: motion planning

For unit quaternions, spherical linear interpolation with [`slerp`](@ref) can be used to interpolate between two rotations with unit quaternions, which can be used to plan motion between two orientations.
Similarly, we can interpolate between unit dual quaternions to plan motion between two rigid poses.
Conveniently, we can do this using the exact same `slerp` implementation.

```@example dualquat
y = sign(randdualquat())
```

```@example dualquat
slerp(x, y, 0) ≈ x
```

```@example dualquat
slerp(x, y, 1) ≈ y
```

```@example dualquat
slerp(x, y, 0.3)
```
