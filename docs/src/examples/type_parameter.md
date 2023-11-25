# The type parameter `T` in `Quaternion{T}`

The type parameter `T <: Real` in `Quaternion{T}` represents the type of real and imaginary parts of a quaternion.

## Lipschitz quaternions
By using this type parameter, some special quaternions such as [**Lipschitz quaternions**](https://en.wikipedia.org/wiki/Hurwitz_quaternion) ``L`` can be represented.

```math
L = \left\{a+bi+cj+dk \in \mathbb{H} \mid a,b,c,d \in \mathbb{Z}\right\}
```

```@setup LipschitzHurwitz
using Quaternions
```

```@repl LipschitzHurwitz
q1 = Quaternion{Int}(1,2,3,4)
q2 = Quaternion{Int}(5,6,7,8)
islipschitz(q::Quaternion) = isinteger(q.s) & isinteger(q.v1) & isinteger(q.v2) & isinteger(q.v3)
islipschitz(q1)
islipschitz(q2)
islipschitz(q1 + q2)
islipschitz(q1 * q2)
islipschitz(q1 / q2)  # Division is not defined on L.
q1 * q2 == q2 * q1  # non-commutative
```

## Hurwitz quaternions
If all coefficients of a quaternion are integers or half-integers, the quaternion is called a [**Hurwitz quaternion**](https://en.wikipedia.org/wiki/Hurwitz_quaternion).
The set of Hurwitz quaternions is defined by

```math
H = \left\{a+bi+cj+dk \in \mathbb{H} \mid a,b,c,d \in \mathbb{Z} \ \text{or} \ a,b,c,d \in \mathbb{Z} + \tfrac{1}{2}\right\}.
```

Hurwitz quaternions can be implemented with [HalfIntegers.jl](https://github.com/sostock/HalfIntegers.jl) package.

```@repl LipschitzHurwitz
using HalfIntegers
q1 = Quaternion{HalfInt}(1, 2, 3, 4)
q2 = Quaternion{HalfInt}(5.5, 6.5, 7.5, 8.5)
q3 = Quaternion{HalfInt}(1, 2, 3, 4.5)  # not Hurwitz quaternion
ishurwitz(q::Quaternion) = (isinteger(q.s) & isinteger(q.v1) & isinteger(q.v2) & isinteger(q.v3)) | (ishalfinteger(q.s) & ishalfinteger(q.v1) & ishalfinteger(q.v2) & ishalfinteger(q.v3))
ishurwitz(q1)
ishurwitz(q2)
ishurwitz(q3)
ishurwitz(q1 + q2)
ishurwitz(q1 * q2)
ishurwitz(q1 / q2)  # Division is not defined on H.
q1 * q2 == q2 * q1  # non-commucative
abs2(q1)  # Squared norm is always an integer.
abs2(q2)  # Squared norm is always an integer.
abs2(q3)  # Squared norm is not an integer because `q3` is not Hurwitz quaternion.
```

## Biquaternions
If all coefficients of a quaternion are complex numbers, the quaternion is called a [**Biquaternion**](https://en.wikipedia.org/wiki/Biquaternion).
However, the type parameter `T` is restricted to `<:Real`, so biquaternions are not supported in this package.
Note that `Base.Complex` has the same type parameter restriction, and [bicomplex numbers](https://en.wikipedia.org/wiki/Bicomplex_number) are not supported in Base.
See [issue#79](https://github.com/JuliaGeometry/Quaternions.jl/issues/79) for more discussion.
