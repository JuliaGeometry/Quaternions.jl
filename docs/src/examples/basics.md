# Basics

## Basic operations for quaternions
Quaternions can be defined with the [`Quaternion`](@ref) constructor or [`quat`](@ref) function.
Note that the order of the arguments is ``w+xi+yj+zk``, not ``xi+yj+zk+w``.

```@repl intro
using Quaternions
q1 = Quaternion(1,2,3,4)
q2 = quat(5,6,7,8)
q3 = quat(9)
```

The multiplication is not commutative.
```@repl intro
q1 * q2
q2 * q1
```

The multiplicative inverse can be calculated with [`Base.inv`](@ref).
```@repl intro
inv(q1)
inv(q1) * q1
```

The division is also not commutative.

```@repl intro
q1 / q2  # Same as `q1*inv(q2)` mathematically.
q2 \ q1  # Same as `inv(q2)*q1` mathematically.
```

A conjugate of a quaternion can be calculated with [`Base.conj`](@ref).
But `Base.imag(::Quaternion)` is not defined because it should return three real values which is not consistent with `imag(::Complex)` and `imag(::Real)`.
Instead, the [`imag_part`](@ref) function can be used to obtain the imaginary part of a quaternion.
See [issue#61](https://github.com/JuliaGeometry/Quaternions.jl/issues/61) for more discussion.

```@repl intro
conj(q1)
imag(q1)  # Not supported.
imag_part(q1)  # Use this instead.
```

## `Quaternion` vs `quat`
The general rule is that [`quat`](@ref) is to [`Quaternion`](@ref) as [`complex`](https://docs.julialang.org/en/v1/base/numbers/#Base.complex-Tuple{Complex}) is to [`Complex`](https://docs.julialang.org/en/v1/base/numbers/#Base.Complex).
`Complex` and `Quaternion` are both constructors so should return an object of the corresponding type, whereas `quat` and `complex` both can operate on types and arrays.

```@setup Quaternion-quat
using Quaternions
```

```@repl Quaternion-quat
Quaternion(1,2,3,4)
quat(1,2,3,4)
Quaternion(Int)  # Similar to `Complex(Int)`.
quat(Int)  # Similar to `complex(Int)`.
```

## The type parameter `T` in `Quaternion{T}`

The type parameter `T <: Real` in `Quaternion{T}` represents the type of real and imaginary parts of a quaternion.

### Lipschitz quaternions
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

### Hurwitz quaternions
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

### Biquaternions
If all coefficients of a quaternion are complex numbers, the quaternion is called a [**Biquaternion**](https://en.wikipedia.org/wiki/Biquaternion).
However, the type parameter `T` is restricted to `<:Real`, so biquaternions are not supported in this package.
Note that `Base.Complex` has the same type parameter restriction, and [bicomplex numbers](https://en.wikipedia.org/wiki/Bicomplex_number) are not supported in Base.
See [issue#79](https://github.com/JuliaGeometry/Quaternions.jl/issues/79) for more discussion.

## Compatibility with `Complex`
There are no natural embedding ``\mathbb{C}\to\mathbb{H}``.
Thus, `quat(w,x,0,0)` is not equal to `complex(w,x)`, i.e.

```math
\mathbb{C} \ni w+ix \ne w+ix+0j+0k \in \mathbb{H}.
```

```@setup complex
using Quaternions
```

```@repl complex
1 + complex(1,2)   # `Complex` is compatible with `Real`
1 + quat(1,2,3,4)  # `Quaternion` is compatible with `Real`
1 + complex(1,2) + quat(1,2,3,4)  # no compatibility
complex(1,2) == quat(1,2,0,0)     # no compatibility
complex(1) == quat(1)             # no compatibility
complex(1) == 1 == quat(1)  # Both `quat(1)` and `complex(1)` are equal to `1`.
```

See [issue#62](https://github.com/JuliaGeometry/Quaternions.jl/issues/62) for more discussion.
