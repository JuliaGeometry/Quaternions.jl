# Basics

## Basic operations for quaternions
Quaternions can be defined with [`Quaternion`](@ref) constructor or [`quat`](@ref) function.
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

A conjugate of a quaternion can be calculated with [`Base.conj`](@ref).
But `Base.imag(::Quaternion)` is not defined because it should return three real values which is not consistent with `imag(::Complex)` and `imag(::Real)`.
Instead, [`imag_part`](@ref) function can be used to obtain the imaginary part of a quaternion.
```@repl intro
conj(q1)
imag_part(q1)
imag(q1)
```

## `Quaternion` vs `quat`
The general rule is that `quat` is to `Quaternion` as `complex` is to `Complex`.
`Complex` and `Quaternion` are both constructors so should return an object of the corresponding type, whereas `quat` and `complex` both can operate on types and arrays.

```@setup Quaternion-quat
using Quaternions
```

```@repl Quaternion-quat
Quaternion(1,2,3,4)
quat(1,2,3,4)
Quaternion(Int)
quat(Int)
```

## The type parameter in `Quaternion{T}`


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
