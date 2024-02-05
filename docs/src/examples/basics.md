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

Unit quaternions can be obtained with [`sign`](@ref).

```@repl intro
sign(q1)
sign(q2)
sign(q3)
sign(quat(0))  # Zero-quaternion will not be normalized.
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

## Compatibility with `Complex`
There is no natural embedding ``\mathbb{C}\to\mathbb{H}``.
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
