# Quaternions.jl
A Julia implementation of quaternions.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGeometry.github.io/Quaternions.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGeometry.github.io/Quaternions.jl/dev)
[![Build Status](https://github.com/JuliaGeometry/Quaternions.jl/workflows/CI/badge.svg)](https://github.com/JuliaGeometry/Quaternions.jl/actions?query=workflow%3ACI+branch%3Amain)
[![codecov](https://codecov.io/gh/JuliaGeometry/Quaternions.jl/branch/main/graph/badge.svg?token=dJBiR91dCD)](https://codecov.io/gh/JuliaGeometry/Quaternions.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

[Quaternions](http://en.wikipedia.org/wiki/Quaternion) are best known for their suitability
as representations of 3D rotational orientation.
They can also be viewed as an extension of complex numbers.

## First example

```julia
julia> using Quaternions

julia> k = quat(0, 0, 0, 1)
Quaternion{Int64}(0, 0, 0, 1)

julia> j = quat(0, 0, 1, 0)
Quaternion{Int64}(0, 0, 1, 0)

julia> i = j*k
Quaternion{Int64}(0, 1, 0, 0)

julia> i^2 == j^2 == k^2 == i*j*k == -1 # Similar to `im^2`.
true

julia> 1 + i + k + j  # Compatible with arithmetic operations as a `Number`.
Quaternion{Int64}(1, 1, 1, 1)
```

Check out [the docs](https://juliageometry.github.io/Quaternions.jl) for further instructions.

In JuliaGeometry organization, there is also [Octonions.jl](https://github.com/JuliaGeometry/Octonions.jl) package.
