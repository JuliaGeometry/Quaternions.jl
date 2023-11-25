# Quaternions.jl

A Julia package implementing [quaternions](https://en.wikipedia.org/wiki/Quaternion).

!!! note "Documentation"
    The documentation is still work in progress.
    For more information, see also
    * [README in the repository](https://github.com/JuliaGeometry/Quaternions.jl)
    * [Tests in the repository](https://github.com/JuliaGeometry/Quaternions.jl/tree/main/test)
    Feel free to [open pull requests](https://github.com/JuliaGeometry/Quaternions.jl/pulls) and improve this document!

## Installation
```
pkg> add Quaternions
```

## First example

```@repl
using Quaternions
k = quat(0, 0, 0, 1)
j = quat(0, 0, 1, 0)
i = j*k
i^2 == j^2 == k^2 == i*j*k == -1 # Similar to `im^2`.
1 + i + k + j  # Compatible with arithmetic operations as a `Number`.
```
