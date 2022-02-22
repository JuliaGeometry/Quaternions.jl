# Quaternions.jl

A Julia module with [quaternion](https://en.wikipedia.org/wiki/Quaternion), [octonion](https://en.wikipedia.org/wiki/Octonion) and [dual-quaternion](https://en.wikipedia.org/wiki/Dual_quaternion) functionality

!!! note "Documentation"
    The documentation is still work in progress.
    For more information, see also
    * [README in the repository](https://github.com/JuliaGeometry/Quaternions.jl)
    * [Tests in the repository](https://github.com/JuliaGeometry/Quaternions.jl/tree/master/test)
    Feel free to [open pull requests](https://github.com/JuliaGeometry/Quaternions.jl/pulls) and improve this document!

## Installation
```
pkg> add Quaternions
```

## First example

```@repl
using Quaternions
q = quat(0.0, 0.0, 0.0, 1.0)
r = quat(0, 0, 1, 0)
q*r
q+r
```
