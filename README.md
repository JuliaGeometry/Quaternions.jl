# Quaternions.jl
A Julia module with quaternion, octonion and dual-quaternion functionality

[![Build Status](https://github.com/JuliaGeometry/Quaternions.jl/workflows/CI/badge.svg)](https://github.com/JuliaGeometry/Quaternions.jl/actions?query=workflow%3ACI+branch%3Amaster)
[![codecov](https://codecov.io/gh/JuliaGeometry/Quaternions.jl/branch/master/graph/badge.svg?token=dJBiR91dCD)](https://codecov.io/gh/JuliaGeometry/Quaternions.jl)

[Quaternions](http://en.wikipedia.org/wiki/Quaternion) are best known for their suitability
as representations of 3D rotational orientation. They can also be viewed as an extension of complex numbers.

Implemented functions are:

    +-*/^
    real
    imag  (a vector)
    conj
    abs
    abs2
    normalize
    normalizea  (return normalized quaternion and absolute value as a pair)
    angleaxis  (taken as an orientation, return the angle and axis (3 vector) as a tuple)
    angle
    axis
    sqrt
    exp
    exp2
    exp10
    expm1
    log2
    log10
    log1p
    cis
    cispi
    sin
    cos
    tan
    asin
    acos
    atan
    sinh
    cosh
    tanh
    asinh
    acosh
    atanh
    csc
    sec
    cot
    acsc
    asec
    acot
    csch
    sech
    coth
    acsch
    asech
    acoth
    sinpi
    cospi
    sincos
    sincospi
    linpol  (interpolate between 2 normalized quaternions)
    rand
    randn

[Dual quaternions](http://en.wikipedia.org/wiki/Dual_quaternion) are an extension, combining quaternions with
[dual numbers](https://github.com/scidom/DualNumbers.jl).
On top of just orientation, they can represent all rigid transformations.

There are two conjugation concepts here

    conj  (quaternion conjugation)
    dconj (dual conjugation)

further implemented here:

    Q0  (the 'real' quaternion)
    Qe  ( the 'dual' part)
    +-*/^
    abs
    abs2
    normalize
    normalizea
    angleaxis
    angle
    axis
    exp
    log
    sqrt
    rand

[Octonions](http://en.wikipedia.org/wiki/Octonion) form the logical next step on the Complex-Quaternion path.
They play a role, for instance, in the mathematical foundation of String theory.

    +-*/^
    real
    imag  (a vector)
    conj
    abs
    abs2
    exp
    log
    normalize
    normalizea  (return normalized octonion and absolute value as a tuple)
    exp
    log
    sqrt
    rand
    randn
