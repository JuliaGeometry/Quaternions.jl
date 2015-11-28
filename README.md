# Quaternions.jl
A Julia module with quaternion, octonion and dual-quaternion functionality

[![Build Status](https://travis-ci.org/JuliaGeometry/Quaternions.jl.png?branch=master)](https://travis-ci.org/JuliaGeometry/Quaternions.jl)
[![Coverage Status](https://coveralls.io/repos/JuliaGeometry/Quaternions.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaGeometry/Quaternions.jl?branch=master)

[Quaternions](http://en.wikipedia.org/wiki/Quaternion) are best known for their suitability
as representations of 3D rotational orientation. They can also be viewed as an extension of complex numbers.

Implemented functions are:  

    +-*/^
    real  
    imag  (a vector)  
    conj  
    abs  
    abs2  
    exp  
    log  
    normalize  
    normalizea  (return normalized quaternion and absolute value as a pair)  
    angleaxis  (taken as an orientation, return the angle and axis (3 vector) as a tuple)  
    angle  
    axis  
    exp  
    log  
    sin  
    cos  
    sqrt  
    linpol  (interpolate between 2 normalized quaternions)  

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
