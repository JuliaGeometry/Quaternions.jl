struct Quaternion{T<:Real} <: Number
    s::T
    v1::T
    v2::T
    v3::T
    norm::Bool
end

const QuaternionF16 = Quaternion{Float16}
const QuaternionF32 = Quaternion{Float32}
const QuaternionF64 = Quaternion{Float64}

Quaternion{T}(x::Real) where {T<:Real} = Quaternion(convert(T, x))
Quaternion{T}(x::Complex) where {T<:Real} = Quaternion(convert(Complex{T}, x))
Quaternion{T}(q::Quaternion) where {T<:Real} = Quaternion{T}(q.s, q.v1, q.v2, q.v3, q.norm)
Quaternion(s::Real, v1::Real, v2::Real, v3::Real, n::Bool = false) =
    Quaternion(promote(s, v1, v2, v3)..., n)
Quaternion(x::Real) = Quaternion(x, zero(x), zero(x), zero(x), abs(x) == one(x))
Quaternion(z::Complex) = Quaternion(z.re, z.im, zero(z.re), zero(z.re), abs(z) == one(z.re))
Quaternion(s::Real, a::Vector) = Quaternion(s, a[1], a[2], a[3])
Quaternion(a::Vector) = Quaternion(0, a[1], a[2], a[3])

promote_rule(::Type{Quaternion{T}}, ::Type{S}) where {T <: Real, S <: Real} = Quaternion{promote_type(T, S)}
promote_rule(::Type{Quaternion{T}}, ::Type{Complex{S}}) where {T <: Real, S <: Real} = Quaternion{promote_type(T, S)}
promote_rule(::Type{Quaternion{T}}, ::Type{Quaternion{S}}) where {T <: Real, S <: Real} = Quaternion{promote_type(T, S)}

quat(p, v1, v2, v3) = Quaternion(p, v1, v2, v3)
quat(p, v1, v2, v3, n) = Quaternion(p, v1, v2, v3, n)
quat(x) = Quaternion(x)
quat(s, a) = Quaternion(s, a)

function show(io::IO, q::Quaternion)
    pm(x) = x < 0 ? " - $(-x)" : " + $x"
    print(io, q.s, pm(q.v1), "im", pm(q.v2), "jm", pm(q.v3), "km")
end

real(::Type{Quaternion{T}}) where {T} = T
real(q::Quaternion) = q.s
imag(q::Quaternion) = [q.v1, q.v2, q.v3]

(/)(q::Quaternion, x::Real) = Quaternion(q.s / x, q.v1 / x, q.v2 / x, q.v3 / x)

conj(q::Quaternion) = Quaternion(q.s, -q.v1, -q.v2, -q.v3, q.norm)
abs(q::Quaternion) = sqrt(q.s * q.s + q.v1 * q.v1 + q.v2 * q.v2 + q.v3 * q.v3)
float(q::Quaternion{T}) where T = convert(Quaternion{float(T)}, q)
abs_imag(q::Quaternion) = sqrt(q.v1 * q.v1 + q.v2 * q.v2 + q.v3 * q.v3)
abs2(q::Quaternion) = q.s * q.s + q.v1 * q.v1 + q.v2 * q.v2 + q.v3 * q.v3
inv(q::Quaternion) = q.norm ? conj(q) : conj(q) / abs2(q)

isreal(q::Quaternion) = iszero(q.v1) & iszero(q.v2) & iszero(q.v3)
isfinite(q::Quaternion) = q.norm | (isfinite(q.s) & isfinite(q.v1) & isfinite(q.v2) & isfinite(q.v3))
iszero(q::Quaternion) = ~q.norm & iszero(real(q)) & iszero(q.v1) & iszero(q.v2) & iszero(q.v3)
isnan(q::Quaternion) = isnan(real(q)) | isnan(q.v1) | isnan(q.v2) | isnan(q.v3)
isinf(q::Quaternion) = ~q.norm & (isinf(q.s) | isinf(q.v1) | isinf(q.v2) | isinf(q.v3))

function normalize(q::Quaternion)
    if (q.norm)
        return q
    end
    q = q / abs(q)
    Quaternion(q.s, q.v1, q.v2, q.v3, true)
end

function normalizea(q::Quaternion)
    if (q.norm)
        return (q, one(q.s))
    end
    a = abs(q)
    q = q / a
    (Quaternion(q.s, q.v1, q.v2, q.v3, true), a)
end

function normalizeq(q::Quaternion)
    a = abs(q)
    if a > 0
        q = q / a
        Quaternion(q.s, q.v1, q.v2, q.v3, true)
    else
        Quaternion(0.0, 1.0, 0.0, 0.0, true)
    end
end

(-)(q::Quaternion) = Quaternion(-q.s, -q.v1, -q.v2, -q.v3, q.norm)

(+)(q::Quaternion, w::Quaternion) =
    Quaternion(q.s + w.s, q.v1 + w.v1, q.v2 + w.v2, q.v3 + w.v3)

(-)(q::Quaternion, w::Quaternion) =
    Quaternion(q.s - w.s, q.v1 - w.v1, q.v2 - w.v2, q.v3 - w.v3)

(*)(q::Quaternion, w::Quaternion) = Quaternion(q.s * w.s - q.v1 * w.v1 - q.v2 * w.v2 - q.v3 * w.v3,
                                               q.s * w.v1 + q.v1 * w.s + q.v2 * w.v3 - q.v3 * w.v2,
                                               q.s * w.v2 - q.v1 * w.v3 + q.v2 * w.s + q.v3 * w.v1,
                                               q.s * w.v3 + q.v1 * w.v2 - q.v2 * w.v1 + q.v3 * w.s,
                                               q.norm && w.norm)
(/)(q::Quaternion, w::Quaternion) = q * inv(w)

(==)(q::Quaternion, w::Quaternion) = (q.s == w.s) & (q.v1 == w.v1) & (q.v2 == w.v2) & (q.v3 == w.v3) # ignore .norm field

angleaxis(q::Quaternion) = angle(q), axis(q)

angle(q::Quaternion) = 2 * atan(√(q.v1^2 + q.v2^2 + q.v3^2), q.s)

function axis(q::Quaternion)
    q = normalize(q)
    s = sin(angle(q) / 2)
    abs(s) > 0 ?
        [q.v1, q.v2, q.v3] / s :
        [1.0, 0.0, 0.0]
end

argq(q::Quaternion) = normalizeq(Quaternion(0, q.v1, q.v2, q.v3))

"""
    extend_analytic(f, q::Quaternion)

Evaluate the extension of the complex analytic function `f` to the quaternions at `q`.

Given ``q = s + a u``, where ``s`` is the real part, ``u`` is a pure unit quaternion,
and ``a \\ge 0`` is the magnitude of the imaginary part of ``q``,

```math
f(q) = \\Re(f(z)) + \\Im(f(z)) u,
```
is the extension of `f` to the quaternions, where ``z = a + s i`` is a complex analog to
``q``.

See Theorem 5 of [^Sudbery1970] for details.

[^Sudbery1970]
    Sudbery (1979). Quaternionic analysis. Mathematical Proceedings of the Cambridge 
    Philosophical Society,85, pp 199­225
    doi:[10.1017/S030500410005563](https://doi.org/10.1017/S0305004100055638)
"""
function extend_analytic(f, q::Quaternion)
    a = abs_imag(q)
    s = q.s
    z = complex(s, a)
    w = f(z)
    wr, wi = reim(w)
    scale = wi / a
    norm = _isexpfun(f) && iszero(s)
    if a > 0
        return Quaternion(wr, scale * q.v1, scale * q.v2, scale * q.v3, norm)
    else
        # q == real(q), so f(real(q)) may be real or complex, i.e. wi may be nonzero.
        # we choose to embed complex numbers in the quaternions by identifying the first
        # imaginary quaternion basis with the complex imaginary basis.
        return Quaternion(wr, oftype(scale, wi), zero(scale), zero(scale), norm)
    end
end

_isexpfun(::Union{typeof(exp),typeof(exp2),typeof(exp10)}) = true
_isexpfun(::Any) = false

"""
    cis(q::Quaternion)

Return ``\\exp(u * q)``, where ``u`` is the normalized imaginary part of `q`.

Let ``q = s + a u``, where ``s`` is the real part, ``u`` is a pure unit quaternion,
and ``a \\ge 0`` is the magnitude of the imaginary part of ``q``.

!!! Note
    This is the extension of `cis(z)` for complex `z` to the quaternions and is not
    equivalent to `exp(im * q)`. As a result, `cis(Quaternion(z)) ≠ cis(z)` when
    `imag(z) < 0`.
"""
cis(q::Quaternion)

if VERSION ≥ v"1.6"
    """
        cispi(q::Quaternion)

    Compute `cis(π * q)` more accurately.

    !!! Note
        This is not equivalent to `exp(π*im*q)`. See [cis(::Quaternion)](@ref) for details.
    """
    Base.cispi(q::Quaternion) = extend_analytic(cispi, q)
end

for f in (
    :sqrt, :exp, :exp2, :exp10, :expm1, :log2, :log10, :log1p, :cis,
    :sin, :cos, :tan, :asin, :acos, :atan, :sinh, :cosh, :tanh, :asinh, :acosh, :atanh,
    :csc, :sec, :cot, :acsc, :asec, :acot, :csch, :sech, :coth, :acsch, :asech, :acoth,
    :sinpi, :cospi,
)
    @eval Base.$f(q::Quaternion) = extend_analytic($f, q)
end

for f in (@static(VERSION ≥ v"1.6" ? (:sincos, :sincospi) : (:sincos,)))
    @eval begin
        function Base.$f(q::Quaternion)
            a = abs_imag(q)
            z = complex(q.s, a)
            s, c = $f(z)
            sr, si = reim(s)
            cr, ci = reim(c)
            sscale = si / a
            cscale = ci / a
            if a > 0
                return (
                    Quaternion(sr, sscale * q.v1, sscale * q.v2, sscale * q.v3),
                    Quaternion(cr, cscale * q.v1, cscale * q.v2, cscale * q.v3),
                )
            else
                return (
                    Quaternion(sr, oftype(sscale, si), zero(sscale), zero(sscale)),
                    Quaternion(cr, oftype(cscale, ci), zero(cscale), zero(cscale)),
                )
            end
        end
    end
end

function log(q::Quaternion)
    a = abs(q)
    M = abs_imag(q)
    theta = atan(M, q.s)
    scale = theta / ifelse(iszero(M), oneunit(M), M)
    return Quaternion(log(a), q.v1 * scale, q.v2 * scale, q.v3 * scale)
end

(^)(q::Quaternion, w::Quaternion) = exp(w * log(q))

function linpol(p::Quaternion, q::Quaternion, t::Real)
    p = normalize(p)
    q = normalize(q)
    qm = -q
    if abs(p - q) > abs(p - qm)
        q = qm
    end
    c = p.s * q.s + p.v1 * q.v1 + p.v2 * q.v2 + p.v3 * q.v3
    if c < 1.0
        o = acos(c)
        s = sin(o)
        sp = sin((1 - t) * o) / s
        sq = sin(t * o) / s
    else
        sp = 1 - t
        sq = t
    end
    Quaternion(sp * p.s  + sq * q.s,
                sp * p.v1 + sq * q.v1,
                sp * p.v2 + sq * q.v2,
                sp * p.v3 + sq * q.v3, true)
end

quatrand(rng = Random.GLOBAL_RNG)  = quat(randn(rng), randn(rng), randn(rng), randn(rng))
nquatrand(rng = Random.GLOBAL_RNG) = normalize(quatrand(rng))

function rand(rng::AbstractRNG, ::Random.SamplerType{Quaternion{T}}) where {T<:Real}
    Quaternion{T}(rand(rng, T), rand(rng, T), rand(rng, T), rand(rng, T), false)
end

function randn(rng::AbstractRNG, ::Type{Quaternion{T}}) where {T<:AbstractFloat}
    Quaternion{T}(
        randn(rng, T) * 1//2,
        randn(rng, T) * 1//2,
        randn(rng, T) * 1//2,
        randn(rng, T) * 1//2,
        false,
    )
end

## Rotations

function qrotation(axis::Vector{T}, theta) where {T <: Real}
    if length(axis) != 3
        error("Must be a 3-vector")
    end
    normaxis = norm(axis)
    if iszero(normaxis)
        normaxis = oneunit(normaxis)
        theta = zero(theta)
    end
    s,c = sincos(theta / 2)
    scaleby = s / normaxis
    Quaternion(c, scaleby * axis[1], scaleby * axis[2], scaleby * axis[3], true)
end

# Variant of the above where norm(rotvec) encodes theta
function qrotation(rotvec::Vector{T}) where {T <: Real}
    if length(rotvec) != 3
        error("Must be a 3-vector")
    end
    theta = norm(rotvec)
    s,c = sincos(theta / 2)
    scaleby = s / (iszero(theta) ? one(theta) : theta)
    Quaternion(c, scaleby * rotvec[1], scaleby * rotvec[2], scaleby * rotvec[3], true)
end

function qrotation(dcm::Matrix{T}) where {T<:Real}
    # See https://arxiv.org/pdf/math/0701759.pdf
    a2 = 1 + dcm[1,1] + dcm[2,2] + dcm[3,3]
    b2 = 1 + dcm[1,1] - dcm[2,2] - dcm[3,3]
    c2 = 1 - dcm[1,1] + dcm[2,2] - dcm[3,3]
    d2 = 1 - dcm[1,1] - dcm[2,2] + dcm[3,3]

    if a2 >= max(b2, c2, d2)
        a = sqrt(a2)/2
        return Quaternion(a, (dcm[3,2]-dcm[2,3])/4a, (dcm[1,3]-dcm[3,1])/4a, (dcm[2,1]-dcm[1,2])/4a)
    elseif b2 >= max(c2, d2)
        b = sqrt(b2)/2
        return Quaternion((dcm[3,2]-dcm[2,3])/4b, b, (dcm[2,1]+dcm[1,2])/4b, (dcm[1,3]+dcm[3,1])/4b)
    elseif c2 >= d2
        c = sqrt(c2)/2
        return Quaternion((dcm[1,3]-dcm[3,1])/4c, (dcm[2,1]+dcm[1,2])/4c, c, (dcm[3,2]+dcm[2,3])/4c)
    else
        d = sqrt(d2)/2
        return Quaternion((dcm[2,1]-dcm[1,2])/4d, (dcm[1,3]+dcm[3,1])/4d, (dcm[3,2]+dcm[2,3])/4d, d)
    end
end

function qrotation(dcm::Matrix{T}, qa::Quaternion) where {T<:Real}
    q = qrotation(dcm)
    abs(q-qa) < abs(q+qa) ? q : -q
end

rotationmatrix(q::Quaternion) = rotationmatrix_normalized(normalize(q))

function rotationmatrix_normalized(q::Quaternion)
    sx, sy, sz = 2q.s * q.v1, 2q.s * q.v2, 2q.s * q.v3
    xx, xy, xz = 2q.v1^2, 2q.v1 * q.v2, 2q.v1 * q.v3
    yy, yz, zz = 2q.v2^2, 2q.v2 * q.v3, 2q.v3^2
    [1 - (yy + zz)     xy - sz     xz + sy;
        xy + sz   1 - (xx + zz)    yz - sx;
        xz - sy      yz + sx  1 - (xx + yy)]
end


function slerp(qa::Quaternion{T}, qb::Quaternion{T}, t::T) where {T}
    # http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/
    coshalftheta = qa.s * qb.s + qa.v1 * qb.v1 + qa.v2 * qb.v2 + qa.v3 * qb.v3;

    if coshalftheta < 0
        qm = -qb
        coshalftheta = -coshalftheta
    else
        qm = qb
    end
    abs(coshalftheta) >= 1.0 && return qa

    halftheta    = acos(coshalftheta)
    sinhalftheta = sqrt(one(T) - coshalftheta * coshalftheta)

    if abs(sinhalftheta) < 0.001
        return Quaternion(
            T(0.5) * (qa.s  + qb.s),
            T(0.5) * (qa.v1 + qb.v1),
            T(0.5) * (qa.v2 + qb.v2),
            T(0.5) * (qa.v3 + qb.v3),
        )
    end

    ratio_a = sin((one(T) - t) * halftheta) / sinhalftheta
    ratio_b = sin(t * halftheta) / sinhalftheta

    Quaternion(
        qa.s  * ratio_a + qm.s  * ratio_b,
        qa.v1 * ratio_a + qm.v1 * ratio_b,
        qa.v2 * ratio_a + qm.v2 * ratio_b,
        qa.v3 * ratio_a + qm.v3 * ratio_b,
    )
end

function sylvester(a::Quaternion{T}, b::Quaternion{T}, c::Quaternion{T}) where {T<:Real}
    isreal(a) && return sylvester(real(a), b, c)
    isreal(b) && return sylvester(a, real(b), c)
    abs2a = abs2(a)
    abs2b = abs2(b)
    if abs2a > abs2b
        inva = conj(a) / abs2a
        d1 = -2real(b) - a - inva * abs2b
        x = d1 \ (c + inva * c * conj(b))
    else
        invb = conj(b) / abs2b
        d2 = -2real(a) - b - invb * abs2a
        x = (c + conj(a) * c * invb) / d2
    end
    return x
end
sylvester(a::Quaternion, b::Quaternion, c::Quaternion) = sylvester(promote(a, b, c)...)
sylvester(a::Quaternion, b::Quaternion, c::Real) = sylvester(promote(a, b, c)...)
# if either a or b commute with x, use a simpler expression
sylvester(a::Real, b::Real, c::Quaternion) = c / -(a + b)
sylvester(a::Real, b::Quaternion, c::Quaternion) = c / -(a + b)
sylvester(a::Quaternion, b::Real, c::Quaternion) = -(a + b) \ c
sylvester(a::Real, b::Quaternion, c::Real) = -c / (a + b)
sylvester(a::Quaternion, b::Real, c::Real) = (a + b) \ -c

function lyap(a::Quaternion{T}, c::Quaternion{T}) where {T<:Real}
    # if a commutes with c, use a simpler expression
    (isreal(a) || isreal(c)) && return c / -2real(a)
    return (c + a \ c * a) / -4real(a)
end
lyap(a::Quaternion, c::Quaternion) = lyap(promote(a, c)...)
# if a commutes with c, use a simpler expression
lyap(a::Real, c::Quaternion) = c / -2a
lyap(a::Quaternion, c::Real) = c / -2real(a)
