"""
    Quaternion{T<:Real} <: Number

Quaternion number type with real and imaginary parts of type `T`.

`QuaternionF16`, `QuaternionF32`, and `QuaternionF64` are aliases for
`Quaternion{Float16}`, `Quaternion{Float32}`, and `Quaternion{Float64}`, respectively.

See also: [`quat`](@ref), [`real`](@ref), [`imag_part`](@ref).
"""
struct Quaternion{T<:Real} <: Number
    s::T
    v1::T
    v2::T
    v3::T
end

const QuaternionF16 = Quaternion{Float16}
const QuaternionF32 = Quaternion{Float32}
const QuaternionF64 = Quaternion{Float64}

Quaternion{T}(x::Real) where {T<:Real} = Quaternion(convert(T, x))
Quaternion{T}(q::Quaternion) where {T<:Real} = Quaternion{T}(q.s, q.v1, q.v2, q.v3)
Quaternion(s::Real, v1::Real, v2::Real, v3::Real) = Quaternion(promote(s, v1, v2, v3)...)
Quaternion(x::Real) = Quaternion(x, zero(x), zero(x), zero(x))

Base.promote_rule(::Type{Quaternion{T}}, ::Type{S}) where {T <: Real, S <: Real} = Quaternion{promote_type(T, S)}
Base.promote_rule(::Type{Quaternion{T}}, ::Type{Quaternion{S}}) where {T <: Real, S <: Real} = Quaternion{promote_type(T, S)}

"""
    quat(w, [x, y, z])

Convert real numbers or arrays to quaternion. `x, y, z` defaults to zero.

# Examples
```jldoctest
julia> quat(7)
Quaternion{Int64}(7, 0, 0, 0)

julia> quat(1.0, 2, 3, 4)
QuaternionF64(1.0, 2.0, 3.0, 4.0)

julia> quat([1, 2, 3])
3-element Vector{Quaternion{Int64}}:
 Quaternion{Int64}(1, 0, 0, 0)
 Quaternion{Int64}(2, 0, 0, 0)
 Quaternion{Int64}(3, 0, 0, 0)
```
"""
quat

quat(q::Quaternion) = q
quat(s::Real) = Quaternion(s)
quat(s::Real, v1::Real, v2::Real, v3::Real) = Quaternion(s, v1, v2, v3)

## Array operations on quaternions ##
quat(A::AbstractArray{<:Quaternion}) = A
function quat(A::AbstractArray{T}) where T
    if !isconcretetype(T)
        error("`quat` not defined on abstractly-typed arrays; please convert to a more specific type")
    end
    convert(AbstractArray{typeof(quat(zero(T)))}, A)
end

"""
    quat(T::Type)

Return an appropriate type that can represent a value of type `T` as a quaternion.
Equivalent to `typeof(quat(zero(T)))`.

# Examples
```jldoctest
julia> quat(Quaternion{Int})
Quaternion{Int64}

julia> quat(Int)
Quaternion{Int64}
```
"""
quat(::Type{T}) where {T<:Real} = Quaternion{T}
quat(::Type{Quaternion{T}}) where {T<:Real} = Quaternion{T}

quat(::Missing) = missing
# Same definitioin as in Base: https://github.com/JuliaLang/julia/blob/v1.9.3/base/missing.jl#L111-L115
quat(::Type{Missing}) = Missing
function quat(::Type{Union{T, Missing}}) where T
    T === Any && throw(MethodError(quat, (Any,)))  # To prevent StackOverflowError
    Union{quat(T), Missing}
end

"""
    real(T::Type{<:Quaternion})

Return the type that represents the real part of a value of type `T`.
e.g: for `T == Quaternion{R}`, returns `R`.
Equivalent to `typeof(real(zero(T)))`.

# Examples
```jldoctest
julia> real(Quaternion{Int})
Int64
```
"""
Base.real(::Type{Quaternion{T}}) where {T} = T

"""
    real(q::Quaternion)

Return the real part of the quaternion `q`.

See also: [`imag_part`](@ref), [`quat`](@ref)

# Examples
```jldoctest
julia> real(quat(1,2,3,4))
1
```
"""
Base.real(q::Quaternion) = q.s

"""
    real(A::AbstractArray{<:Quaternion})

Return an array containing the real part of each quaternion in `A`.

# Examples
```jldoctest
julia> real([quat(5,6,7,8), 9])
2-element Vector{Int64}:
 5
 9
```
"""
Base.real(::AbstractArray{<:Quaternion})

"""
    imag_part(q::Quaternion{T}) -> NTuple{3, T}

Return the imaginary part of the quaternion `q`.

Note that this function is different from `Base.imag`, which returns `Real` for complex numbers.

See also: [`real`](@ref), [`conj`](@ref).

# Examples
```jldoctest
julia> imag_part(Quaternion(1,2,3,4))
(2, 3, 4)
```
"""
imag_part(q::Quaternion) = (q.v1, q.v2, q.v3)

Base.:/(q::Quaternion, x::Real) = Quaternion(q.s / x, q.v1 / x, q.v2 / x, q.v3 / x)
Base.:*(q::Quaternion, x::Real) = Quaternion(q.s * x, q.v1 * x, q.v2 * x, q.v3 * x)
Base.:*(x::Real, q::Quaternion) = q * x

"""
    conj(q::Quaternion)

Compute the quaternion conjugate of a quaternion `q`.

# Examples
```jldoctest
julia> conj(Quaternion(1,2,3,4))
Quaternion{Int64}(1, -2, -3, -4)
```
"""
Base.conj(q::Quaternion) = Quaternion(q.s, -q.v1, -q.v2, -q.v3)
function Base.abs(q::Quaternion)
    a = max(abs(q.s), abs(q.v1), abs(q.v2), abs(q.v3))
    if isnan(a) && isinf(q)
        return typeof(a)(Inf)
    elseif iszero(a) || isinf(a)
        return a
    else
        return sqrt(abs2(q / a)) * a
    end
end
Base.float(q::Quaternion{T}) where T = convert(Quaternion{float(T)}, q)
function abs_imag(q::Quaternion)
    a = max(abs(q.v1), abs(q.v2), abs(q.v3))
    if isnan(a) && (isinf(q.v1) | isinf(q.v2) | isinf(q.v3))
        return oftype(a, Inf)
    elseif iszero(a) || isinf(a)
        return a
    else
        return sqrt((q.v1 / a)^2 + (q.v2 / a)^2 + (q.v3 / a)^2) * a
    end
end
Base.abs2(q::Quaternion) = RealDot.realdot(q,q)

"""
    inv(q::Quaternion)

Return the multiplicative inverse of `q::Quaternion`, such that `q*inv(q)` or `inv(q)*q`
yields `one(q)` (the multiplicative identity) up to roundoff errors.

# Examples
```jldoctest
julia> inv(quat(1))
QuaternionF64(1.0, -0.0, -0.0, -0.0)

julia> inv(quat(1, 2, 0, 0))
QuaternionF64(0.2, -0.4, -0.0, -0.0)

julia> inv(quat(2//3))
Quaternion{Rational{Int64}}(3//2, 0//1, 0//1, 0//1)
```
"""
function Base.inv(q::Quaternion)
    if isinf(q)
        return quat(
            copysign(zero(q.s), q.s),
            flipsign(-zero(q.v1), q.v1),
            flipsign(-zero(q.v2), q.v2),
            flipsign(-zero(q.v3), q.v3),
        )
    end
    a = max(abs(q.s), abs(q.v1), abs(q.v2), abs(q.v3))
    p = q / a
    iq = conj(p) / (a * abs2(p))
    return iq
end

Base.isreal(q::Quaternion) = iszero(q.v1) & iszero(q.v2) & iszero(q.v3)
Base.isfinite(q::Quaternion) = isfinite(q.s) & isfinite(q.v1) & isfinite(q.v2) & isfinite(q.v3)
Base.iszero(q::Quaternion) = iszero(real(q)) & iszero(q.v1) & iszero(q.v2) & iszero(q.v3)
Base.isnan(q::Quaternion) = isnan(real(q)) | isnan(q.v1) | isnan(q.v2) | isnan(q.v3)
Base.isinf(q::Quaternion) = isinf(q.s) | isinf(q.v1) | isinf(q.v2) | isinf(q.v3)
Base.isinteger(q::Quaternion) = isinteger(real(q)) & isreal(q)

# included strictly for documentation; the base implementation is sufficient
"""
    sign(q::Quaternion) -> Quaternion

Return zero if `q==0` and ``q/|q|`` otherwise.

# Examples
```jldoctest
julia> sign(Quaternion(4, 0, 0, 0))
QuaternionF64(1.0, 0.0, 0.0, 0.0)

julia> sign(Quaternion(1, 0, 1, 0))
QuaternionF64(0.7071067811865475, 0.0, 0.7071067811865475, 0.0)
```
"""
sign(::Quaternion)

Base.:-(q::Quaternion) = Quaternion(-q.s, -q.v1, -q.v2, -q.v3)

Base.:+(q::Quaternion, w::Quaternion) =
    Quaternion(q.s + w.s, q.v1 + w.v1, q.v2 + w.v2, q.v3 + w.v3)

Base.:-(q::Quaternion, w::Quaternion) =
    Quaternion(q.s - w.s, q.v1 - w.v1, q.v2 - w.v2, q.v3 - w.v3)

function Base.:*(q::Quaternion, w::Quaternion)
    s  = (q.s * w.s - q.v2 * w.v2) - (q.v1 * w.v1 + q.v3 * w.v3)
    v1 = (q.s * w.v1 + q.v1 * w.s) + (q.v2 * w.v3 - q.v3 * w.v2)
    v2 = (q.s * w.v2 + q.v2 * w.s) + (q.v3 * w.v1 - q.v1 * w.v3)
    v3 = (q.s * w.v3 + q.v3 * w.s) + (q.v1 * w.v2 - q.v2 * w.v1)
    return Quaternion(s, v1, v2, v3)
end

function Base.:/(q::Quaternion{T}, w::Quaternion{T}) where T
    # handle over/underflow while matching the behavior of /(a::Complex, b::Complex)
    a = max(abs(w.s), abs(w.v1), abs(w.v2), abs(w.v3))
    if isinf(w)
        if isfinite(q)
            return quat(
                zero(T)*sign(q.s)*sign(w.s),
                -zero(T)*sign(q.v1)*sign(w.v1),
                -zero(T)*sign(q.v2)*sign(w.v2),
                -zero(T)*sign(q.v3)*sign(w.v3),
            )
        end
        return quat(T(NaN), T(NaN), T(NaN), T(NaN))
    end
    p = w / a
    return (q * conj(p)) / RealDot.realdot(w, p)
end

Base.://(x::Quaternion, y::Real) = quat(real(x)//y, imag_part(x).//y...)
Base.://(x::Number, y::Quaternion) = x*conj(y)//abs2(y)

Base.:(==)(q::Quaternion, w::Quaternion) = (q.s == w.s) & (q.v1 == w.v1) & (q.v2 == w.v2) & (q.v3 == w.v3)
function Base.isequal(q::Quaternion, w::Quaternion)
    isequal(q.s, w.s) & isequal(q.v1, w.v1) & isequal(q.v2, w.v2) & isequal(q.v3, w.v3)
end

if UInt === UInt64
    const h_imags = 0xfa29df508725c1bf, 0x5e6c4b26a400626d, 0x8dfd375bac9f9403
else
    const h_imags = 0x62802719, 0x5a9b072e, 0x209019b1
end
const hash_0_imags = hash.(0, h_imags)

function Base.hash(z::Quaternion, h::UInt)
    hash(real(z), h ⊻ xor(xor.(hash.(imag_part(z), h_imags), hash_0_imags)...))
end

"""
    extend_analytic(f, q::Quaternion)

Evaluate the extension of the complex analytic function `f` to the quaternions at `q`.

Given ``q = s + a u``, where ``s`` is the real part, ``u`` is a pure unit quaternion,
and ``a \\ge 0`` is the magnitude of the imaginary part of ``q``,

```math
f(q) = \\Re(f(z)) + \\Im(f(z)) u,
```
is the extension of `f` to the quaternions, where ``z = s + a i`` is a complex analog to
``q``.

See Theorem 5 of [^Sudbery1970] for details.

[^Sudbery1970]:
    Sudbery (1979). Quaternionic analysis. Mathematical Proceedings of the Cambridge
    Philosophical Society,85, pp 199225
    doi:[10.1017/S030500410005563](https://doi.org/10.1017/S0305004100055638)
"""
function extend_analytic(f, q::Quaternion)
    a = abs_imag(q)
    s = q.s
    z = complex(s, a)
    w = f(z)
    wr, wi = reim(w)
    scale = wi / a
    if a > 0
        return Quaternion(wr, scale * q.v1, scale * q.v2, scale * q.v3)
    else
        # q == real(q), so f(real(q)) may be real or complex, i.e. wi may be nonzero.
        # we choose to embed complex numbers in the quaternions by identifying the first
        # imaginary quaternion basis with the complex imaginary basis.
        return Quaternion(wr, oftype(scale, wi), zero(scale), zero(scale))
    end
end

for f in (
    :sqrt, :exp, :exp2, :exp10, :expm1, :log2, :log10, :log1p,
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

# this implementation is roughly 2x as fast as extend_analytic(log, q)
function Base.log(q::Quaternion)
    a = abs_imag(q)
    theta = atan(a, q.s)
    scale = theta / a
    if a > 0
        return Quaternion(log(abs(q)), scale * q.v1, scale * q.v2, scale * q.v3)
    else
        # q == real(q), so f(real(q)) may be real or complex.
        # we choose to embed complex numbers in the quaternions by identifying the first
        # imaginary quaternion basis with the complex imaginary basis.
        return Quaternion(log(abs(q.s)), oftype(scale, theta), zero(scale), zero(scale))
    end
end

Base.:^(q::Quaternion, w::Quaternion) = exp(w * log(q))

Base.@deprecate quatrand(rng::AbstractRNG=Random.GLOBAL_RNG) randn(rng, QuaternionF64) * 2
Base.@deprecate nquatrand(rng = Random.GLOBAL_RNG) sign(randn(rng, QuaternionF64))

function Base.rand(rng::AbstractRNG, ::Random.SamplerType{Quaternion{T}}) where {T<:Real}
    Quaternion{T}(rand(rng, T), rand(rng, T), rand(rng, T), rand(rng, T))
end

function Base.randn(rng::AbstractRNG, ::Type{Quaternion{T}}) where {T<:AbstractFloat}
    Quaternion{T}(
        randn(rng, T) / 2,
        randn(rng, T) / 2,
        randn(rng, T) / 2,
        randn(rng, T) / 2,
    )
end

"""
    slerp(qa::Quaternion, qb::Quaternion, t::Real)

Spherical linear interpolation (Slerp) between the inputs `qa` and `qb`.
Since the input is normalized inside the function, the absolute value of the return value will be 1.

# Examples
```jldoctest
julia> using Quaternions

julia> qa = Quaternion(1,0,0,0)
Quaternion{Int64}(1, 0, 0, 0)

julia> qb = Quaternion(0,1,0,0)
Quaternion{Int64}(0, 1, 0, 0)

julia> slerp(qa, qb, 0.6)
QuaternionF64(0.5877852522924731, 0.8090169943749475, 0.0, 0.0)

julia> ans ≈ Quaternion(cospi(0.3), sinpi(0.3), 0, 0)
true
```
"""
@inline function slerp(qa0::Quaternion{T}, qb0::Quaternion{T}, t::T) where T<:Real
    # http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/
    iszero(qa0) && throw(DomainError(qa0, "The input quaternion must be non-zero."))
    iszero(qb0) && throw(DomainError(qb0, "The input quaternion must be non-zero."))
    qa = qa0 / abs(qa0)
    qb = qb0 / abs(qb0)
    coshalftheta = RealDot.realdot(qa, qb)

    if coshalftheta < 0
        qb = -qb
        coshalftheta = -coshalftheta
    end

    if coshalftheta < 1
        halftheta    = acos(coshalftheta)
        sinhalftheta = sqrt(1 - coshalftheta^2)

        ratio_a = sin((1 - t) * halftheta) / sinhalftheta
        ratio_b = sin(t * halftheta) / sinhalftheta
    else
        ratio_a = float(1 - t)
        ratio_b = float(t)
    end

    return Quaternion(
        qa.s  * ratio_a + qb.s  * ratio_b,
        qa.v1 * ratio_a + qb.v1 * ratio_b,
        qa.v2 * ratio_a + qb.v2 * ratio_b,
        qa.v3 * ratio_a + qb.v3 * ratio_b,
    )
end

function slerp(qa::Quaternion{Ta}, qb::Quaternion{Tb}, t::T) where {Ta, Tb, T}
    S = promote_type(Ta,Tb,T)
    return slerp(Quaternion{S}(qa),Quaternion{S}(qb),S(t))
end

function LinearAlgebra.sylvester(a::Quaternion{T}, b::Quaternion{T}, c::Quaternion{T}) where {T<:Real}
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
LinearAlgebra.sylvester(a::Quaternion, b::Quaternion, c::Quaternion) = sylvester(promote(a, b, c)...)
LinearAlgebra.sylvester(a::Quaternion, b::Quaternion, c::Real) = sylvester(promote(a, b, c)...)
# if either a or b commute with x, use a simpler expression
LinearAlgebra.sylvester(a::Real, b::Real, c::Quaternion) = c / -(a + b)
LinearAlgebra.sylvester(a::Real, b::Quaternion, c::Quaternion) = c / -(a + b)
LinearAlgebra.sylvester(a::Quaternion, b::Real, c::Quaternion) = -(a + b) \ c
LinearAlgebra.sylvester(a::Real, b::Quaternion, c::Real) = -c / (a + b)
LinearAlgebra.sylvester(a::Quaternion, b::Real, c::Real) = (a + b) \ -c

function LinearAlgebra.lyap(a::Quaternion{T}, c::Quaternion{T}) where {T<:Real}
    # if a commutes with c, use a simpler expression
    (isreal(a) || isreal(c)) && return c / -2real(a)
    return (c + a \ c * a) / -4real(a)
end
LinearAlgebra.lyap(a::Quaternion, c::Quaternion) = lyap(promote(a, c)...)
# if a commutes with c, use a simpler expression
LinearAlgebra.lyap(a::Real, c::Quaternion) = c / -2a
LinearAlgebra.lyap(a::Quaternion, c::Real) = c / -2real(a)

## RealDot
# ordering chosen so that real(p'q) == real(q'p) == realdot(p, q) == realdot(q, p), i.e. exact equality
@inline RealDot.realdot(p::Quaternion, q::Quaternion) = (p.s * q.s + p.v2 * q.v2) + (p.v1 * q.v1 + p.v3 * q.v3)
Base.widen(::Type{Quaternion{T}}) where {T} = Quaternion{widen(T)}

Base.flipsign(x::Quaternion, y::Real) = ifelse(signbit(y), -x, x)

function Base.read(io::IO, ::Type{Quaternion{T}}) where T<:Real
    return Quaternion{T}(ntuple(_ -> read(io, T), Val(4))...)
end
Base.write(io::IO, q::Quaternion) = write(io, real(q), imag_part(q)...)

Base.big(::Type{Quaternion{T}}) where {T<:Real} = Quaternion{big(T)}
Base.big(z::Quaternion{T}) where {T<:Real} = Quaternion{big(T)}(z)

"""
    round(q::Quaternion[, RoundingModeReal, [RoundingModeImaginary]]; kwargs...)
    round(q::Quaternion, RoundingModeReal,
          RoundingModeImaginary1, RoundingModeImaginary2, RoundingModeImaginary3; kwargs...)

Return the nearest integral value of the same type as the quaternion-valued `q` to `q`,
breaking ties using the specified `RoundingMode`s.

The first `RoundingMode` is used for rounding the real part while the second is used
for rounding the imaginary parts. Alternatively, a `RoundingMode` may be provided for each
part.

The `kwargs` are the same as those for `round(::Real[, RoundingMode]; kwargs...)`.

# Example
```jldoctest
julia> round(quat(3.14, 4.5, 8.3, -2.8))
QuaternionF64(3.0, 4.0, 8.0, -3.0)
```
"""
function Base.round(
    q::Quaternion,
    rs::RoundingMode=RoundNearest,
    rv::RoundingMode=rs;
    kwargs...,
)
    return round(q, rs, rv, rv, rv; kwargs...)
end
function Base.round(
    q::Quaternion,
    rs::RoundingMode,
    rv1::RoundingMode,
    rv2::RoundingMode,
    rv3::RoundingMode;
    kwargs...,
)
    return Quaternion(
        round(real(q), rs; kwargs...),
        round(q.v1, rv1; kwargs...),
        round(q.v2, rv2; kwargs...),
        round(q.v3, rv3; kwargs...),
    )
end
