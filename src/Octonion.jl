struct Octonion{T<:Real} <: Number
  s::T
  v1::T
  v2::T
  v3::T
  v4::T
  v5::T
  v6::T
  v7::T
  norm::Bool
end

Octonion{T}(x::Real) where {T<:Real} = Octonion(convert(T, x))
Octonion{T}(x::Complex) where {T<:Real} = Octonion(convert(Complex{T}, x))
Octonion{T}(q::Quaternion) where {T<:Real} = Octonion(convert(Quaternion{T}, q))
Octonion{T}(o::Octonion) where {T<:Real} =
  Octonion{T}(o.s, o.v1, o.v2, o.v3, o.v4, o.v5, o.v6, o.v7, o.norm)

Octonion(s::Real, v1::Real, v2::Real, v3::Real, v4::Real, v5::Real, v6::Real, v7::Real, n::Bool = false) =
  Octonion(promote(s, v1, v2, v3, v4, v5, v6, v7)..., n)
Octonion(x::Real) = Octonion(x, zero(x), zero(x), zero(x), zero(x), zero(x), zero(x), zero(x), abs(x) == one(x))
Octonion(z::Complex) = Octonion(z.re, z.im, zero(z.re), zero(z.re), zero(z.re), zero(z.re), zero(z.re), zero(z.re), abs(z) == one(z.re))
Octonion(q::Quaternion) = Octonion(q.s, q.v1, q.v2, q.v3, zero(q.s), zero(q.s), zero(q.s), zero(q.s), q.norm)
Octonion(s::Real, a::Vector) = Octonion(s, a[1], a[2], a[3], a[4], a[5], a[6], a[7])
Octonion(a::Vector) = Octonion(0, a[1], a[2], a[3], a[4], a[5], a[6], a[7])

const OctonionF16 = Octonion{Float16}
const OctonionF32 = Octonion{Float32}
const OctonionF64 = Octonion{Float64}

promote_rule(::Type{Octonion{T}}, ::Type{S}) where {T <: Real, S <: Real} = Octonion{promote_type(T, S)}
promote_rule(::Type{Octonion{T}}, ::Type{Complex{S}}) where {T <: Real, S <: Real} = Octonion{promote_type(T, S)}
promote_rule(::Type{Octonion{T}}, ::Type{Quaternion{S}}) where {T <: Real, S <: Real} = Octonion{promote_type(T, S)}
promote_rule(::Type{Octonion{T}}, ::Type{Octonion{S}}) where {T <: Real, S <: Real} = Octonion{promote_type(T, S)}

octo(p, v1, v2, v3, v4, v5, v6, v7) = Octonion(p, v1, v2, v3, v4, v5, v6, v7)
octo(p, v1, v2, v3, v4, v5, v6, v7, n) = Octonion(p, v1, v2, v3, v4, v5, v6, v7, n)
octo(x) = Octonion(x)
octo(s, a) = Octonion(s, a)

function show(io::IO, o::Octonion)
  pm(x) = x < 0 ? " - $(-x)" : " + $x"
  print(io, o.s, pm(o.v1), "im", pm(o.v2), "jm", pm(o.v3), "km", pm(o.v4), "ilm", pm(o.v5), "jlm", pm(o.v6), "klm", pm(o.v7), "lm")
end

real(o::Octonion) = o.s
imag_part(o::Octonion) = (o.v1, o.v2, o.v3, o.v4, o.v5, o.v6, o.v7)
@deprecate imag(o::Octonion) collect(imag_part(o)) false

(/)(o::Octonion, x::Real) = Octonion(o.s / x, o.v1 / x, o.v2 / x, o.v3 / x, o.v4 / x, o.v5 / x, o.v6 / x, o.v7 / x)
(*)(o::Octonion, x::Real) = Octonion(o.s * x, o.v1 * x, o.v2 * x, o.v3 * x, o.v4 * x, o.v5 * x, o.v6 * x, o.v7 * x)
(*)(x::Real, o::Octonion) = o * x

conj(o::Octonion) = Octonion(o.s, -o.v1, -o.v2, -o.v3, -o.v4, -o.v5, -o.v6, -o.v7, o.norm)
abs(o::Octonion) = sqrt(abs2(o))
float(q::Octonion{T}) where T = convert(Octonion{float(T)}, q)
abs_imag(o::Octonion) = sqrt((o.v4 * o.v4 + (o.v2 * o.v2 + o.v6 * o.v6)) + ((o.v1 * o.v1 + o.v5 * o.v5) + (o.v3 * o.v3 + o.v7 * o.v7))) # ordered to match abs2
abs2(o::Octonion) = ((o.s * o.s + o.v4 * o.v4) + (o.v2 * o.v2 + o.v6 * o.v6)) + ((o.v1 * o.v1  + o.v5 * o.v5) + (o.v3 * o.v3 + o.v7 * o.v7))
inv(o::Octonion) = o.norm ? conj(o) : conj(o) / abs2(o)

isreal(o::Octonion) = iszero(o.v1) & iszero(o.v2) & iszero(o.v3) & iszero(o.v4) & iszero(o.v5) & iszero(o.v6) & iszero(o.v7)
isfinite(o::Octonion) = o.norm | (isfinite(real(o)) & isfinite(o.v1) & isfinite(o.v2) & isfinite(o.v3) & isfinite(o.v4) & isfinite(o.v5) & isfinite(o.v6) & isfinite(o.v7))
iszero(o::Octonion) = ~o.norm & iszero(real(o)) & iszero(o.v1) & iszero(o.v2) & iszero(o.v3) & iszero(o.v4) & iszero(o.v5) & iszero(o.v6) & iszero(o.v7)
isnan(o::Octonion) = isnan(real(o)) | isnan(o.v1) | isnan(o.v2) | isnan(o.v3) | isnan(o.v4) | isnan(o.v5) | isnan(o.v6) | isnan(o.v7)
isinf(o::Octonion) = ~o.norm & (isinf(real(o)) | isinf(o.v1) | isinf(o.v2) | isinf(o.v3) | isinf(o.v4) | isinf(o.v5) | isinf(o.v6) | isinf(o.v7))

function normalize(o::Octonion)
  if (o.norm)
    return o
  end
  o = o / abs(o)
  Octonion(o.s, o.v1, o.v2, o.v3, o.v4, o.v5, o.v6, o.v7, true)
end

function normalizea(o::Octonion)
  if (o.norm)
    return (o, one(o.s))
  end
  a = abs(o)
  o = o / a
  (Octonion(o.s, o.v1, o.v2, o.v3, o.v4, o.v5, o.v6, o.v7, true), a)
end

(-)(o::Octonion) = Octonion(-o.s, -o.v1, -o.v2, -o.v3, -o.v4, -o.v5, -o.v6, -o.v7, o.norm)

(+)(o::Octonion, w::Octonion) = Octonion(o.s + w.s,
                                            o.v1 + w.v1,
                                            o.v2 + w.v2,
                                            o.v3 + w.v3,
                                            o.v4 + w.v4,
                                            o.v5 + w.v5,
                                            o.v6 + w.v6,
                                            o.v7 + w.v7)

(-)(o::Octonion, w::Octonion) = Octonion(o.s - w.s,
                                            o.v1 - w.v1,
                                            o.v2 - w.v2,
                                            o.v3 - w.v3,
                                            o.v4 - w.v4,
                                            o.v5 - w.v5,
                                            o.v6 - w.v6,
                                            o.v7 - w.v7)

function (*)(o::Octonion, w::Octonion)
    s  = ((o.s * w.s - o.v4 * w.v4) - (o.v2 * w.v2 + o.v6 * w.v6)) - ((o.v1 * w.v1 + o.v5 * w.v5) + (o.v3 * w.v3 + o.v7 * w.v7))
    v1 = ((o.s * w.v1 + o.v1 * w.s) + (o.v6 * w.v5 - o.v5 * w.v6)) + ((o.v2 * w.v3 - o.v3 * w.v2) + (o.v7 * w.v4 - o.v4 * w.v7))
    v2 = ((o.s * w.v2 + o.v2 * w.s) + (o.v4 * w.v6 - o.v6 * w.v4)) + ((o.v3 * w.v1 - o.v1 * w.v3) + (o.v7 * w.v5 - o.v5 * w.v7))
    v3 = ((o.s * w.v3 + o.v3 * w.s) + (o.v5 * w.v4 - o.v4 * w.v5)) + ((o.v1 * w.v2 - o.v2 * w.v1) + (o.v7 * w.v6 - o.v6 * w.v7))
    v4 = ((o.s * w.v4 + o.v4 * w.s) + (o.v3 * w.v5 - o.v5 * w.v3)) + ((o.v1 * w.v7 - o.v7 * w.v1) + (o.v6 * w.v2 - o.v2 * w.v6))
    v5 = ((o.s * w.v5 + o.v5 * w.s) + (o.v2 * w.v7 - o.v7 * w.v2)) + ((o.v1 * w.v6 - o.v6 * w.v1) + (o.v4 * w.v3 - o.v3 * w.v4))
    v6 = ((o.s * w.v6 + o.v6 * w.s) + (o.v3 * w.v7 - o.v7 * w.v3)) + ((o.v2 * w.v4 - o.v4 * w.v2) + (o.v5 * w.v1 - o.v1 * w.v5))
    v7 = ((o.s * w.v7 + o.v7 * w.s) + (o.v5 * w.v2 - o.v2 * w.v5)) + ((o.v4 * w.v1 - o.v1 * w.v4) + (o.v6 * w.v3 - o.v3 * w.v6))
    return Octonion(s, v1, v2, v3, v4, v5, v6, v7, o.norm & w.norm)
end

(/)(o::Octonion, w::Octonion) = o * inv(w)

(==)(q::Octonion, w::Octonion) = (q.s == w.s) & (q.v1 == w.v1) & (q.v2 == w.v2) & (q.v3 == w.v3) &
                                 (q.v4 == w.v4) & (q.v5 == w.v5) & (q.v6 == w.v6) & (q.v7 == w.v7) # ignore .norm field

function exp(o::Octonion)
  s = o.s
  se = exp(s)
  scale = se
  th = abs_imag(o)
  if th > 0
    scale *= sin(th) / th
  end
  Octonion(se * cos(th),
            scale * o.v1,
            scale * o.v2,
            scale * o.v3,
            scale * o.v4,
            scale * o.v5,
            scale * o.v6,
            scale * o.v7,
            abs(s) == 0)
end

function log(o::Octonion)
  o, a = normalizea(o)
  s = o.s
  M = abs_imag(o)
  th = atan(M, s)
  if M > 0
    M = th / M
    return Octonion(log(a),
                     o.v1 * M,
                     o.v2 * M,
                     o.v3 * M,
                     o.v4 * M,
                     o.v5 * M,
                     o.v6 * M,
                     o.v7 * M)
  else
    return Octonion(complex(log(a), ifelse(iszero(a), zero(th), th)))
  end
end

(^)(o::Octonion, w::Octonion) = exp(w * log(o))

function sqrt(o::Octonion)
  exp(0.5 * log(o))
end

octorand() = octo(randn(), randn(), randn(), randn(), randn(), randn(), randn(), randn())

function rand(rng::AbstractRNG, ::Random.SamplerType{Octonion{T}}) where {T<:Real}
  Octonion{T}(rand(rng, T), rand(rng, T), rand(rng, T), rand(rng, T),
              rand(rng, T), rand(rng, T), rand(rng, T), rand(rng, T), false)
end

function randn(rng::AbstractRNG, ::Type{Octonion{T}}) where {T<:AbstractFloat}
  Octonion{T}(
      randn(rng, T) * INV_SQRT_EIGHT,
      randn(rng, T) * INV_SQRT_EIGHT,
      randn(rng, T) * INV_SQRT_EIGHT,
      randn(rng, T) * INV_SQRT_EIGHT,
      randn(rng, T) * INV_SQRT_EIGHT,
      randn(rng, T) * INV_SQRT_EIGHT,
      randn(rng, T) * INV_SQRT_EIGHT,
      randn(rng, T) * INV_SQRT_EIGHT,
      false,
  )
end
