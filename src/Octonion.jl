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

Octonion(s::Real, v1::Real, v2::Real, v3::Real, v4::Real, v5::Real, v6::Real, v7::Real, n::Bool = false) =
  Octonion(promote(s, v1, v2, v3, v4, v5, v6, v7)..., n)
Octonion(x::Real) = Octonion(x, zero(x), zero(x), zero(x), zero(x), zero(x), zero(x), zero(x), abs(x) == one(x))
Octonion(z::Complex) = Octonion(z.re, z.im, zero(z.re), zero(z.re), zero(z.re), zero(z.re), zero(z.re), zero(z.re), abs(z) == one(z.re))
Octonion(q::Quaternion) = Octonion(q.s, q.v1, q.v2, q.v3, zero(q.s), zero(q.s), zero(q.s), zero(q.s), q.norm)
Octonion(s::Real, a::Vector) = Octonion(s, a[1], a[2], a[3], a[4], a[5], a[6], a[7])
Octonion(a::Vector) = Octonion(0, a[1], a[2], a[3], a[4], a[5], a[6], a[7])

convert(::Type{Octonion{T}}, x::Real) where {T} =
  Octonion(convert(T, x), convert(T, 0), convert(T, 0), convert(T, 0), convert(T, 0), convert(T, 0), convert(T, 0), convert(T, 0))
convert(::Type{Octonion{T}}, z::Complex) where {T} =
  Octonion(convert(T, real(z)), convert(T, imag(z)), convert(T, 0), convert(T, 0), convert(T, 0), convert(T, 0), convert(T, 0), convert(T, 0))
convert(::Type{Octonion{T}}, q::Quaternion) where {T} =
  Octonion(convert(T, real(q)), convert(T, q.v1), convert(T, q.v2), convert(T, q.v3), convert(T, 0), convert(T, 0), convert(T, 0), convert(T, 0))
convert(::Type{Octonion{T}}, o::Octonion{T}) where {T <: Real} = o
convert(::Type{Octonion{T}}, o::Octonion) where {T} =
  Octonion(convert(T, o.s), convert(T, o.v1), convert(T, o.v2), convert(T, o.v3), convert(T, o.v4), convert(T, o.v5), convert(T, o.v6), convert(T, o.v7), o.norm)

promote_rule(::Type{Octonion{T}}, ::Type{T}) where {T <: Real} = Octonion{T}
promote_rule(::Type{Octonion}, ::Type{T}) where {T <: Real} = Octonion
promote_rule(::Type{Octonion{T}}, ::Type{S}) where {T <: Real, S <: Real} = Octonion{promote_type(T, S)}
promote_rule(::Type{Complex{T}}, ::Type{Octonion{S}}) where {T <: Real, S <: Real} = Octonion{promote_type(T, S)}
promote_rule(::Type{Quaternion{T}}, ::Type{Octonion{S}}) where {T <: Real, S <: Real} = Octonion{promote_type(T, S)}
promote_rule(::Type{Octonion{T}}, ::Type{Octonion{S}}) where {T <: Real, S <: Real} = Octonion{promote_type(T, S)}

octo(p, v1, v2, v3, v4, v5, v6, v7, n = false) = Octonion(p, v1, v2, v3, v4, v5, v6, v7, n)
octo(x) = Octonion(x)
octo(s, a) = Octonion(s, a)

function show(io::IO, o::Octonion)
  pm(x) = x < 0 ? " - $(-x)" : " + $x"
  print(io, o.s, pm(o.v1), "im", pm(o.v2), "jm", pm(o.v3), "km", pm(o.v4), "ilm", pm(o.v5), "jlm", pm(o.v6), "klm", pm(o.v7), "lm")
end

real(o::Octonion) = o.s
imag(o::Octonion) = [o.v1, o.v2, o.v3, o.v4, o.v5, o.v6, o.v7]

(/)(o::Octonion, x::Real) = Octonion(o.s / x, o.v1 / x, o.v2 / x, o.v3 / x, o.v4 / x, o.v5 / x, o.v6 / x, o.v7 / x)

conj(o::Octonion) = Octonion(o.s, -o.v1, -o.v2, -o.v3, -o.v4, -o.v5, -o.v6, -o.v7, o.norm)
abs(o::Octonion) = sqrt(o.s * o.s + o.v1 * o.v1 + o.v2 * o.v2 + o.v3 * o.v3 + o.v4 * o.v4 + o.v5 * o.v5 + o.v6 * o.v6 + o.v7 * o.v7)
abs2(o::Octonion) = o.s * o.s + o.v1 * o.v1 + o.v2 * o.v2  + o.v3 * o.v3 + o.v4 * o.v4 + o.v5 * o.v5 + o.v6 * o.v6 + o.v7 * o.v7
inv(o::Octonion) = o.norm ? conj(o) : conj(o) / abs2(o)

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


macro ocmul(i, j)
  si = Symbol("v$i")
  sj = Symbol("v$j")
  esc(:(o.$si * w.$sj - o.$sj * w.$si))
end

(*)(o::Octonion, w::Octonion) =
  Octonion(
            o.s * w.s - o.v1 * w.v1 - o.v2 * w.v2 - o.v3 * w.v3 - o.v4 * w.v4 - o.v5 * w.v5 - o.v6 * w.v6 - o.v7 * w.v7,
            o.s * w.v1 + o.v1 * w.s + @ocmul( 2, 3 ) + @ocmul( 6, 5 ) + @ocmul( 7, 4 ),
            o.s * w.v2 + o.v2 * w.s + @ocmul( 3, 1 ) + @ocmul( 4, 6 ) + @ocmul( 7, 5 ),
            o.s * w.v3 + o.v3 * w.s + @ocmul( 1, 2 ) + @ocmul( 5, 4 ) + @ocmul( 7, 6 ),
            o.s * w.v4 + o.v4 * w.s + @ocmul( 1, 7 ) + @ocmul( 3, 5 ) + @ocmul( 6, 2 ),
            o.s * w.v5 + o.v5 * w.s + @ocmul( 1, 6 ) + @ocmul( 2, 7 ) + @ocmul( 4, 3 ),
            o.s * w.v6 + o.v6 * w.s + @ocmul( 2, 4 ) + @ocmul( 3, 7 ) + @ocmul( 5, 1 ),
            o.s * w.v7 + o.v7 * w.s + @ocmul( 4, 1 ) + @ocmul( 5, 2 ) + @ocmul( 6, 3 )
          )

(/)(o::Octonion, w::Octonion) = o * inv(w)

function exp(o::Octonion)
  s = o.s
  se = exp(s)
  scale = se
  th = abs(Octonion(imag(o)))
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
  M = abs(Octonion(imag(o)))
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
    return Octonion(log(a), th, 0.0, 0.0)
  end
end

(^)(o::Octonion, w::Octonion) = exp(w * log(o))

function sqrt(o::Octonion)
  exp(0.5 * log(o))
end

octorand() = octo(randn(), randn(), randn(), randn(), randn(), randn(), randn(), randn())
