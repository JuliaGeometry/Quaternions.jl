using DualNumbers

struct DualQuaternion{T<:Real} <: Number
  q0::Quaternion{T}
  qe::Quaternion{T}
  norm::Bool
end

DualQuaternion(q0::Quaternion, qe::Quaternion, n::Bool = false) =
  DualQuaternion(promote(q0, qe)..., n)

DualQuaternion(d1::Dual, d2::Dual, d3::Dual, d4::Dual, n::Bool = false) =
  DualQuaternion(Quaternion(d1.re, d2.re, d3.re, d4.re, n),
                  Quaternion(d1.du, d2.du, d3.du, d4.du), n)

DualQuaternion(x::Real) = DualQuaternion(Quaternion(x), Quaternion(zero(x)), abs(x) == one(x))

DualQuaternion(d::Dual) = DualQuaternion(d, zero(d), zero(d), zero(d), abs(d) == one(d.re))

DualQuaternion(q::Quaternion) = DualQuaternion(q, zero(q), q.norm)

DualQuaternion(a::Vector) = DualQuaternion(zero(Quaternion{typeof(a[1])}), Quaternion(a))

convert(::Type{DualQuaternion{T}}, x::Real) where {T} =
    DualQuaternion(convert(Quaternion{T}, x), convert(Quaternion{T}, 0))

convert(::Type{DualQuaternion{T}}, d::Dual) where {T} =
    DualQuaternion(convert(Dual{T}, d), convert(Dual{T}, 0), convert(Dual{T}, 0), convert(Dual{T}, 0))

convert(::Type{DualQuaternion{T}}, q::Quaternion) where {T} =
    DualQuaternion(convert(Quaternion{T}, q), convert(Quaternion{T}, 0), q.norm)

convert(::Type{DualQuaternion{T}}, q::DualQuaternion{T}) where {T <: Real} = q

convert(::Type{DualQuaternion{T}}, dq::DualQuaternion) where {T} =
    DualQuaternion(convert(Quaternion{T}, dq.q0), convert(Quaternion{T}, dq.qe), dq.norm)

promote_rule(::Type{DualQuaternion{T}}, ::Type{T}) where {T <: Real} = DualQuaternion{T}
promote_rule(::Type{DualQuaternion}, ::Type{T}) where {T <: Real} = DualQuaternion
promote_rule(::Type{DualQuaternion{T}}, ::Type{S}) where {T <: Real, S <: Real} = DualQuaternion{promote_type(T, S)}
promote_rule(::Type{Quaternion{T}}, ::Type{DualQuaternion{S}}) where {T <: Real, S <: Real} = DualQuaternion{promote_type(T, S)}
promote_rule(::Type{DualQuaternion{T}}, ::Type{DualQuaternion{S}}) where {T <: Real, S <: Real} = DualQuaternion{promote_type(T, S)}

dualquat(q1, q2, n=false) = DualQuaternion(q1, q2, n)
dualquat(d1, d2, d3, d4, n=false) = DualQuaternion(d1, d2, d3, d4, n)
dualquat(x) = DualQuaternion(x)

function show(io::IO, dq::DualQuaternion)
  show(io, dq.q0)
  print(io, " + ( ")
  show(io, dq.qe)
  print(io, " )du")
end

Q0(dq::DualQuaternion) = dq.q0
Qe(dq::DualQuaternion) = dq.qe

(/)(dq::DualQuaternion, x::Real) = DualQuaternion(dq.q0 / x, dq.qe / x)

(/)(dq::DualQuaternion, d::Dual) =
  DualQuaternion(dual(dq.q0.s, dq.qe.s) / d,
                  dual(dq.q0.v1, dq.qe.v1) / d,
                  dual(dq.q0.v2, dq.qe.v2) / d,
                  dual(dq.q0.v3, dq.qe.v3) / d)

abs2(dq::DualQuaternion) = dq.norm ? dual(one(dq.q0.s)) :
  dual(abs2(dq.q0),
        2.0 * (dq.q0.s  * dq.qe.s  +
                dq.q0.v1 * dq.qe.v1 +
                dq.q0.v2 * dq.qe.v2 +
                dq.q0.v3 * dq.qe.v3))

abs(dq::DualQuaternion) = dq.norm ? dual(one(dq.q0.s)) : sqrt(abs2(dq))

conj(dq::DualQuaternion) = DualQuaternion(conj(dq.q0), conj(dq.qe), dq.norm)
dconj(dq::DualQuaternion) = DualQuaternion(dq.q0, -dq.qe, dq.norm)

inv(dq::DualQuaternion) = dq.norm ? conj(dq) : conj(dq) / abs2(dq)

function normalize(dq::DualQuaternion)
  if (dq.norm)
    return dq
  end
  a = abs(dq)
  if abs(a) > 0
    qa = dq / a
    dualquat(qa.q0, qa.qe, true)
  else
    dq
  end
end

function normalizea(dq::DualQuaternion)
  if (dq.norm)
    return (dq, one(dual))
  end
  a = abs(dq)
  if abs(a) > 0
    qa = dq / a
    dualquat(qa.q0, qa.qe, true), a
  else
    dq, zero(dual)
  end
end

(-)(dq::DualQuaternion) = DualQuaternion(-dq.q0, -dq.qe, dq.norm)

(+)(dq::DualQuaternion, dw::DualQuaternion) = DualQuaternion(dq.q0 + dw.q0, dq.qe + dw.qe)
(-)(dq::DualQuaternion, dw::DualQuaternion) = DualQuaternion(dq.q0 - dw.q0, dq.qe - dw.qe)
(*)(dq::DualQuaternion, dw::DualQuaternion) = DualQuaternion(dq.q0 * dw.q0,
                                                              dq.q0 * dw.qe + dq.qe * dw.q0,
                                                              dq.norm && dw.norm)
(/)(dq::DualQuaternion, dw::DualQuaternion) = dq * inv(dw)

function angleaxis(dq::DualQuaternion)
  tq = dq.qe * conj(dq.q0)
  t = [2.0 * tq.v1, 2.0 * tq.v2, 2.0 * tq.v3]
  q0s = dq.q0.s
  th0, s0 = angleaxis(dq.q0)
  sq0 = quat(0.0, s0)
  if abs(abs(q0s) - one(q0s)) == 0
    th = dual(th0, 0.5 * abs(quat(0, t)))
    th, dualquat(sq0)
  else
    th = dual(th0, 0.5 * dot(t, s0))
    s0c1 = cross(s0, t)
    tanth = tan(th0)
    s0c2 = (s0c1 / tanth + t) * 0.5
    sqev = cross(s0c2, s0)
    th, dualquat(sq0, quat(0.0, sqev))
  end
end

function angle(dq::DualQuaternion)
  th, ax = angleaxis(dq)
  th
end

function axis(dq::DualQuaternion)
  th, ax = angleaxis(dq)
  ax
end

function exp(dq::DualQuaternion)
  se = dual(dq.q0.s, dq.qe.s)
  se = exp(se)
  dq = dualquat(quat(0.0, imag(dq.q0)), quat(0.0, imag(dq.qe)))
  dq, th = normalizea(dq)
  if dq.norm
    dualquat(se) * (dualquat(cos(th)) + dq * dualquat(sin(th)))
  else
    dualquat(se)
  end
end

function log(dq::DualQuaternion)
  dq, a = normalizea(dq)
  sl = log(a)
  th, s = angleaxis(dq)
  s * dualquat(th) + dualquat(sl)
end

(^)(dq::DualQuaternion, dw::DualQuaternion) = exp(dw * log(dq))

function sqrt(dq::DualQuaternion)
  exp(0.5 * log(dq))
end

dualquatrand() = dualquat(quatrand(), quatrand())
ndualquatrand() = normalize(dualquatrand())
