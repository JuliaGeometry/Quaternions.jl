using DualNumbers

struct DualQuaternion{T<:Real} <: Number
  q0::Quaternion{T}
  qe::Quaternion{T}
  norm::Bool
end

DualQuaternion{T}(dq::DualQuaternion) where {T<:Real} = DualQuaternion{T}(dq.q0, dq.qe, dq.norm)
function DualQuaternion{T}(d1::Dual, d2::Dual, d3::Dual, d4::Dual, n::Bool=false) where {T<:Real}
  return DualQuaternion{T}(
    Quaternion(DualNumbers.value(d1), DualNumbers.value(d2), DualNumbers.value(d3), DualNumbers.value(d4), n),
    Quaternion(DualNumbers.epsilon(d1), DualNumbers.epsilon(d2), DualNumbers.epsilon(d3), DualNumbers.epsilon(d4)),
    n,
  )
end
function DualQuaternion{T}(q0::Quaternion) where {T<:Real}
  return DualQuaternion{T}(convert(Quaternion{T}, q0), zero(Quaternion{T}), q0.norm)
end
function DualQuaternion{T}(d::Dual) where {T<:Real}
  return DualQuaternion(
    Quaternion{T}(DualNumbers.value(d)),
    Quaternion{T}(DualNumbers.epsilon(d)),
    (DualNumbers.value(d)==one(DualNumbers.value(d))) & iszero(DualNumbers.epsilon(d)))
end
function DualQuaternion{T}(x::Real) where {T<:Real}
  return DualQuaternion(convert(Quaternion{T}, x), zero(Quaternion{T}), abs(x) == one(x))
end

DualQuaternion(q0::Quaternion, qe::Quaternion) = DualQuaternion(promote(q0, qe)..., false)
DualQuaternion(d1::Dual, d2::Dual, d3::Dual, d4::Dual, n::Bool=false) =
  DualQuaternion(Quaternion(DualNumbers.value(d1), DualNumbers.value(d2), DualNumbers.value(d3), DualNumbers.value(d4), n),
                  Quaternion(DualNumbers.epsilon(d1), DualNumbers.epsilon(d2), DualNumbers.epsilon(d3), DualNumbers.epsilon(d4)), n)
DualQuaternion(x::Real) = DualQuaternion(Quaternion(x), Quaternion(zero(x)), abs(x) == one(x))
DualQuaternion(d::Dual) = DualQuaternion(Quaternion(DualNumbers.value(d)), Quaternion(DualNumbers.epsilon(d)), (DualNumbers.value(d)==one(DualNumbers.value(d))) & iszero(DualNumbers.epsilon(d)))
DualQuaternion(q::Quaternion) = DualQuaternion(q, zero(q), q.norm)
DualQuaternion(a::Vector) = DualQuaternion(zero(Quaternion{typeof(a[1])}), Quaternion(a))

const DualQuaternionF16 = DualQuaternion{Float16}
const DualQuaternionF32 = DualQuaternion{Float32}
const DualQuaternionF64 = DualQuaternion{Float64}

promote_rule(::Type{DualQuaternion{T}}, ::Type{S}) where {T <: Real, S <: Real} = DualQuaternion{promote_type(T, S)}
promote_rule(::Type{DualQuaternion{T}}, ::Type{Dual{S}}) where {T <: Real, S <: Real} = DualQuaternion{promote_type(T, S)}
promote_rule(::Type{DualQuaternion{T}}, ::Type{Quaternion{S}}) where {T <: Real, S <: Real} = DualQuaternion{promote_type(T, S)}
promote_rule(::Type{DualQuaternion{T}}, ::Type{DualQuaternion{S}}) where {T <: Real, S <: Real} = DualQuaternion{promote_type(T, S)}

dualquat(q1, q2) = DualQuaternion(q1, q2)
dualquat(q1, q2, n) = DualQuaternion(q1, q2, n)
dualquat(d1, d2, d3, d4) = DualQuaternion(d1, d2, d3, d4)
dualquat(d1, d2, d3, d4, n) = DualQuaternion(d1, d2, d3, d4, n)
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
float(dq::DualQuaternion{T}) where T = convert(DualQuaternion{float(T)}, dq)

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
(==)(q::DualQuaternion, w::DualQuaternion) = (q.q0 == w.q0) & (q.qe == w.qe) # ignore .norm field

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

function rand(rng::AbstractRNG, ::Random.SamplerType{DualQuaternion{T}}) where {T<:Real}
    return DualQuaternion{T}(rand(rng, Quaternion{T}), rand(rng, Quaternion{T}), false)
end
