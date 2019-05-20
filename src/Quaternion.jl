struct Quaternion{T<:Real} <: Number
    s::T
    v1::T
    v2::T
    v3::T
    norm::Bool
end

(::Type{Quaternion{T}})(x::Real) where {T} = Quaternion{T}(x, 0, 0, 0, false)
(::Type{Quaternion{T}})(q::Quaternion{T}) where {T<:Real} = q
(::Type{Quaternion{T}})(q::Quaternion) where {T<:Real} = Quaternion{T}(q.s, q.v1, q.v2, q.v3, q.norm)
Quaternion(s::Real, v1::Real, v2::Real, v3::Real, n::Bool = false) =
    Quaternion(promote(s, v1, v2, v3)..., n)
Quaternion(x::Real) = Quaternion(x, zero(x), zero(x), zero(x), abs(x) == one(x))
Quaternion(z::Complex) = Quaternion(z.re, z.im, zero(z.re), zero(z.re), abs(z) == one(z.re))
Quaternion(s::Real, a::Vector) = Quaternion(s, a[1], a[2], a[3])
Quaternion(a::Vector) = Quaternion(0, a[1], a[2], a[3])

convert(::Type{Quaternion{T}}, x::Real) where {T} =
    Quaternion(convert(T, x), convert(T, 0), convert(T, 0), convert(T, 0))
convert(::Type{Quaternion{T}}, z::Complex) where {T} =
    Quaternion(convert(T, real(z)), convert(T, imag(z)), convert(T, 0), convert(T, 0))
convert(::Type{Quaternion{T}}, q::Quaternion{T}) where {T <: Real} = q
convert(::Type{Quaternion{T}}, q::Quaternion) where {T} =
    Quaternion(convert(T, q.s), convert(T, q.v1), convert(T, q.v2), convert(T, q.v3), q.norm)

promote_rule(::Type{Quaternion{T}}, ::Type{T}) where {T <: Real} = Quaternion{T}
promote_rule(::Type{Quaternion}, ::Type{T}) where {T <: Real} = Quaternion
promote_rule(::Type{Quaternion{T}}, ::Type{S}) where {T <: Real, S <: Real} = Quaternion{promote_type(T, S)}
promote_rule(::Type{Complex{T}}, ::Type{Quaternion{S}}) where {T <: Real, S <: Real} = Quaternion{promote_type(T, S)}
promote_rule(::Type{Quaternion{T}}, ::Type{Quaternion{S}}) where {T <: Real, S <: Real} = Quaternion{promote_type(T, S)}

quat(p, v1, v2, v3, n = false) = Quaternion(p, v1, v2, v3, n)
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

isfinite(q::Quaternion) = q.norm ? true : (isfinite(q.s) && isfinite(q.v1) && isfinite(q.v2) && isfinite(q.v3))

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

function exp(q::Quaternion)
    s = q.s
    se = exp(s)
    scale = se
    th = abs_imag(q)
    if th > 0
        scale *= sin(th) / th
    end
    Quaternion(se * cos(th), scale * q.v1, scale * q.v2, scale * q.v3, abs(s) < eps(typeof(s)))
end

function log(q::Quaternion)
    q, a = normalizea(q)
    s = q.s
    M = abs_imag(q)
    th = atan(M, s)
    if M > 0
        M = th / M
        return Quaternion(log(a), q.v1 * M, q.v2 * M, q.v3 * M)
    else
        return Quaternion(log(a), th, 0.0, 0.0)
    end
end

function sin(q::Quaternion)
    L = argq(q)
    return (exp(L * q) - exp(-L * q)) / (2 * L)
end

function cos(q::Quaternion)
    L = argq(q)
    return (exp(L * q) + exp(-L * q)) / 2
end

(^)(q::Quaternion, w::Quaternion) = exp(w * log(q))

sqrt(q::Quaternion) = exp(0.5 * log(q))

function linpol(p::Quaternion, q::Quaternion, t::Real)
    p = normalize(p)
    q = normalize(q)
    qm = -q
    if abs(p - q) > abs(p - qm)
        q = qm
    end
    c = p.s * q.s + p.v1 * q.v1 + p.v2 * q.v2 + p.v3 * q.v3
    if c > - 1.0
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
    else
        s  =  p.v3
        v1 = -p.v2
        v2 =  p.v1
        v3 = -p.s
        sp = sin((0.5 - t) * pi)
        sq = sin(t * pi)
        Quaternion(s,
                   sp * p.v1 + sq * v1,
                   sp * p.v2 + sq * v2,
                   sp * p.v3 + sq * v3, true)
    end
end

quatrand()  = quat(randn(), randn(), randn(), randn())
nquatrand() = normalize(quatrand())

## Rotations

function qrotation(axis::Vector{T}, theta) where {T <: Real}
    if length(axis) != 3
        error("Must be a 3-vector")
    end
    u = normalize(axis)
    thetaT = convert(eltype(u), theta)
    s = sin(thetaT / 2)
    Quaternion(cos(thetaT / 2), s * u[1], s * u[2], s * u[3], true)
end

# Variant of the above where norm(rotvec) encodes theta
function qrotation(rotvec::Vector{T}) where {T <: Real}
    if length(rotvec) != 3
        error("Must be a 3-vector")
    end
    theta = norm(rotvec)
    if theta > 0
        s = sin(theta / 2) / theta  # divide by theta to make rotvec a unit vector
        return Quaternion(cos(theta / 2), s * rotvec[1], s * rotvec[2], s * rotvec[3], true)
    end
    Quaternion(one(T), zero(T), zero(T), zero(T), true)
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

function normalize(v::Vector{T}) where {T}
    nv = norm(v)
    if nv > 0
        return v / nv
    end
    zeros(promote_type(T, typeof(nv)), length(v))
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
