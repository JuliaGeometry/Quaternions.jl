immutable Quaternion{T<:Real} <: Number
  s::T
  v1::T
  v2::T
  v3::T
  norm::Bool
end

Quaternion( s::Real, v1::Real, v2::Real, v3::Real, n::Bool = false) =
  Quaternion( promote( s, v1, v2, v3 )..., n )
Quaternion( x::Real ) = Quaternion( x, zero(x), zero(x), zero(x), abs(x) == one(x) )
Quaternion( z::Complex ) = Quaternion( z.re, z.im, zero(z.re), zero(z.re), abs(z) == one(z.re) )
Quaternion( s::Real, a::Vector ) = Quaternion( s, a[1], a[2], a[3] )
Quaternion( a::Vector ) = Quaternion( 0, a[1], a[2], a[3] )

convert{T}(::Type{Quaternion{T}}, x::Real) =
  Quaternion(convert(T,x), convert(T,0), convert(T,0), convert(T,0))
convert{T}(::Type{Quaternion{T}}, z::Complex) =
  Quaternion(convert(T,real(z)), convert(T,imag(z)), convert(T,0), convert(T,0))
convert{T<:Real}(::Type{Quaternion{T}}, q::Quaternion{T}) = q
convert{T}(::Type{Quaternion{T}}, q::Quaternion) =
  Quaternion(convert(T,q.s), convert(T,q.v1), convert(T,q.v2), convert(T,q.v3), q.norm)

promote_rule{T<:Real}(::Type{Quaternion{T}}, ::Type{T}) = Quaternion{T}
promote_rule{T<:Real}(::Type{Quaternion}, ::Type{T}) = Quaternion
promote_rule{T<:Real,S<:Real}(::Type{Quaternion{T}}, ::Type{S}) = Quaternion{promote_type(T,S)}
promote_rule{T<:Real,S<:Real}(::Type{Complex{T}}, ::Type{Quaternion{S}}) = Quaternion{promote_type(T,S)}
promote_rule{T<:Real,S<:Real}(::Type{Quaternion{T}}, ::Type{Quaternion{S}}) = Quaternion{promote_type(T,S)}

quat( p, v1, v2, v3, n = false ) = Quaternion( p, v1, v2, v3, n )
quat( x ) = Quaternion( x )
quat( s, a ) = Quaternion( s, a )

function show(io::IO, q::Quaternion)
  pm(x) = x < 0 ? " - $(-x)" : " + $x"
  print(io, q.s, pm(q.v1), "im", pm(q.v2), "jm", pm(q.v3), "km")
end

real( q::Quaternion ) = q.s
imag( q::Quaternion ) = [ q.v1, q.v2, q.v3 ]

(/)( q::Quaternion, x::Real ) = Quaternion( q.s/x, q.v1/x, q.v2/x, q.v3/x )

conj( q::Quaternion ) = Quaternion( q.s, -q.v1, -q.v2, -q.v3, q.norm )
abs( q::Quaternion ) = sqrt( q.s*q.s + q.v1*q.v1 + q.v2*q.v2 + q.v3*q.v3 )
abs2( q::Quaternion ) = q.s*q.s + q.v1*q.v1 + q.v2*q.v2 + q.v3*q.v3
inv( q::Quaternion ) = q.norm ? conj(q) : conj(q)/abs2(q)

function normalize( q::Quaternion )
  if ( q.norm )
    return q
  end
  q = q / abs( q )
  Quaternion( q.s, q.v1, q.v2, q.v3, true )
end

function normalizea( q::Quaternion )
  if ( q.norm )
    return ( q, one( q.s ) )
  end
  a = abs( q )
  q = q / a
  ( Quaternion( q.s, q.v1, q.v2, q.v3, true ), a )
end

function normalizeq( q::Quaternion )
  a = abs( q )
  if a > 0
    q = q / a
    Quaternion( q.s, q.v1, q.v2, q.v3, true )
  else
    Quaternion( 0.0, 1.0, 0.0, 0.0, true )
  end
end

(-)(q::Quaternion) = Quaternion( -q.s, -q.v1, -q.v2, -q.v3, q.norm )

(+)(q::Quaternion, w::Quaternion) = Quaternion( q.s + w.s,
                                                q.v1 + w.v1, q.v2 + w.v2, q.v3 + w.v3 )

(-)(q::Quaternion, w::Quaternion) = Quaternion( q.s - w.s,
                                                q.v1 - w.v1, q.v2 - w.v2, q.v3 - w.v3 )

(*)(q::Quaternion, w::Quaternion) = Quaternion( q.s*w.s - q.v1*w.v1 - q.v2*w.v2 - q.v3*w.v3,
                                                q.s*w.v1 + q.v1*w.s + q.v2*w.v3 - q.v3*w.v2,
                                                q.s*w.v2 - q.v1*w.v3 + q.v2*w.s + q.v3*w.v1,
                                                q.s*w.v3 + q.v1*w.v2 - q.v2*w.v1 + q.v3*w.s,
                                                q.norm && w.norm )
(/)( q::Quaternion, w::Quaternion ) = q * inv(w)

function angleaxis( q::Quaternion )
  q = normalize(q)
  angle = atan2( abs( Quaternion( 0, q.v1, q.v2, q.v3 ) ), q.s )
  s = sin( angle )
  axis = abs( s ) > 0 ?
    [ q.v1, q.v2, q.v3 ] / s :
    [ 1.0, 0.0, 0.0 ]
  angle, axis
end

function angle( q::Quaternion )
  q = normalize(q)
  atan2( abs( Quaternion( 0, q.v1, q.v2, q.v3 ) ), q.s )
end

function axis( q:: Quaternion )
  q = normalize(q)
  s = sin( angle( q ) )
  abs( s ) > 0 ?
    [ q.v1, q.v2, q.v3 ] / s :
    [ 1.0, 0.0, 0.0 ]
end

function argq( q::Quaternion )
  q = normalize( q )
  q = Quaternion( 0.0, q.v1, q.v2, q.v3 )
  normalizeq( q )
end

function exp( q::Quaternion )
  s = q.s
  se = exp( s )
  scale = se
  th = abs( Quaternion( imag( q ) ) )
  if th > 0
    scale *= sin( th ) / th
  end
  Quaternion( se * cos( th ), scale * q.v1, scale * q.v2, scale * q.v3, abs( s ) < eps( typeof( s ) ) )
end

function log( q::Quaternion )
  q, a = normalizea( q )
  s = q.s
  M = abs( Quaternion( imag( q ) ) )
  th = atan2( M, s )
  if M > 0
    M = th / M
    return Quaternion( log( a ), q.v1 * M, q.v2 * M, q.v3 * M )
  else
    return Quaternion( log( a ), th , 0.0, 0.0 )
  end
end

function sin( q::Quaternion )
  L = argq( q )
  ( exp( L * q ) - exp( -L * q ) ) / (2*L)
end

function cos( q::Quaternion )
  L = argq( q )
  ( exp( L * q ) + exp( -L * q ) ) / 2.0
end

(^)(q::Quaternion, w::Quaternion) = exp( w * log( q ) )

function sqrt( q::Quaternion )
  exp( 0.5 * log( q ) )
end

function linpol( p::Quaternion, q::Quaternion, t::Real )
  p = normalize( p )
  q = normalize( q )
  qm = -q
  if  abs( p - q ) > abs( p - qm )
    q = qm
  end
  c = p.s * q.s + p.v1 * q.v1 + p.v2 * q.v2 + p.v3 * q.v3
  if c > - 1.0 
    if c < 1.0
      o = acos( c )
      s = sin( o )
      sp = sin( ( 1.0 - t ) * o  ) / s
      sq = sin( t * o ) / s
    else
      sp = 1.0 - t
      sq = t
    end
    Quaternion( sp * p.s  + sq * q.s,
                sp * p.v1 + sq * q.v1,
                sp * p.v2 + sq * q.v2,
                sp * p.v3 + sq * q.v3, true )
  else
    s  =  p.v3
    v1 = -p.v2
    v2 =  p.v1
    v3 = -p.s
    sp = sin( (0.5 - t ) * pi )
    sq = sin( t * pi )
    Quaternion( s,
                sp * p.v1 + sq * v1,
                sp * p.v2 + sq * v2,
                sp * p.v3 + sq * v3, true )
  end
  
end

quatrand() = quat( rand(), rand(), rand(), rand() )
