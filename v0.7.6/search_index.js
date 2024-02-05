var documenterSearchIndex = {"docs":
[{"location":"examples/rotations/#Rotations-with-quaternions","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"","category":"section"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"One of the most useful application of quaternions is representation of 3D-rotations. See also Rotations.jl documentation","category":"page"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"using Quaternions\nusing LinearAlgebra","category":"page"},{"location":"examples/rotations/#Basics","page":"Rotations with quaternions","title":"Basics","text":"","category":"section"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"A 3D rotation can be represented by a unit quaternion (versor). For example, a 90° rotation around the y-axis is q = frac1sqrt2 + 0i + frac1sqrt2 j + 0k. Rotations with quaternions have the following properties:","category":"page"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"A unit quaternion (4 real numbers) is more efficient for representing a rotation than a rotation matrix (9 real numbers).\nThis results in higher computational performance in terms of time, memory usage, and accuracy.\nThe negative of a unit quaternion represents the same rotation.\nThe conjugate of a unit quaternion represents the inverse rotation.\nThe quaternion has unit length, so conjugate and multiplicative inverse is the same.\nThe set of unit quaternion leftw + ix + jy + kz in mathbbH    x y z in mathbbR right = U(1mathbbH) simeq S^3 forms a group, and the group is homomorphic to the following groups.\nSU(2) = R in mathcalM(2mathbbC)    R R^* = I is isomorphic to U(1mathbbH).\nSO(3) = R in mathcalM(3mathbbR)    R R^top = I is homomorphic to U(1mathbbH), and the mapping U(1mathbbH) to SO(3) is double covering.","category":"page"},{"location":"examples/rotations/#Rotation-around-a-vector","page":"Rotations with quaternions","title":"Rotation around a vector","text":"","category":"section"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"A theta rotation around a unit vector v = (v_x v_y v_z) can be obtained as","category":"page"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"q = cos(theta2) + sin(theta2)(iv_x + jv_y + kv_z)","category":"page"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"function quat_from_axisangle(axis::AbstractVector, theta::Real)\n    if length(axis) != 3\n        error(\"Must be a 3-vector\")\n    end\n    s, c = sincos(theta / 2)\n    axis = normalize(axis)\n    return Quaternion(c, s*axis[1], s*axis[2], s*axis[3])\nend\nnothing  # hide","category":"page"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"q1 = quat_from_axisangle([0,1,0], deg2rad(90))  # 90° rotation around y-axis\nq2 = quat_from_axisangle([1,1,1], deg2rad(120))\nq3 = -q2  # additive inverse quaternion represents the same rotation","category":"page"},{"location":"examples/rotations/#Rotate-a-vector-with-a-quaternion","page":"Rotations with quaternions","title":"Rotate a vector with a quaternion","text":"","category":"section"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"A vector u = (u_x u_y u_z) can be rotated by a unit quaternion q. The rotated vector v = (v_x v_y v_z) can be obtained as","category":"page"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"beginaligned\nq_u = iu_x + ju_y + ku_z \nq_v = q q_u barq = 0 + iv_x + jv_y + kv_z \nv = (v_x v_y v_z)\nendaligned","category":"page"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"function rotate_vector(q::Quaternion, u::AbstractVector)\n    if length(u) != 3\n        error(\"Must be a 3-vector\")\n    end\n    q_u = Quaternion(0, u[1], u[2], u[3])\n    q_v = q*q_u*conj(q)\n    return [imag_part(q_v)...]\nend\nnothing  # hide","category":"page"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"rotate_vector(q1, [1,2,3])\nrotate_vector(q2, [1,2,3])\nrotate_vector(q3, [1,2,3])  # Same as q2","category":"page"},{"location":"examples/rotations/#Convert-a-quaternion-to-a-rotation-matrix","page":"Rotations with quaternions","title":"Convert a quaternion to a rotation matrix","text":"","category":"section"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"A unit quaternion can be converted to a rotation matrix.","category":"page"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"function rotmatrix_from_quat(q::Quaternion)\n    sx, sy, sz = 2q.s * q.v1, 2q.s * q.v2, 2q.s * q.v3\n    xx, xy, xz = 2q.v1^2, 2q.v1 * q.v2, 2q.v1 * q.v3\n    yy, yz, zz = 2q.v2^2, 2q.v2 * q.v3, 2q.v3^2\n    r = [1 - (yy + zz)     xy - sz     xz + sy;\n            xy + sz   1 - (xx + zz)    yz - sx;\n            xz - sy      yz + sx  1 - (xx + yy)]\n    return r\nend\nnothing  # hide","category":"page"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"m1 = rotmatrix_from_quat(q1)\nm2 = rotmatrix_from_quat(q2)\nm3 = rotmatrix_from_quat(q3)  # Same as q2","category":"page"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"This function does not return StaticMatrix, so the implementation is not much effective. If you need more performance, please consider using Rotations.jl.","category":"page"},{"location":"examples/rotations/#Convert-a-rotation-matrix-to-a-quaternion","page":"Rotations with quaternions","title":"Convert a rotation matrix to a quaternion","text":"","category":"section"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"A rotation matrix can be converted to a unit quaternion. The following implementation is based on https://arxiv.org/pdf/math/0701759.pdf. Note that the following mapping SO(3) to U(1mathbbH) is not surjective.","category":"page"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"function quat_from_rotmatrix(dcm::AbstractMatrix{T}) where {T<:Real}\n    a2 = 1 + dcm[1,1] + dcm[2,2] + dcm[3,3]\n    a = sqrt(a2)/2\n    b,c,d = (dcm[3,2]-dcm[2,3])/4a, (dcm[1,3]-dcm[3,1])/4a, (dcm[2,1]-dcm[1,2])/4a\n    return Quaternion(a,b,c,d)\nend\nnothing  # hide","category":"page"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"quat_from_rotmatrix(m1)\nquat_from_rotmatrix(m2)\nquat_from_rotmatrix(m3)\nquat_from_rotmatrix(m1) ≈ q1\nquat_from_rotmatrix(m2) ≈ q2\nquat_from_rotmatrix(m3) ≈ q3  # q2 == -q3","category":"page"},{"location":"examples/rotations/#Interpolate-two-rotations-(slerp)","page":"Rotations with quaternions","title":"Interpolate two rotations (slerp)","text":"","category":"section"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"Slerp (spherical linear interpolation) is a method to interpolate between two unit quaternions. This function slerp equates antipodal points, and interpolates the shortest path. Therefore, the output slerp(q1, q2, 1) may be different from q2. (slerp(q1, q2, 0) is always equal to q1.)","category":"page"},{"location":"examples/rotations/","page":"Rotations with quaternions","title":"Rotations with quaternions","text":"slerp(q1, q2, 0) ≈ q1\nslerp(q1, q2, 1) ≈ q2\nslerp(q1, q3, 1) ≈ q3\nslerp(q1, q3, 1) ≈ -q3\nr = slerp(q1, q2, 1/2)\nabs(q1-r) ≈ abs(q2-r)  # Same distance\nabs(r)  # Interpolates on the unit sphere S³","category":"page"},{"location":"api/#API","page":"APIs","title":"API","text":"","category":"section"},{"location":"api/","page":"APIs","title":"APIs","text":"Quaternion","category":"page"},{"location":"api/#Quaternions.Quaternion","page":"APIs","title":"Quaternions.Quaternion","text":"Quaternion{T<:Real} <: Number\n\nQuaternion number type with real and imaginary parts of type T.\n\nQuaternionF16, QuaternionF32, and QuaternionF64 are aliases for Quaternion{Float16}, Quaternion{Float32}, and Quaternion{Float64}, respectively.\n\nSee also: quat, real, imag_part.\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"APIs","title":"APIs","text":"quat","category":"page"},{"location":"api/#Quaternions.quat","page":"APIs","title":"Quaternions.quat","text":"quat(w, [x, y, z])\n\nConvert real numbers or arrays to quaternion. x, y, z defaults to zero.\n\nExamples\n\njulia> quat(7)\nQuaternion{Int64}(7, 0, 0, 0)\n\njulia> quat(1.0, 2, 3, 4)\nQuaternionF64(1.0, 2.0, 3.0, 4.0)\n\njulia> quat([1, 2, 3])\n3-element Vector{Quaternion{Int64}}:\n Quaternion{Int64}(1, 0, 0, 0)\n Quaternion{Int64}(2, 0, 0, 0)\n Quaternion{Int64}(3, 0, 0, 0)\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"APIs","title":"APIs","text":"real(::Quaternion)\nreal(::AbstractArray{<:Quaternion})\nreal(::Type{Quaternion{T}}) where {T}","category":"page"},{"location":"api/#Base.real-Tuple{Quaternion}","page":"APIs","title":"Base.real","text":"real(q::Quaternion)\n\nReturn the real part of the quaternion q.\n\nSee also: imag_part, quat\n\nExamples\n\njulia> real(quat(1,2,3,4))\n1\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.real-Tuple{AbstractArray{<:Quaternion}}","page":"APIs","title":"Base.real","text":"real(A::AbstractArray{<:Quaternion})\n\nReturn an array containing the real part of each quaternion in A.\n\nExamples\n\njulia> real([quat(5,6,7,8), 9])\n2-element Vector{Int64}:\n 5\n 9\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.real-Union{Tuple{Type{Quaternion{T}}}, Tuple{T}} where T","page":"APIs","title":"Base.real","text":"real(T::Type{<:Quaternion})\n\nReturn the type that represents the real part of a value of type T. e.g: for T == Quaternion{R}, returns R. Equivalent to typeof(real(zero(T))).\n\nExamples\n\njulia> real(Quaternion{Int})\nInt64\n\n\n\n\n\n","category":"method"},{"location":"api/","page":"APIs","title":"APIs","text":"imag_part","category":"page"},{"location":"api/#Quaternions.imag_part","page":"APIs","title":"Quaternions.imag_part","text":"imag_part(q::Quaternion{T}) -> NTuple{3, T}\n\nReturn the imaginary part of the quaternion q.\n\nNote that this function is different from Base.imag, which returns Real for complex numbers.\n\nSee also: real, conj.\n\nExamples\n\njulia> imag_part(Quaternion(1,2,3,4))\n(2, 3, 4)\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"APIs","title":"APIs","text":"round(::Quaternion)","category":"page"},{"location":"api/#Base.round-Tuple{Quaternion}","page":"APIs","title":"Base.round","text":"round(q::Quaternion[, RoundingModeReal, [RoundingModeImaginary]]; kwargs...)\nround(q::Quaternion, RoundingModeReal,\n      RoundingModeImaginary1, RoundingModeImaginary2, RoundingModeImaginary3; kwargs...)\n\nReturn the nearest integral value of the same type as the quaternion-valued q to q, breaking ties using the specified RoundingModes.\n\nThe first RoundingMode is used for rounding the real part while the second is used for rounding the imaginary parts. Alternatively, a RoundingMode may be provided for each part.\n\nThe kwargs are the same as those for round(::Real[, RoundingMode]; kwargs...).\n\nExample\n\njulia> round(quat(3.14, 4.5, 8.3, -2.8))\nQuaternionF64(3.0, 4.0, 8.0, -3.0)\n\n\n\n\n\n","category":"method"},{"location":"api/","page":"APIs","title":"APIs","text":"conj","category":"page"},{"location":"api/#Base.conj","page":"APIs","title":"Base.conj","text":"conj(q::Quaternion)\n\nCompute the quaternion conjugate of a quaternion q.\n\nExamples\n\njulia> conj(Quaternion(1,2,3,4))\nQuaternion{Int64}(1, -2, -3, -4)\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"APIs","title":"APIs","text":"inv","category":"page"},{"location":"api/#Base.inv","page":"APIs","title":"Base.inv","text":"inv(q::Quaternion)\n\nReturn the multiplicative inverse of q::Quaternion, such that q*inv(q) or inv(q)*q yields one(q) (the multiplicative identity) up to roundoff errors.\n\nExamples\n\njulia> inv(quat(1))\nQuaternionF64(1.0, -0.0, -0.0, -0.0)\n\njulia> inv(quat(1, 2, 0, 0))\nQuaternionF64(0.2, -0.4, -0.0, -0.0)\n\njulia> inv(quat(2//3))\nQuaternion{Rational{Int64}}(3//2, 0//1, 0//1, 0//1)\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"APIs","title":"APIs","text":"sign","category":"page"},{"location":"api/#Base.sign","page":"APIs","title":"Base.sign","text":"sign(q::Quaternion) -> Quaternion\n\nReturn zero if q==0 and qq otherwise.\n\nExamples\n\njulia> sign(Quaternion(4, 0, 0, 0))\nQuaternionF64(1.0, 0.0, 0.0, 0.0)\n\njulia> sign(Quaternion(1, 0, 1, 0))\nQuaternionF64(0.7071067811865475, 0.0, 0.7071067811865475, 0.0)\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"APIs","title":"APIs","text":"slerp","category":"page"},{"location":"api/#Quaternions.slerp","page":"APIs","title":"Quaternions.slerp","text":"slerp(qa::Quaternion, qb::Quaternion, t::Real)\n\nSpherical linear interpolation (Slerp) between the inputs qa and qb. Since the input is normalized inside the function, the absolute value of the return value will be 1.\n\nExamples\n\njulia> using Quaternions\n\njulia> qa = Quaternion(1,0,0,0)\nQuaternion{Int64}(1, 0, 0, 0)\n\njulia> qb = Quaternion(0,1,0,0)\nQuaternion{Int64}(0, 1, 0, 0)\n\njulia> slerp(qa, qb, 0.6)\nQuaternionF64(0.5877852522924731, 0.8090169943749475, 0.0, 0.0)\n\njulia> ans ≈ Quaternion(cospi(0.3), sinpi(0.3), 0, 0)\ntrue\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"APIs","title":"APIs","text":"Quaternions.extend_analytic","category":"page"},{"location":"api/#Quaternions.extend_analytic","page":"APIs","title":"Quaternions.extend_analytic","text":"extend_analytic(f, q::Quaternion)\n\nEvaluate the extension of the complex analytic function f to the quaternions at q.\n\nGiven q = s + a u, where s is the real part, u is a pure unit quaternion, and a ge 0 is the magnitude of the imaginary part of q,\n\nf(q) = Re(f(z)) + Im(f(z)) u\n\nis the extension of f to the quaternions, where z = s + a i is a complex analog to q.\n\nSee Theorem 5 of [Sudbery1970] for details.\n\n[Sudbery1970]: Sudbery (1979). Quaternionic analysis. Mathematical Proceedings of the Cambridge Philosophical Society,85, pp 199225 doi:10.1017/S030500410005563\n\n\n\n\n\n","category":"function"},{"location":"examples/type_parameter/#The-type-parameter-T-in-Quaternion{T}","page":"The type parameter T in Quaternion{T}","title":"The type parameter T in Quaternion{T}","text":"","category":"section"},{"location":"examples/type_parameter/","page":"The type parameter T in Quaternion{T}","title":"The type parameter T in Quaternion{T}","text":"The type parameter T <: Real in Quaternion{T} represents the type of real and imaginary parts of a quaternion.","category":"page"},{"location":"examples/type_parameter/#Lipschitz-quaternions","page":"The type parameter T in Quaternion{T}","title":"Lipschitz quaternions","text":"","category":"section"},{"location":"examples/type_parameter/","page":"The type parameter T in Quaternion{T}","title":"The type parameter T in Quaternion{T}","text":"By using this type parameter, some special quaternions such as Lipschitz quaternions L can be represented.","category":"page"},{"location":"examples/type_parameter/","page":"The type parameter T in Quaternion{T}","title":"The type parameter T in Quaternion{T}","text":"L = lefta+bi+cj+dk in mathbbH mid abcd in mathbbZright","category":"page"},{"location":"examples/type_parameter/","page":"The type parameter T in Quaternion{T}","title":"The type parameter T in Quaternion{T}","text":"using Quaternions","category":"page"},{"location":"examples/type_parameter/","page":"The type parameter T in Quaternion{T}","title":"The type parameter T in Quaternion{T}","text":"q1 = Quaternion{Int}(1,2,3,4)\nq2 = Quaternion{Int}(5,6,7,8)\nislipschitz(q::Quaternion) = isinteger(q.s) & isinteger(q.v1) & isinteger(q.v2) & isinteger(q.v3)\nislipschitz(q1)\nislipschitz(q2)\nislipschitz(q1 + q2)\nislipschitz(q1 * q2)\nislipschitz(q1 / q2)  # Division is not defined on L.\nq1 * q2 == q2 * q1  # non-commutative","category":"page"},{"location":"examples/type_parameter/#Hurwitz-quaternions","page":"The type parameter T in Quaternion{T}","title":"Hurwitz quaternions","text":"","category":"section"},{"location":"examples/type_parameter/","page":"The type parameter T in Quaternion{T}","title":"The type parameter T in Quaternion{T}","text":"If all coefficients of a quaternion are integers or half-integers, the quaternion is called a Hurwitz quaternion. The set of Hurwitz quaternions is defined by","category":"page"},{"location":"examples/type_parameter/","page":"The type parameter T in Quaternion{T}","title":"The type parameter T in Quaternion{T}","text":"H = lefta+bi+cj+dk in mathbbH mid abcd in mathbbZ  textor  abcd in mathbbZ + tfrac12right","category":"page"},{"location":"examples/type_parameter/","page":"The type parameter T in Quaternion{T}","title":"The type parameter T in Quaternion{T}","text":"Hurwitz quaternions can be implemented with HalfIntegers.jl package.","category":"page"},{"location":"examples/type_parameter/","page":"The type parameter T in Quaternion{T}","title":"The type parameter T in Quaternion{T}","text":"using HalfIntegers\nq1 = Quaternion{HalfInt}(1, 2, 3, 4)\nq2 = Quaternion{HalfInt}(5.5, 6.5, 7.5, 8.5)\nq3 = Quaternion{HalfInt}(1, 2, 3, 4.5)  # not Hurwitz quaternion\nishalfodd(x::Number) = isodd(twice(x))  # Should be defined in HalfIntegers.jl (HalfIntegers.jl#59)\nishurwitz(q::Quaternion) = (isinteger(q.s) & isinteger(q.v1) & isinteger(q.v2) & isinteger(q.v3)) | (ishalfodd(q.s) & ishalfodd(q.v1) & ishalfodd(q.v2) & ishalfodd(q.v3))\nishurwitz(q1)\nishurwitz(q2)\nishurwitz(q3)\nishurwitz(q1 + q2)\nishurwitz(q1 * q2)\nishurwitz(q1 / q2)  # Division is not defined on H.\nq1 * q2 == q2 * q1  # non-commucative\nabs2(q1)  # Squared norm is always an integer.\nabs2(q2)  # Squared norm is always an integer.\nabs2(q3)  # Squared norm is not an integer because `q3` is not Hurwitz quaternion.","category":"page"},{"location":"examples/type_parameter/#Biquaternions","page":"The type parameter T in Quaternion{T}","title":"Biquaternions","text":"","category":"section"},{"location":"examples/type_parameter/","page":"The type parameter T in Quaternion{T}","title":"The type parameter T in Quaternion{T}","text":"If all coefficients of a quaternion are complex numbers, the quaternion is called a Biquaternion. However, the type parameter T is restricted to <:Real, so biquaternions are not supported in this package. Note that Base.Complex has the same type parameter restriction, and bicomplex numbers are not supported in Base. See issue#79 for more discussion.","category":"page"},{"location":"examples/basics/#Basics","page":"Basics","title":"Basics","text":"","category":"section"},{"location":"examples/basics/#Basic-operations-for-quaternions","page":"Basics","title":"Basic operations for quaternions","text":"","category":"section"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"Quaternions can be defined with the Quaternion constructor or quat function. Note that the order of the arguments is w+xi+yj+zk, not xi+yj+zk+w.","category":"page"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"using Quaternions\nq1 = Quaternion(1,2,3,4)\nq2 = quat(5,6,7,8)\nq3 = quat(9)","category":"page"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"The multiplication is not commutative.","category":"page"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"q1 * q2\nq2 * q1","category":"page"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"The multiplicative inverse can be calculated with Base.inv.","category":"page"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"inv(q1)\ninv(q1) * q1","category":"page"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"The division is also not commutative.","category":"page"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"q1 / q2  # Same as `q1*inv(q2)` mathematically.\nq2 \\ q1  # Same as `inv(q2)*q1` mathematically.","category":"page"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"A conjugate of a quaternion can be calculated with Base.conj. But Base.imag(::Quaternion) is not defined because it should return three real values which is not consistent with imag(::Complex) and imag(::Real). Instead, the imag_part function can be used to obtain the imaginary part of a quaternion. See issue#61 for more discussion.","category":"page"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"conj(q1)\nimag(q1)  # Not supported.\nimag_part(q1)  # Use this instead.","category":"page"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"Unit quaternions can be obtained with sign.","category":"page"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"sign(q1)\nsign(q2)\nsign(q3)\nsign(quat(0))  # Zero-quaternion will not be normalized.","category":"page"},{"location":"examples/basics/#Quaternion-vs-quat","page":"Basics","title":"Quaternion vs quat","text":"","category":"section"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"The general rule is that quat is to Quaternion as complex is to Complex. Complex and Quaternion are both constructors so should return an object of the corresponding type, whereas quat and complex both can operate on types and arrays.","category":"page"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"using Quaternions","category":"page"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"Quaternion(1,2,3,4)\nquat(1,2,3,4)\nQuaternion(Int)  # Similar to `Complex(Int)`.\nquat(Int)  # Similar to `complex(Int)`.","category":"page"},{"location":"examples/basics/#Compatibility-with-Complex","page":"Basics","title":"Compatibility with Complex","text":"","category":"section"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"There is no natural embedding mathbbCtomathbbH. Thus, quat(w,x,0,0) is not equal to complex(w,x), i.e.","category":"page"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"mathbbC ni w+ix ne w+ix+0j+0k in mathbbH","category":"page"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"using Quaternions","category":"page"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"1 + complex(1,2)   # `Complex` is compatible with `Real`\n1 + quat(1,2,3,4)  # `Quaternion` is compatible with `Real`\n1 + complex(1,2) + quat(1,2,3,4)  # no compatibility\ncomplex(1,2) == quat(1,2,0,0)     # no compatibility\ncomplex(1) == quat(1)             # no compatibility\ncomplex(1) == 1 == quat(1)  # Both `quat(1)` and `complex(1)` are equal to `1`.","category":"page"},{"location":"examples/basics/","page":"Basics","title":"Basics","text":"See issue#62 for more discussion.","category":"page"},{"location":"examples/dual_quaternions/#Dual-quaternions","page":"Dual quaternions","title":"Dual quaternions","text":"","category":"section"},{"location":"examples/dual_quaternions/#Introduction","page":"Dual quaternions","title":"Introduction","text":"","category":"section"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"The dual quaternions are an example of \"biquaternions.\" They can be represented equivalently either as a dual number where both both the \"primal\" and \"tangent\" part are quaternions","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"d = q_0 + q_e epsilon = (s_0 + a_0 i + b_0 j + c_0 k) + (s_e + a_e i + b_e j + c_e k) epsilon","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"or as a quaternion where the scalar part and three imaginary parts are all dual numbers","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"d = s + ai + bj + ck = (s_0 + s_e epsilon) + (a_0 + a_e epsilon) i + (b_0 + b_e epsilon) j + (c_0 + c_e epsilon) k","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"Like unit quaternions can compactly representation rotations in 3D space, dual quaternions can compactly represent rigid transformations (rotation with translation).","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"Without any special glue code, we can construct a dual quaternion by composing ForwardDiff.Dual and Quaternion; this uses the second representation described above:","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"note: Note\nPreviously this package contained a specialized DualQuaternion type. This was removed in v0.6.0 because it offered nothing extra over composing ForwardDiff and Quaternions.","category":"page"},{"location":"examples/dual_quaternions/#Utility-functions","page":"Dual quaternions","title":"Utility functions","text":"","category":"section"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"First let's load the packages:","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"using Quaternions, ForwardDiff, Random","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"Then we'll create some utility types/functions:","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"const DualQuaternion{T} = Quaternion{ForwardDiff.Dual{Nothing,T,1}}\n\npurequat(p::AbstractVector) = quat(false, @views(p[begin:begin+2])...)\n\ndual(x::Real, v::Real) = ForwardDiff.Dual(x, v)\n\nfunction dualquat(_q0::Union{Real,Quaternion}, _qe::Union{Real,Quaternion})\n    q0 = quat(_q0)\n    qe = quat(_qe)\n    Quaternion(\n        dual(real(q0), real(qe)),\n        dual.(imag_part(q0), imag_part(qe))...,\n    )\nend\n\nfunction primal(d::DualQuaternion)\n    return Quaternion(\n        ForwardDiff.value(real(d)),\n        ForwardDiff.value.(imag_part(d))...,\n    )\nend\n\nfunction tangent(d::DualQuaternion)\n    return Quaternion(\n        ForwardDiff.partials(real(d), 1),\n        ForwardDiff.partials.(imag_part(d), 1)...,\n    )\nend\n\nfunction dualconj(d::DualQuaternion)\n    de = tangent(d)\n    return dualquat(conj(primal(d)), quat(-real(de), imag_part(de)...))\nend\n\nrotation_part(d::DualQuaternion) = primal(d)\n\ntranslation_part(d::DualQuaternion) = dualquat(true, conj(rotation_part(d)) * tangent(d))\n\n# first=true returns the translation performed before the rotation: R(p+t)\n# first=false returns the translation performed after the rotation: R(p)+t\nfunction translation(d::DualQuaternion; first::Bool=true)\n    v = first ? primal(d)' * tangent(d) : tangent(d) * primal(d)'\n    return collect(2 .* imag_part(v))\nend\n\nfunction transform(d::DualQuaternion, p::AbstractVector)\n    dp = dualquat(true, purequat(p))\n    dpnew = d * dp * dualconj(d)\n    pnew_parts = imag_part(tangent(dpnew))\n    pnew = similar(p, eltype(pnew_parts))\n    pnew .= pnew_parts\n    return pnew\nend\n\nfunction rotmatrix_from_quat(q::Quaternion)\n    sx, sy, sz = 2q.s * q.v1, 2q.s * q.v2, 2q.s * q.v3\n    xx, xy, xz = 2q.v1^2, 2q.v1 * q.v2, 2q.v1 * q.v3\n    yy, yz, zz = 2q.v2^2, 2q.v2 * q.v3, 2q.v3^2\n    r = [1 - (yy + zz)     xy - sz     xz + sy;\n            xy + sz   1 - (xx + zz)    yz - sx;\n            xz - sy      yz + sx  1 - (xx + yy)]\n    return r\nend\n\nfunction transformationmatrix(d::DualQuaternion)\n    R = rotmatrix_from_quat(rotation_part(d))\n    t = translation(d; first=false)\n    T = similar(R, 4, 4)\n    T[1:3, 1:3] .= R\n    T[1:3, 4] .= t\n    T[4, 1:3] .= 0\n    T[4, 4] = 1\n    return T\nend\n\nranddualquat(rng::AbstractRNG,T=Float64) = dualquat(rand(rng, Quaternion{T}), rand(rng, Quaternion{T}))\nranddualquat(T=Float64) = randdualquat(Random.GLOBAL_RNG,T)\nnothing  # hide","category":"page"},{"location":"examples/dual_quaternions/#Example:-transforming-a-point","page":"Dual quaternions","title":"Example: transforming a point","text":"","category":"section"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"Now we'll create a unit dual quaternion.","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"x = sign(randdualquat())","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"sign(q) == q / abs(q) both normalizes the primal part of the dual quaternion and makes the tangent part perpendicular to it.","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"abs(primal(x)) ≈ 1\nisapprox(real(primal(x)' * tangent(x)), 0; atol=1e-10)","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"Here's how we use dual quaternions to transform a point:","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"p = randn(3)","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"transform(x, p)","category":"page"},{"location":"examples/dual_quaternions/#Example:-homomorphism-from-unit-dual-quaternions-to-the-transformation-matrices","page":"Dual quaternions","title":"Example: homomorphism from unit dual quaternions to the transformation matrices","text":"","category":"section"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"Each unit dual quaternion can be mapped to an affine transformation matrix T. T can be used to transform a vector p like this:","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"T beginpmatrix p  1endpmatrix = beginpmatrix R  t  0^mathrmT  1endpmatrix beginpmatrix p  1endpmatrix = beginpmatrix Rp + t  1endpmatrix","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"where R is a rotation matrix, and t is a translation vector. Our helper function transformationmatrix maps from a unit dual quaternion to such an affine matrix.","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"y = sign(randdualquat())","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"X = transformationmatrix(x)\nY = transformationmatrix(y)\nXY = transformationmatrix(x*y)\nX*Y ≈ XY","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"We can check that our transformation using the unit dual quaternion gives the same result as transforming with an affine transformation matrix:","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"transform(x, p) ≈ (X * vcat(p, 1))[1:3]","category":"page"},{"location":"examples/dual_quaternions/#Example:-motion-planning","page":"Dual quaternions","title":"Example: motion planning","text":"","category":"section"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"For unit quaternions, spherical linear interpolation with slerp can be used to interpolate between two rotations with unit quaternions, which can be used to plan motion between two orientations. Similarly, we can interpolate between unit dual quaternions to plan motion between two rigid poses. Conveniently, we can do this using the exact same slerp implementation.","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"slerp(x, y, 0) ≈ x","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"slerp(x, y, 1) ≈ y","category":"page"},{"location":"examples/dual_quaternions/","page":"Dual quaternions","title":"Dual quaternions","text":"slerp(x, y, 0.3)","category":"page"},{"location":"#Quaternions.jl","page":"Home","title":"Quaternions.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A Julia package implementing quaternions.","category":"page"},{"location":"","page":"Home","title":"Home","text":"note: Documentation\nThe documentation is still work in progress. For more information, see alsoREADME in the repository\nTests in the repositoryFeel free to open pull requests and improve this document!","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"pkg> add Quaternions","category":"page"},{"location":"#First-example","page":"Home","title":"First example","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using Quaternions\nk = quat(0, 0, 0, 1)\nj = quat(0, 0, 1, 0)\ni = j*k\ni^2 == j^2 == k^2 == i*j*k == -1 # Similar to `im^2`.\n1 + i + k + j  # Compatible with arithmetic operations as a `Number`.","category":"page"}]
}
