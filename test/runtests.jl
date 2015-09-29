using Base.Test
using Quaternions

qx = qrotation([1,0,0], pi/4)
@test_approx_eq qx*qx qrotation([1,0,0], pi/2)
@test_approx_eq qx^2 qrotation([1,0,0], pi/2)
theta = pi/8
qx = qrotation([1,0,0], theta)
c = cos(theta); s = sin(theta)
Rx = [1 0 0; 0 c -s; 0 s c]
@test_approx_eq rotationmatrix(qx) Rx
theta = pi/6
qy = qrotation([0,1,0], theta)
c = cos(theta); s = sin(theta)
Ry = [c 0 s; 0 1 0; -s 0 c]
@test_approx_eq rotationmatrix(qy) Ry
theta = 4pi/3
qz = qrotation([0,0,1], theta)
c = cos(theta); s = sin(theta)
Rz = [c -s 0; s c 0; 0 0 1]
@test_approx_eq rotationmatrix(qz) Rz

@test_approx_eq rotationmatrix(qx*qy*qz) Rx*Ry*Rz
@test_approx_eq rotationmatrix(qy*qx*qz) Ry*Rx*Rz
@test_approx_eq rotationmatrix(qz*qx*qy) Rz*Rx*Ry

a, b = qrotation([0,0,1], deg2rad(0)), qrotation([0,0,1], deg2rad(180))
@test_approx_eq slerp(a,b,0.0) a
@test_approx_eq slerp(a,b,1.0) b
@test_approx_eq slerp(a,b,0.5) qrotation([0,0,1], deg2rad(90))

@test_approx_eq angle(qrotation([1,0,0], 0)) 0
@test_approx_eq angle(qrotation([0,1,0], pi/4)) pi/4
@test_approx_eq angle(qrotation([0,0,1], pi/2)) pi/2
