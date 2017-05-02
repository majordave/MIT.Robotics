from sympy import *
import matplotlib.pyplot as plt


# Find position Matrix from parameters by forward kinematics
def MatPos(th0, th1, th2, d1, a1, a2, unit='d'):
    if unit == 'd':
        th0, th1, th2 = mpmath.radians(th0), mpmath.radians(th1), mpmath.radians(th2)

    th0 *= -1
    th1 *= -1
    th2 -= pi/2
    th2 *= -1
    T = Matrix([[cos(th0)*cos(th1 + th2), -sin(th1 + th2)*cos(th0), -sin(th0), (a1*cos(th1) + a2*cos(th1 + th2))*cos(th0)],
                [sin(th0)*cos(th1 + th2), -sin(th0)*sin(th1 + th2), cos(th0), (a1*cos(th1) + a2*cos(th1 + th2))*sin(th0)],
                [-sin(th1 + th2), -cos(th1 + th2), 0, -a1*sin(th1) - a2*sin(th1 + th2) + d1],
                [0, 0, 0, 1]])
    return T


# Find th0 angle from desired position Matrix by inverse kinematics
def Th0(pos):
    nx, ox, ax, px = pos[0,:]
    ny, oy, ay, py = pos[1,:]

    if nx != 0:
        th =  -atan(ny/nx)
    elif ox != 0:
        th = -atan(oy/ox)
    elif ax != 0:
        th =  acot(ay/ax)
    elif px != 0:
        th = -atan(py/px)
    else:
        th = atan(ax/ay)

    if th < 0:
        th += pi

    return th

# Find th1 angle from desired position Matrix and th0 value by inverse
# kinematics
def Th1(pos, th0, d1, a2):
    th0 *= -1
    nx, ny, nz = pos[:3,0]
    px, py, pz = pos[:3,-1]
    th = -2.0*atan((-713.0*nx*tan(0.5*th0)**2 + 713.0*nx + 1426.0*ny*tan(0.5*th0) + 20.0*px*tan(0.5*th0)**2 - 20.0*px - 40.0*py*tan(0.5*th0) + sqrt(400.0*d1**2*tan(0.5*th0)**4 + 800.0*d1**2*tan(0.5*th0)**2 + 400.0*d1**2 + 28520.0*d1*nz*tan(0.5*th0)**4 + 57040.0*d1*nz*tan(0.5*th0)**2 + 28520.0*d1*nz - 800.0*d1*pz*tan(0.5*th0)**4 - 1600.0*d1*pz*tan(0.5*th0)**2 - 800.0*d1*pz + 508369.0*nx**2*tan(0.5*th0)**4 - 1016738.0*nx**2*tan(0.5*th0)**2 + 508369.0*nx**2 - 2033476.0*nx*ny*tan(0.5*th0)**3 + 2033476.0*nx*ny*tan(0.5*th0) - 28520.0*nx*px*tan(0.5*th0)**4 + 57040.0*nx*px*tan(0.5*th0)**2 - 28520.0*nx*px + 57040.0*nx*py*tan(0.5*th0)**3 - 57040.0*nx*py*tan(0.5*th0) + 2033476.0*ny**2*tan(0.5*th0)**2 + 57040.0*ny*px*tan(0.5*th0)**3 - 57040.0*ny*px*tan(0.5*th0) - 114080.0*ny*py*tan(0.5*th0)**2 + 508369.0*nz**2*tan(0.5*th0)**4 + 1016738.0*nz**2*tan(0.5*th0)**2 + 508369.0*nz**2 - 28520.0*nz*pz*tan(0.5*th0)**4 - 57040.0*nz*pz*tan(0.5*th0)**2 - 28520.0*nz*pz + 400.0*px**2*tan(0.5*th0)**4 - 800.0*px**2*tan(0.5*th0)**2 + 400.0*px**2 - 1600.0*px*py*tan(0.5*th0)**3 + 1600.0*px*py*tan(0.5*th0) + 1600.0*py**2*tan(0.5*th0)**2 + 400.0*pz**2*tan(0.5*th0)**4 + 800.0*pz**2*tan(0.5*th0)**2 + 400.0*pz**2))*cos(0.5*th0)**2/(20.0*d1 + 713.0*nz - 20.0*pz))
    return th

# Find th2 angle from desired position Matrix, th0 and th1 values by inverse
# kinematics
def Th2(pos, th0, th1, d1, a2):
    th0 *= -1
    th1 *= -1
    nx, ny, nz = pos[0:3,0]
    th = asin(sin(th1)*(nx*cos(th0)+ny*sin(th0))+nz*cos(th1))+pi/2
    return th

# Validates position Matrix derived by forward kinematics, plotting res in
# each predicted position
def validate(angles, pos, d1, a1, a2, unit='d', show=False):
    rows = len(angles[:,0])
    res = zeros(rows, 3)
    error = zeros(rows, 3)
    mats = Matrix()
    for i in range(0, rows):
        th0, th1, th2 = angles[i,:]
        mat = MatPos(th0, th1, th2, d1, a1, a2, unit)
        mats  = Matrix([mats, [mat]])
        resy = mat[:3, -1]
        err1, err2, err3 = pos[i,:] - resy.T
        error[i,:] = [[abs(err1), abs(err2), abs(err3)]]
        res[i,:] = resy.T

    xAxis = range(1, rows+1)
    avg = []
    for i in xAxis:
        ex, ey, ez = error[i-1,:]
        avg.append((ex + ey + ez)/3)

    px0, = plt.plot(xAxis, pos[:,0], label=r'real $p_x$')
    py0, = plt.plot(xAxis, pos[:,1], label=r'real $p_y$')
    pz0, = plt.plot(xAxis, pos[:,2], label=r'real $p_z$')
    px1, = plt.plot(xAxis, res[:,0], label=r'calculated $p_x$')
    py1, = plt.plot(xAxis, res[:,1], label=r'calculated $p_y$')
    pz1, = plt.plot(xAxis, res[:,2], label=r'calculated $p_z$')
    plt.legend(handles=[px0, py0, pz0, px1, py1, pz1])
    plt.xlabel("Sample")
    plt.ylabel("Position [mm]")
    plt.figure()

    [ex], [ey] = plt.plot(xAxis, error[:,0], label=r'$p_x$'), plt.plot(xAxis, error[:,1], label=r'$p_y$')
    [ez], [avg] = plt.plot(xAxis, error[:,2], label=r'$p_z$'), plt.plot(xAxis, avg, label='average')
    plt.legend(handles=[ex, ey, ez, avg])
    plt.xlabel("Sample")
    plt.ylabel("Error [mm]")

    if show: plt.show()
    return mats

# Validates rotation angles derived by inverse kinematics, plotting the error in
# each predicted angle
def valinv(mats, angles, d1, a2, unit='d'):
    rows = len(mats[:,0])
    error = zeros(rows, 3)
    res = zeros(rows, 3)
    for i in range(0, rows):
        mat, = mats[i,:]
        th0 = Th0(mat)
        th1 = Th1(mat, th0, d1, a2)
        th2 = Th2(mat, th0, th1, d1, a2)
        res[i,:] = [[th0, th1, th2]]
        a0, a1, a2 = angles[i, :]

        if unit == "d":
            a0, a1, a2 = mpmath.radians(a0), mpmath.radians(a1), mpmath.radians(a2)
            angles[i, :] = [[a0, a1, a2]]

        errow = [abs(a0 - th0), abs(a1 - th1), abs(a2 - th2)]
        error[i,:] = [errow]
    
    xAxis = range(1, rows + 1)
    
    plt.figure()
    rth0, = plt.plot(xAxis, angles[:,0], label=r'real $\theta_0$')
    rth1, = plt.plot(xAxis, angles[:,1], label=r'real $\theta_1$')
    rth2, = plt.plot(xAxis, angles[:,2], label=r'real $\theta_2$')
    cth0, = plt.plot(xAxis, res[:,0], label=r'calculated $\theta_0$')
    cth1, = plt.plot(xAxis, res[:,1], label=r'calculated $\theta_1$')
    cth2, = plt.plot(xAxis, res[:,2], label=r'calculated $\theta_2$')
    plt.legend(handles=[rth0, rth1, rth2, cth0, cth1, cth2])
    plt.xlabel("Sample")
    plt.ylabel("Angle [rad]")
    
    plt.figure()
    et0, = plt.plot(xAxis, error[:,0], label=r'$\theta_0$')
    et1, = plt.plot(xAxis, error[:,1], label=r'$\theta_1$')
    et2, = plt.plot(xAxis, error[:,2], label=r'$\theta_2$')
    plt.legend(handles=[et0, et1, et2])
    plt.xlabel("Sample")
    plt.ylabel("Error [rad]")
    plt.show()