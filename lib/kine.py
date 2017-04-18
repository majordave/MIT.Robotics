import sympy as sym
import matplotlib.pyplot as plt
import math


# Find position Matrix from parameters by forward kinematics
def MatPos(th0, th1, th2, unit, d1, a1, a2):
    if unit == 'd':
            th0, th1, th2 = math.radians(th0), math.radians(th1), math.radians(th2)
    th1-= sym.pi/2
    T = sym.Matrix([[sym.cos(th0)*sym.cos(th1 + th2), -sym.sin(th1 + th2)*sym.cos(th0), -sym.sin(th0), 
                 (a1*sym.cos(th1) + a2*sym.cos(th1 + th2))*sym.cos(th0)],
                [sym.sin(th0)*sym.cos(th1 + th2), -sym.sin(th0)*sym.sin(th1 + th2), sym.cos(th0), 
                 (a1*sym.cos(th1) + a2*sym.cos(th1 + th2))*sym.sin(th0)],
                [-sym.sin(th1 + th2), -sym.cos(th1 + th2), 0,
                 -a1*sym.sin(th1) - a2*sym.sin(th1 + th2) + d1],
                [0, 0, 0, 1]])
    return T


# Find th0 angle from desired position Matrix by inverse kinematics
def Th0(pos):
    nx, ox, ax, px = pos[0,:]
    ny, oy, ay, py = pos[1,:]

    if nx != 0:
        th =  sym.atan(ny/nx)
    elif ox != 0:
        th = sym.atan(oy/ox)
    elif ax != 0:
        th =  -sym.acot(ay/ax)
    elif px != 0:
        th = sym.atan(py/px)
    else:
        th = -sym.atan(ax/ay)

    if th < 0:
        th += sym.pi
    
    return th

# Find th1 angle from desired position Matrix and th0 value by inverse
# kinematics
def Th1(pos, th0, d1, a2):
    nx, ny, nz = pos[:3,0]
    px, py, pz = pos[:3,-1]
    th1 = 2.0*sym.atan((-713.0*nx*sym.tan(0.5*th0)**2 + 713.0*nx + 1426.0*ny*sym.tan(0.5*th0) + 20.0*px*sym.tan(0.5*th0)**2 - 20.0*px - 40.0*py*sym.tan(0.5*th0) + sym.sqrt(400.0*d1**2*sym.tan(0.5*th0)**4 + 800.0*d1**2*sym.tan(0.5*th0)**2 + 400.0*d1**2 + 28520.0*d1*nz*sym.tan(0.5*th0)**4 + 57040.0*d1*nz*sym.tan(0.5*th0)**2 + 28520.0*d1*nz - 800.0*d1*pz*sym.tan(0.5*th0)**4 - 1600.0*d1*pz*sym.tan(0.5*th0)**2 - 800.0*d1*pz + 508369.0*nx**2*sym.tan(0.5*th0)**4 - 1016738.0*nx**2*sym.tan(0.5*th0)**2 + 508369.0*nx**2 - 2033476.0*nx*ny*sym.tan(0.5*th0)**3 + 2033476.0*nx*ny*sym.tan(0.5*th0) - 28520.0*nx*px*sym.tan(0.5*th0)**4 + 57040.0*nx*px*sym.tan(0.5*th0)**2 - 28520.0*nx*px + 57040.0*nx*py*sym.tan(0.5*th0)**3 - 57040.0*nx*py*sym.tan(0.5*th0) + 2033476.0*ny**2*sym.tan(0.5*th0)**2 + 57040.0*ny*px*sym.tan(0.5*th0)**3 - 57040.0*ny*px*sym.tan(0.5*th0) - 114080.0*ny*py*sym.tan(0.5*th0)**2 + 508369.0*nz**2*sym.tan(0.5*th0)**4 + 1016738.0*nz**2*sym.tan(0.5*th0)**2 + 508369.0*nz**2 - 28520.0*nz*pz*sym.tan(0.5*th0)**4 - 57040.0*nz*pz*sym.tan(0.5*th0)**2 - 28520.0*nz*pz + 400.0*px**2*sym.tan(0.5*th0)**4 - 800.0*px**2*sym.tan(0.5*th0)**2 + 400.0*px**2 - 1600.0*px*py*sym.tan(0.5*th0)**3 + 1600.0*px*py*sym.tan(0.5*th0) + 1600.0*py**2*sym.tan(0.5*th0)**2 + 400.0*pz**2*sym.tan(0.5*th0)**4 + 800.0*pz**2*sym.tan(0.5*th0)**2 + 400.0*pz**2))*sym.cos(0.5*th0)**2/(20.0*d1 + 713.0*nz - 20.0*pz))
    return th1+sym.pi/2
    
# Find th2 angle from desired position Matrix, th0 and th1 values by inverse
# kinematics
def Th2(pos, th1, d1, a2):
    nz = pos[2,0]
    return -th1 - sym.asin(nz) + sym.pi/2

# Validates position Matrix derived by forward kinematics, plotting res in
# each predicted position
def validate(angles, pos, unit, d1, a1, a2):
    res = sym.Matrix()
    error = sym.Matrix()
    mats = sym.Matrix()
    for i in range(0, len(angles[:,0])):
        th0, th1, th2 = angles[i,:]
        mat = MatPos(th0, th1, th2, unit, d1, a1, a2)
        mats  = sym.Matrix([mats, [mat]])
        resy = mat[:3, -1]
        err1, err2, err3 = pos[i,:] - resy.T
        error = sym.Matrix([error, [abs(err1), abs(err2), abs(err3)]])
        res = sym.Matrix([res, resy.T])

    xAxis = range(1, len(res[:,0])+1)
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
    plt.ylim(0, 50)
    plt.show()
    return mats

# Validates rotation angles derived by inverse kinematics, plotting the error in
# each predicted angle
def valinv(mats, angles, unit, d1, a2):
    error = sym.Matrix()
    res = sym.Matrix()
    for i in range(0, len(mats[:,0])):
        mat, = mats[i,:]
        th0 = Th0(mat)
        th1 = Th1(mat, th0, d1, a2)
        th2 = Th2(mat, th1, d1, a2)
        res = sym.Matrix([res, [th0, th1, th2]])
        a0, a1, a2 = angles[i, :]
        
        if unit == "d":
            a0, a1, a2 = math.radians(a0), math.radians(a1), math.radians(a2)            
            angles[i, :] = [[a0, a1, a2]]
        
        errow = [abs(a0 - th0), abs(a1 - th1), abs(a2 - th2)]
        error = sym.Matrix([error, errow])
    
    xAxis = range(1, len(res[:,0])+1)
    
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