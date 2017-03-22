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
    nx, ox, ax, px = pos.row(0)
    ny, oy, ay, py = pos.row(1)

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

# Validates position sym.Matrix derived by forward kinematics, plotting res in
# each predicted position
def validate(angles, pos, unit, d1, a1, a2):
    res = sym.Matrix()
    error = sym.Matrix()
    mats = sym.Matrix()
    for i in range(0, len(angles.col(0))):
        th0, th1, th2 = angles.row(i)
        mat = MatPos(th0, th1, th2, unit, d1, a1, a2)
        mats  = sym.Matrix([mats, [mat]])
        resy = mat[:3, -1]
        err1, err2, err3 = pos.row(i) - resy.T
        error = sym.Matrix([error, [abs(err1), abs(err2), abs(err3)]])
        res = sym.Matrix([res, resy.T])

    xAxis = range(1, len(res.col(0))+1)
    avg = []
    for i in xAxis:
        ex, ey, ez = error.row(i-1)
        avg.append((ex + ey + ez)/3)

    px0, = plt.plot(xAxis, pos.col(0), label='px0')
    py0, = plt.plot(xAxis, pos.col(1), label='py0')
    pz0, = plt.plot(xAxis, pos.col(2), label='pz0')
    px1, = plt.plot(xAxis, res.col(0), label='px1')
    py1, = plt.plot(xAxis, res.col(1), label='py1')
    pz1, = plt.plot(xAxis, res.col(2), label='pz1')
    plt.legend(handles=[px0, py0, pz0, px1, py1, pz1])
    plt.xlabel("Sample")
    plt.ylabel("Position [mm]")
    plt.figure()

    ex, = plt.plot(xAxis, error.col(0), label='ex')
    ey, = plt.plot(xAxis, error.col(1), label='ey')
    ez, = plt.plot(xAxis, error.col(2), label='ez')
    avg, = plt.plot(xAxis, avg, label='average')
    plt.legend(handles=[ex, ey, ez, avg])
    plt.xlabel("Sample")
    plt.ylabel("Error [mm]")
    plt.ylim(0, 50)
    plt.show()
    return mats

# Validates rotation angles derived by inverse kinematics, plotting the error in
# each predicted angle
def valinv(mats, angles, d1, a2):
    error = sym.Matrix()
    for i in range(0, len(mats.col(0))):
        mat, = mats.row(i)
        th0 = Th0(mat)
        th1 = Th1(mat, th0, d1, a2)
        th2 = Th2(mat, th1, d1, a2)
        r0, r1, r2 = angles.row(i)
        errow = [abs(math.radians(r0) - th0), abs(math.radians(r1) - th1), abs(math.radians(r2) - th2)]
        error = sym.Matrix([error, errow])
    pt0, = plt.plot(range(1, len(error.col(0))+1), error.col(0), label='th0')
    pt1, = plt.plot(range(1, len(error.col(0))+1), error.col(1), label='th1')
    pt2, = plt.plot(range(1, len(error.col(0))+1), error.col(2), label='th2')
    plt.legend(handles=[pt0, pt1, pt2])
    plt.xlabel("Sample")
    plt.ylabel("Error [rads]")
    plt.ylim(0, 3.1459)
    plt.show()