from sympy import Matrix, cos, sin, tan, sqrt, asin, atan, acot, pi
from matplotlib.pyplot import plot, show, legend, xlabel, ylabel, figure, ylim
from math import radians


# Find position matrix from parameters by forward kinematics
def MatPos(th0, th1, th2, d1, a1, a2):
    th1-= pi/2
    T = Matrix([[cos(th0)*cos(th1 + th2), -sin(th1 + th2)*cos(th0),
                 -sin(th0), (a1*cos(th1) + a2*cos(th1 + th2))*cos(th0)],
                [sin(th0)*cos(th1 + th2), -sin(th0)*sin(th1 + th2),
                 cos(th0), (a1*cos(th1) + a2*cos(th1 + th2))*sin(th0)],
                [-sin(th1 + th2), -cos(th1 + th2), 0,
                 -a1*sin(th1) - a2*sin(th1 + th2) + d1],
                [0, 0, 0, 1]])
    return T


# Find th0 angle from desired position matrix by inverse kinematics
def Th0(pos):
    nx, ox, ax, px = pos.row(0)
    ny, oy, ay, py = pos.row(1)

    if nx != 0:
        th =  atan(ny/nx)
    elif ox != 0:
        th = atan(oy/ox)
    elif ax != 0:
        th =  -acot(ay/ax)
    elif px != 0:
        th = atan(py/px)
    else:
        th = -atan(ax/ay)

    if th < 0:
        th += pi
    
    return th

# Find th1 angle from desired position matrix and th0 value by inverse
# kinematics
def Th1(pos, th0, d1, a2):
    nx, ny, nz = pos[:3,0]
    px, py, pz = pos[:3,-1]
    th1 = 2.0*atan((-713.0*nx*tan(0.5*th0)**2 + 713.0*nx + 1426.0*ny*tan(0.5*th0) + 20.0*px*tan(0.5*th0)**2 - 20.0*px - 40.0*py*tan(0.5*th0) + sqrt(400.0*d1**2*tan(0.5*th0)**4 + 800.0*d1**2*tan(0.5*th0)**2 + 400.0*d1**2 + 28520.0*d1*nz*tan(0.5*th0)**4 + 57040.0*d1*nz*tan(0.5*th0)**2 + 28520.0*d1*nz - 800.0*d1*pz*tan(0.5*th0)**4 - 1600.0*d1*pz*tan(0.5*th0)**2 - 800.0*d1*pz + 508369.0*nx**2*tan(0.5*th0)**4 - 1016738.0*nx**2*tan(0.5*th0)**2 + 508369.0*nx**2 - 2033476.0*nx*ny*tan(0.5*th0)**3 + 2033476.0*nx*ny*tan(0.5*th0) - 28520.0*nx*px*tan(0.5*th0)**4 + 57040.0*nx*px*tan(0.5*th0)**2 - 28520.0*nx*px + 57040.0*nx*py*tan(0.5*th0)**3 - 57040.0*nx*py*tan(0.5*th0) + 2033476.0*ny**2*tan(0.5*th0)**2 + 57040.0*ny*px*tan(0.5*th0)**3 - 57040.0*ny*px*tan(0.5*th0) - 114080.0*ny*py*tan(0.5*th0)**2 + 508369.0*nz**2*tan(0.5*th0)**4 + 1016738.0*nz**2*tan(0.5*th0)**2 + 508369.0*nz**2 - 28520.0*nz*pz*tan(0.5*th0)**4 - 57040.0*nz*pz*tan(0.5*th0)**2 - 28520.0*nz*pz + 400.0*px**2*tan(0.5*th0)**4 - 800.0*px**2*tan(0.5*th0)**2 + 400.0*px**2 - 1600.0*px*py*tan(0.5*th0)**3 + 1600.0*px*py*tan(0.5*th0) + 1600.0*py**2*tan(0.5*th0)**2 + 400.0*pz**2*tan(0.5*th0)**4 + 800.0*pz**2*tan(0.5*th0)**2 + 400.0*pz**2))*cos(0.5*th0)**2/(20.0*d1 + 713.0*nz - 20.0*pz))
    return th1+pi/2
    
# Find th2 angle from desired position matrix, th0 and th1 values by inverse
# kinematics
def Th2(pos, th1, d1, a2):
    nz = pos[2,0]
    return -th1 - asin(nz) + pi/2

# Validates position matrix derived by forward kinematics, plotting res in
# each predicted position
def validate(angles, pos, unit, d1, a1, a2):
    res = Matrix()
    error = Matrix()
    mats = Matrix()
    for i in range(0, len(angles.col(0))):
        th0, th1, th2 = angles.row(i)
        if unit == 'd':
            th0, th1, th2 = radians(th0), radians(th1), radians(th2)
        mat = MatPos(th0, th1, th2, d1, a1, a2)
        mats  = Matrix([mats, [mat]])
        resy = mat[:3, -1]
        err1, err2, err3 = pos.row(i) - resy.T
        error = Matrix([error, [abs(err1), abs(err2), abs(err3)]])
        res = Matrix([res, resy.T])

    xAxis = range(1, len(res.col(0))+1)
    avg = []
    for i in xAxis:
        ex, ey, ez = error.row(i-1)
        avg.append((ex + ey + ez)/3)

    px0, = plot(xAxis, pos.col(0), label='px0')
    py0, = plot(xAxis, pos.col(1), label='py0')
    pz0, = plot(xAxis, pos.col(2), label='pz0')
    px1, = plot(xAxis, res.col(0), label='px1')
    py1, = plot(xAxis, res.col(1), label='py1')
    pz1, = plot(xAxis, res.col(2), label='pz1')
    legend(handles=[px0, py0, pz0, px1, py1, pz1])
    xlabel("Sample")
    ylabel("Position [mm]")
    figure()

    ex, = plot(xAxis, error.col(0), label='ex')
    ey, = plot(xAxis, error.col(1), label='ey')
    ez, = plot(xAxis, error.col(2), label='ez')
    avg, = plot(xAxis, avg, label='average')
    legend(handles=[ex, ey, ez, avg])
    xlabel("Sample")
    ylabel("Error [mm]")
    ylim(0, 50)
    show()
    return mats

# Validates rotation angles derived by inverse kinematics, plotting the error in
# each predicted angle
def valinv(mats, angles, d1, a2):
    error = Matrix()
    for i in range(0, len(mats.col(0))):
        mat, = mats.row(i)
        th0 = Th0(mat)
        th1 = Th1(mat, th0, d1, a2)
        th2 = Th2(mat, th1, d1, a2)
        r0, r1, r2 = angles.row(i)
        errow = [abs(radians(r0) - th0), abs(radians(r1) - th1), abs(radians(r2) - th2)]
        error = Matrix([error, errow])
    pt0, = plot(range(1, len(error.col(0))+1), error.col(0), label='th0')
    pt1, = plot(range(1, len(error.col(0))+1), error.col(1), label='th1')
    pt2, = plot(range(1, len(error.col(0))+1), error.col(2), label='th2')
    legend(handles=[pt0, pt1, pt2])
    xlabel("Sample")
    ylabel("Error [rads]")
    ylim(0, 3.1459)
    show()