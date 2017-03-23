from sympy import Matrix, var, cos, sin, asin, atan, acot, pi, solve
from matplotlib.pyplot import plot, show, legend, xlabel, ylabel, figure
from math import radians


# Find position matrix from parameters by forward kinematics
def MatPos(th0, th1, th2, d1, a1, a2):
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
        return atan(ny/nx)
    elif ox != 0:
        return atan(oy/ox)
    elif ax != 0:
        return -acot(ay/ax)
    elif px != 0:
        return atan(py/px)
    else:
        return -atan(ax/ay)


# Find th1 angle from desired position matrix and th0 value by inverse
# kinematics
def Th1(pos, th0, d1, a2):
    th1 = var('th1')
    nx, ny, nz, dum = pos.col(0)
    px, py, pz, dum = pos.col(-1)
    eq1 = -asin(sin(th1)*(nx*cos(th0)+ny*sin(th0))+nz*cos(th1))
    eq2 = asin((cos(th1)*(d1-pz)-sin(th1)*(px*cos(th0)+py*sin(th0)))/a2)
    th = solve(eq1 - eq2, th1)
    for val in th:
        if val > pi:
            th.remove(val)
            th.append(equAngle(val))
    return th
    # Unused equivalent th1 equation: return asin((d1-pz-a2*sin(th1+th2))/a1)


# Find th2 angle from desired position matrix, th0 and th1 values by inverse
# kinematics
def Th2(pos, th0, th1, d1, a2):
    th = []
    nx, ny, nz, dum = pos.col(0)
    px, py, pz, dum = pos.col(-1)
    for val in th1:
        # Unused equivalent th2 equation: th.append(-asin(nz)-val)
        th.append(asin((cos(val)*(d1-pz)-sin(val)*(px*cos(th0)+py*sin(th0)))
                  / a2))
    return th


# Find equivalent angle of inverse direction rotation
def equAngle(angle):
    angle = angle/pi*180
    equ = angle - 360
    return equ/180*pi


# Validates position matrix derived by forward kinematics, plotting res in
# each predicted position
def validate(angles, pos, unit, d1, a1, a2):
    res = Matrix()
    error = Matrix()
    for i in range(0, len(angles.col(0))):
        th0, th1, th2 = angles.row(i)
        if unit == 'd':
            th0, th1, th2 = radians(th0), radians(th1), radians(th2)
        mat = MatPos(th0, th1, th2, d1, a1, a2)
        resy = mat.col(-1)
        resy.row_del(-1)
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
    xlabel("iteration")
    ylabel("position[mm]")
    figure()

    ex, = plot(xAxis, error.col(0), label='ex')
    ey, = plot(xAxis, error.col(1), label='ey')
    ez, = plot(xAxis, error.col(2), label='ez')
    avg, = plot(xAxis, avg, label='average')
    legend(handles=[ex, ey, ez, avg])
    xlabel("iteration")
    ylabel("error [mm]")
    show()
