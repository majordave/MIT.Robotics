from sympy import *
import matplotlib.pyplot as plt


# Find position matrix from parameters by forward kinematics
def MatPos(angles, d1, a1, a2):
    [th0, th1, th2] = angles
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
    [nx, ox, ax, px] = pos.row(0)
    [ny, oy, ay, py] = pos.row(1)

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
def Th1(pos, th0):
    var('th1')
    [nx, ny, nz, dum] = pos.col(0)
    [px, py, pz, dum] = pos.col(-1)
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
def Th2(pos, th0, th1):
    th = []
    [nx, ny, nz, dum] = pos.col(0)
    [px, py, pz, dum] = pos.col(-1)
    for val in th1:
        # Unused equivalent th2 equation: th.append(-asin(nz)-val)
        th.append(asin((cos(val)*(d1-pz)-sin(val)*(px*cos(th0)+py*sin(th0)))/a2))
    return th


# Find equivalent angle of inverse direction rotation
def equAngle(angle):
    angle = angle/pi*180
    equ = angle - 360
    return equ/180*pi


# Validates position matrix derived by forward kinematics, plotting error in
# each predicted position
def validate(angles, pos):
    error = Matrix()
    for i in range(0, len(angles.col(0))):
        mat = MatPos(angles.row(i), d1, a1, a2)
        res = mat.col(-1).T
        res.col_del(-1)
        [err1, err2, err3] = pos.row(i) - res
        error = Matrix([error, [abs(err1), abs(err2), abs(err3)]])
    xAxis = range(0, len(error.col(0)))
    plt.plot(xAxis, error.col(0), 'ro')
    plt.plot(xAxis, error.col(1), 'go')
    plt.plot(xAxis, error.col(2), 'bo')
    plt.show()


# Testing paramters and method calls
a1 = 10
a2 = 10
d1 = 10
'''
test = MatPos(Matrix([0, -pi/2, 0]), d1, a1, a2)
pprint(test)
th0 = Th0(test)
th1 = Th1(test, th0)
th2 = Th2(test, th0, th1)
print("th0 = ", th0)
print("th1 = ", th1)
print("th2 = ", th2)
'''
test = Matrix([[0, -pi/2, 0],
               [0, -pi/2, 0]])
pos = Matrix([[0, 0, 30],
              [0, 0, -30]])
validate(test, pos)
