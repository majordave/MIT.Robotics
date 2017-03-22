import sympy as sym
import matplotlib.pyplot as plt
import math

def MatJac(delta0, delta1, unit, d1, a1, a2):
    if unit == 'd':
        delta0[1] = math.radians(delta0[1])
        delta0[2] = math.radians(delta0[2])
        delta0[3] = math.radians(delta0[3])
        delta1[1] = math.radians(delta1[1])
        delta1[2] = math.radians(delta1[2])
        delta1[3] = math.radians(delta1[3])
    
    delta0[2]-= sym.pi/2 
    delta1[2]-= sym.pi/2
    th0, th1, th2 = delta1[1:]
    dtime, dth0, dth1, dth2 = delta1 - delta0
    ths = sym.Matrix([dth0/dtime, dth1/dtime, dth2/dtime])
    jacob = sym.Matrix([[sym.sin(th0)*(-a2*sym.cos(th1 + th2) - a1*sym.cos(th1)), 
                     -sym.cos(th0)*(a2*sym.sin(th1 + th2) + a1*sym.sin(th1)), 
                     -a2*sym.cos(th0)*sym.sin(th1 + th2)],
                    [sym.cos(th0)*(a2*sym.cos(th1 + th2) + a1*sym.cos(th1)), 
                     -sym.sin(th0)*(a2*sym.sin(th1 + th2) + a1*sym.sin(th1)), 
                     -a2*sym.sin(th0)*sym.sin(th1 + th2)],
                    [0, -a1*sym.cos(th1) - a2*sym.cos(th1 + th2), -a2*sym.cos(th1 + th2)],
                    [0, -sym.sin(th0), -sym.sin(th0)],
                    [0, sym.cos(th0), sym.cos(th0)],
                    [1, 0, 0]])
    return jacob*ths

def validate(angles, veloc, unit, d1, a1, a2):
    res = sym.Matrix()
    error = sym.Matrix()
    mats = sym.Matrix()
    for i in range(1, len(angles.col(0))):
        delta1 = angles.row(i)
        delta0 = angles.row(i-1)
        mat = MatJac(delta0, delta1, unit, d1, a1, a2)
        mats  = sym.Matrix([mats, [mat.T]])
        resy = mat[:3, :]
        err1, err2, err3 = veloc.row(i-1) - resy.T
        error = sym.Matrix([error, [abs(err1), abs(err2), abs(err3)]])
        res = sym.Matrix([res, resy.T])

    xAxis = range(1, len(res.col(0))+1)
    avg = []
    for i in xAxis:
        ex, ey, ez = error.row(i-1)
        avg.append((ex + ey + ez)/3)

    vx0, = plt.plot(xAxis, veloc.col(0), label='vx0')
    vy0, = plt.plot(xAxis, veloc.col(1), label='vy0')
    vz0, = plt.plot(xAxis, veloc.col(2), label='vz0')
    vx1, = plt.plot(xAxis, res.col(0), label='vx1')
    vy1, = plt.plot(xAxis, res.col(1), label='vy1')
    vz1, = plt.plot(xAxis, res.col(2), label='vz1')
    plt.legend(handles=[vx0, vy0, vz0, vx1, vy1, vz1])
    plt.xlabel("Sample")
    plt.ylabel("Velocity [mm/s]")
    plt.figure()

    ex, = plt.plot(xAxis, error.col(0), label='ex')
    ey, = plt.plot(xAxis, error.col(1), label='ey')
    ez, = plt.plot(xAxis, error.col(2), label='ez')
    avg, = plt.plot(xAxis, avg, label='average')
    plt.legend(handles=[ex, ey, ez, avg])
    plt.xlabel("Sample")
    plt.ylabel("Error [mm/s]")
    plt.ylim(0, 50)
    plt.show()
    return mats