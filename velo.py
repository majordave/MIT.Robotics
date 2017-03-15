from sympy import Matrix, sin, cos
from math import radians
from matplotlib.pyplot import plot, show, legend, xlabel, ylabel, figure, ylim

def MatJac(th0, th1, th2, dt0, dt1, dt2, d1, a1, a2):
    ths = Matrix([dt0, dt1, dt2])    
    jacob = Matrix([[sin(th0)*(-a2*cos(th1 + th2) + a1*cos(th1)), 
                     -cos(th0)*(a2*sin(th1 + th2) + a1*sin(th1)), 
                     -a2*cos(th0)*sin(th1 + th2)],
                    [cos(th0)*(a2*cos(th1 + th2) + a1*cos(th1)), 
                     -sin(th0)*(a2*sin(th1 + th2) + a1*sin(th1)), 
                     -a2*sin(th0)*sin(th1 + th2)],
                    [0, -a1*cos(th1) - a2*cos(th1 + th2), -a2*cos(th1 + th2)],
                    [0, -sin(th0), -sin(th0)],
                    [0, cos(th0), cos(th0)],
                    [1, 0, 0]])
    return jacob*ths

def validate(angles, veloc, unit, d1, a1, a2):
    res = Matrix()
    error = Matrix()
    mats = Matrix()
    for i in range(0, len(angles.col(0))):
        th0, th1, th2, dt0, dt1, dt2 = angles.row(i)
        if unit == 'd':
            th0, th1, th2, dt0, dt1, dt2 = radians(th0), radians(th1), 
            radians(th2), radians(dt0), radians(dt1), radians(dt2)
        mat = MatJac(th0, th1, th2, dt0, dt1, dt2, d1, a1, a2)
        mats  = Matrix([mats, [mat.T]])
        resy = mat[:3, :]
        err1, err2, err3 = veloc.row(i) - resy.T
        error = Matrix([error, [abs(err1), abs(err2), abs(err3)]])
        res = Matrix([res, resy.T])

    xAxis = range(1, len(res.col(0))+1)
    avg = []
    for i in xAxis:
        ex, ey, ez = error.row(i-1)
        avg.append((ex + ey + ez)/3)

    vx0, = plot(xAxis, veloc.col(0), label='vx0')
    vy0, = plot(xAxis, veloc.col(1), label='vy0')
    vz0, = plot(xAxis, veloc.col(2), label='vz0')
    vx1, = plot(xAxis, res.col(0), label='vx1')
    vy1, = plot(xAxis, res.col(1), label='vy1')
    vz1, = plot(xAxis, res.col(2), label='vz1')
    legend(handles=[vx0, vy0, vz0, vx1, vy1, vz1])
    xlabel("Sample")
    ylabel("Velocity [mm/s]")
    figure()

    ex, = plot(xAxis, error.col(0), label='ex')
    ey, = plot(xAxis, error.col(1), label='ey')
    ez, = plot(xAxis, error.col(2), label='ez')
    avg, = plot(xAxis, avg, label='average')
    legend(handles=[ex, ey, ez, avg])
    xlabel("Sample")
    ylabel("Error [mm/s]")
    #ylim(0, 50)
    show()
    return mats