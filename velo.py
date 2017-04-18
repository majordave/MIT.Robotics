import sympy as sym
import matplotlib.pyplot as plt
import math

def MatJac(delta0, delta1, unit, d1, a1, a2):
    if unit == 'd':
        delta0[1], delta0[2] = math.radians(delta0[1]), math.radians(delta0[2])
        delta0[3], delta1[1] = math.radians(delta0[3]), math.radians(delta1[1])
        delta1[2], delta1[3] = math.radians(delta1[2]), math.radians(delta1[3])
    
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
    for i in range(1, len(angles[:,0])):
        delta1 = angles[i,:]
        delta0 = angles[i-1,:]
        resy = MatJac(delta0, delta1, unit, d1, a1, a2)
        err1, err2, err3, err4, err5, err6 = veloc[i-1,:] - resy.T
        error = sym.Matrix([error, [abs(err1), abs(err2), abs(err3), abs(err4), abs(err5), abs(err6)]])
        res = sym.Matrix([res, resy.T])

    xAxis = range(1, len(res[:,0])+1)
    avgV = []
    avgW = []
    for i in xAxis:
        evx, evy, evz, ewx, ewy, ewz = error[i-1,:]
        avgV.append((evx + evy + evz)/3)
        avgW.append((ewx + ewy + ewz)/3)

    vx0, = plt.plot(xAxis, veloc[:,0], label=r'real $\nu_x$')
    vy0, = plt.plot(xAxis, veloc[:,1], label=r'real $\nu_y$')
    vz0, = plt.plot(xAxis, veloc[:,2], label=r'real $\nu_z$')
    vx1, = plt.plot(xAxis, res[:,0], label=r'calculated $\nu_x$')
    vy1, = plt.plot(xAxis, res[:,1], label=r'calculated $\nu_y$')
    vz1, = plt.plot(xAxis, res[:,2], label=r'calculated $\nu_z$')
    plt.legend(handles=[vx0, vy0, vz0, vx1, vy1, vz1])
    plt.xlabel("Sample")
    plt.ylabel("Linear Velocity [mm/s]")
    plt.figure()

    evx, = plt.plot(xAxis, error[:,0], label=r'$\nu_x$') 
    evy, = plt.plot(xAxis, error[:,1], label=r'$\nu_y$')
    [evz], [avgV] = plt.plot(xAxis, error[:,2], label=r'$\nu_z$'), plt.plot(xAxis, avgV, label='average')
    plt.legend(handles=[evx, evy, evz, avgV])
    plt.xlabel("Sample")
    plt.ylabel("Error [mm/s]")
    plt.ylim(0, 50)
    plt.figure()

    wx0, = plt.plot(xAxis, veloc[:,3], label=r'real $\omega_x$')
    wy0, = plt.plot(xAxis, veloc[:,4], label=r'real $\omega_y$')
    wz0, = plt.plot(xAxis, veloc[:,5], label=r'real $\omega_z$')
    wx1, = plt.plot(xAxis, res[:,3], label=r'calculated $\omega_x$')
    wy1, = plt.plot(xAxis, res[:,4], label=r'calculated $\omega_y$')
    wz1, = plt.plot(xAxis, res[:,5], label=r'calculated $\omega_z$')
    plt.legend(handles=[wx0, wy0, wz0, wx1, wy1, wz1])
    plt.xlabel("Sample")
    plt.ylabel("Angular Velocity [rad/s]")
    plt.figure()
    
    ewx, = plt.plot(xAxis, error[:,3], label=r'$\omega_x$')
    ewy, = plt.plot(xAxis, error[:,4], label=r'$\omega_y$')
    ewz, = plt.plot(xAxis, error[:,5], label=r'$\omega_z$')
    avgW, = plt.plot(xAxis, avgW, label='average')
    plt.legend(handles=[ewx, ewy, ewz, avgW])
    plt.xlabel("Sample")
    plt.ylabel("Error [rad/s]")
    plt.ylim(0, 50)
    plt.show()

    return res