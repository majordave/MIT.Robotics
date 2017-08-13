import matplotlib.pyplot as plt
from numpy import *


def MatJac(delta0, delta1, a1, a2, unit='d'):
    if unit == 'd':
        delta0[:, 1:] = deg2rad(delta0[:, 1:])
        delta1[:, 1:] = deg2rad(delta1[:, 1:])

    delta0[:, 2] -= pi / 2
    delta1[:, 2] -= pi / 2

    deltaT = delta1 - delta0
    dtime, dth0, dth1, dth2 = deltaT[:, 0], deltaT[:, 1], deltaT[:, 2], deltaT[:, 3]

    ths = transpose([dth0 / dtime, dth1 / dtime, dth2 / dtime])

    levels = len(delta1[:, 1])
    res = zeros((levels, 6))

    th0, th1, th2 = delta1[:, 1], delta1[:, 2], delta1[:, 3]

    jacob = array(
        [[sin(th0) * (-a2 * cos(th1 + th2) - a1 * cos(th1)), -cos(th0) * (a2 * sin(th1 + th2) + a1 * sin(th1)),
          -a2 * cos(th0) * sin(th1 + th2)],
         [cos(th0) * (a2 * cos(th1 + th2) + a1 * cos(th1)), -sin(th0) * (a2 * sin(th1 + th2) + a1 * sin(th1)),
          -a2 * sin(th0) * sin(th1 + th2)],
         [zeros(levels), -a1 * cos(th1) - a2 * cos(th1 + th2), -a2 * cos(th1 + th2)],
         [zeros(levels), -sin(th0), -sin(th0)],
         [zeros(levels), cos(th0), cos(th0)],
         [ones(levels), zeros(levels), zeros(levels)]]).transpose(2, 0, 1)

    for i in range(0, levels):
        res[i, :] = dot(jacob[i], ths[i, :].T)

    return res


def validate(angles, veloc, a1, a2, unit='d'):
    rows = len(angles[:, 0])
    x_axis = range(1, rows)

    delta1 = angles[x_axis, :]
    delta0 = angles[0:rows - 1, :]
    res = MatJac(delta0, delta1, a1, a2, unit)
    error = abs(veloc - res)
    avgV = (error[:, 0] + error[:, 1] + error[:, 2]) / 3
    avgW = (error[:, 3] + error[:, 4] + error[:, 5]) / 3

    vx0, = plt.plot(x_axis, veloc[:, 0], label=r'real $\nu_x$')
    vy0, = plt.plot(x_axis, veloc[:, 1], label=r'real $\nu_y$')
    vz0, = plt.plot(x_axis, veloc[:, 2], label=r'real $\nu_z$')
    vx1, = plt.plot(x_axis, res[:, 0], label=r'calculated $\nu_x$')
    vy1, = plt.plot(x_axis, res[:, 1], label=r'calculated $\nu_y$')
    vz1, = plt.plot(x_axis, res[:, 2], label=r'calculated $\nu_z$')
    plt.legend(handles=[vx0, vy0, vz0, vx1, vy1, vz1])
    plt.xlabel("Sample")
    plt.ylabel("Linear Velocity [mm/s]")
    plt.figure()

    evx, = plt.plot(x_axis, error[:, 0], label=r'$\nu_x$')
    evy, = plt.plot(x_axis, error[:, 1], label=r'$\nu_y$')
    [evz], [avgV] = plt.plot(x_axis, error[:, 2], label=r'$\nu_z$'), plt.plot(x_axis, avgV, label='average')
    plt.legend(handles=[evx, evy, evz, avgV])
    plt.xlabel("Sample")
    plt.ylabel("Error [mm/s]")
    plt.ylim(0, 50)
    plt.figure()

    wx0, = plt.plot(x_axis, veloc[:, 3], label=r'real $\omega_x$')
    wy0, = plt.plot(x_axis, veloc[:, 4], label=r'real $\omega_y$')
    wz0, = plt.plot(x_axis, veloc[:, 5], label=r'real $\omega_z$')
    wx1, = plt.plot(x_axis, res[:, 3], label=r'calculated $\omega_x$')
    wy1, = plt.plot(x_axis, res[:, 4], label=r'calculated $\omega_y$')
    wz1, = plt.plot(x_axis, res[:, 5], label=r'calculated $\omega_z$')
    plt.legend(handles=[wx0, wy0, wz0, wx1, wy1, wz1])
    plt.xlabel("Sample")
    plt.ylabel("Angular Velocity [rad/s]")
    plt.figure()

    ewx, = plt.plot(x_axis, error[:, 3], label=r'$\omega_x$')
    ewy, = plt.plot(x_axis, error[:, 4], label=r'$\omega_y$')
    ewz, = plt.plot(x_axis, error[:, 5], label=r'$\omega_z$')
    avgW, = plt.plot(x_axis, avgW, label='average')
    plt.legend(handles=[ewx, ewy, ewz, avgW])
    plt.xlabel("Sample")
    plt.ylabel("Error [rad/s]")
    plt.show()

    return res
