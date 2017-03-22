import kine
import sympy as sym
import csv

# Testing forward kinematics validation method
d1 = 37.25
a1 = 110.32
a2 = 35.65
print(d1 + a1 + a2)
theta = csv.reader(open('theta1.csv'), delimiter=',')
pos = csv.reader(open('pos.csv'), delimiter=',')
theta = sym.Matrix(list(theta))
pos = sym.Matrix(list(pos))
mats = kine.validate(theta, pos, 'd', d1, a1, a2)
kine.valinv(mats, theta, d1, a2)