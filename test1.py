import kine
import sympy as sym
import csv

# Testing forward kinematics validation method
d1 = 46.25
a1 = 82.25
a2 = 34
print(d1 + a1 + a2)
theta = csv.reader(open('theta3.csv'), delimiter=',')
pos = csv.reader(open('pos1.csv'), delimiter=',')
theta = sym.Matrix(list(theta))
pos = sym.Matrix(list(pos))
mats = kine.validate(theta, pos, 'd', d1, a1, a2)
kine.valinv(mats, theta, d1, a2)