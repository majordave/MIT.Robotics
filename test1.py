import csv
import sympy as sym
from lib import kine

# Testing forward kinematics validation method
d1, a1, a2 = 46.30, 82.30, 33.5
print(d1 + a1 + a2)
theta = csv.reader(open('data/theta1.csv'), delimiter=',')
pos = csv.reader(open('data/pos.csv'), delimiter=',')
theta = sym.Matrix(list(theta))
pos = sym.Matrix(list(pos))
mats = kine.validate(theta, pos, 'd', d1, a1, a2)
kine.valinv(mats, theta, 'd', d1, a2)