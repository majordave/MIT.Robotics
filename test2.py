import velo
import sympy as sym
import csv
import vinv
import math

# Testing forward velocity validation method
d1, a1, a2 = 46.30, 82.30, 33.5
print(d1 + a1 + a2)
theta = csv.reader(open('data/theta2.csv'), delimiter=',')
veloc = csv.reader(open('data/veloc.csv'), delimiter=',')
theta = sym.Matrix(list(theta))
veloc = sym.Matrix(list(veloc))
mats = velo.validate(theta, veloc, 'd', d1, a1, a2)
th0, th1, th2 = theta[0,1:]
th0, th1, th2 = math.radians(th0), math.radians(th1), math.radians(th2)
IJ = vinv.JacInv(th0, th1, th2, d1, a1, a2)
sym.pprint(IJ)
