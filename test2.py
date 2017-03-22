import velo
import sympy as sym
import csv

# Testing forward velocity validation method
d1 = 37.25
a1 = 110.32
a2 = 35.65
print(d1 + a1 + a2)
theta = csv.reader(open('theta2.csv'), delimiter=',')
veloc = csv.reader(open('veloc.csv'), delimiter=',')
theta = sym.Matrix(list(theta))
veloc = sym.Matrix(list(veloc))
mats = velo.validate(theta, veloc, 'd', d1, a1, a2)