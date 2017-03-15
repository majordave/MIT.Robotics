from velo import validate, Matrix
import csv

d1 = 37.25
a1 = 110.32
a2 = 35.65
print(d1 + a1 + a2)
theta = csv.reader(open('theta.csv'), delimiter=',')
veloc = csv.reader(open('veloc.csv'), delimiter=',')
theta = Matrix(list(theta))
veloc = Matrix(list(veloc))
mats = validate(theta, veloc, 'd', d1, a1, a2)