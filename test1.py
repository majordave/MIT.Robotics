from kine import Matrix, validate, valinv
import csv
# Testing forward kinematics validation method

d1 = 37.25
a1 = 110.32
a2 = 35.65
print(d1 + a1 + a2)
theta = csv.reader(open('theta.csv'), delimiter=',')
pos = csv.reader(open('pos.csv'), delimiter=',')
theta = Matrix(list(theta))
pos = Matrix(list(pos))
mats = validate(theta, pos, 'd', d1, a1, a2)
valinv(mats, theta, d1, a2)