from kine import Matrix, validate
import csv

# Testing forward kinematics validation method

d1 = 37.25
a1 = 138.44
a2 = 35.65
theta = csv.reader(open('theta.csv'), delimiter=',')
pos = csv.reader(open('pos.csv'), delimiter=',')
theta = Matrix(list(theta))
pos = Matrix(list(pos))
validate(theta, pos, 'd', d1, a1, a2)
