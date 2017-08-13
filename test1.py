import csv
from lib import kine

# Testing forward kinematics validation method
d1, a1, a2 = 46.25, 82.25, 34
print(d1 + a1 + a2)
theta = csv.reader(open('data/theta1.csv'), delimiter=',')
pos = csv.reader(open('data/pos.csv'), delimiter=',')
theta = kine.array(list(theta), dtype=float)
pos = kine.array(list(pos), dtype=float)
mats = kine.validate(theta, pos, d1, a1, a2)
kine.valinv(mats, theta, d1)
