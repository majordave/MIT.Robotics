import csv
from lib import velo, vinv

# Testing forward velocity validation method
d1, a1, a2 = 46.25, 82.25, 34
print(d1 + a1 + a2)
theta = csv.reader(open('data/theta2.csv'), delimiter=',')
veloc = csv.reader(open('data/veloc.csv'), delimiter=',')
theta = velo.array(list(theta), dtype=float)
veloc = velo.array(list(veloc), dtype=float)
mats = velo.validate(theta, veloc, a1, a2)
th0, th1, th2 = velo.deg2rad(theta[0, 1:])
IJ = vinv.JacInv(th0, th1, th2, d1, a1, a2)
vinv.pprint(IJ)
