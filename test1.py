from kine import *
from sympy import pprint

# Testing simple forward and inverse kinematics

d1 = 10
a1 = 10
a2 = 10
test = MatPos(0, -pi/2, 0, d1, a1, a2)
pprint(test)
th0 = Th0(test)
th1 = Th1(test, th0, d1, a2)
th2 = Th2(test, th0, th1, d1, a2)
print("th0 = ", th0)
print("th1 = ", th1)
print("th2 = ", th2)
