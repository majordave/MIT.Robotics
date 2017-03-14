from kine import MatPos, Th0, Th1, Th2
from sympy import pi, pprint

# Testing simple forward and inverse kinematics

d1 = 10
a1 = 10
a2 = 10
test = MatPos(0, -pi/2, 0, d1, a1, a2)
pprint(test)
th0 = Th0(test)
print("th0 = ", th0)
th1 = Th1(test, th0, d1, a2)
print("th1 = ", th1)
th2 = Th2(test, th1, d1, a2)
print("th2 = ", th2)
