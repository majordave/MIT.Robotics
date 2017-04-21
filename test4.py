from lib import duino, kine
import math

d1, a1, a2 = 46.30, 82.30, 33.5
ardu = duino.start()
pos = kine.MatPos(90, 180, 0, d1, a1, a2)
th0 = kine.Th0(pos)
th1 = kine.Th1(pos, th0, d1, a2)
th2 = kine.Th2(pos, th1, d1, a2)
th0 = str(int(math.degrees(th0)))
th1 = str(int(math.degrees(th1)))
th2 = str(int(math.degrees(th2)))
duino.move(ardu, th0, th1, th2)
duino.stop(ardu)
