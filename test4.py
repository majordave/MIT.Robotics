from lib import duino, kine, routes
import math

d1, a1, a2 = 46.25, 82.25, 34
ardu = duino.start()
route1 = routes.routine1
cont = 0

for pos in route1:
    if(cont == 4):
        th0 = kine.Th0(pos)
        th1 = kine.Th1(pos, th0, d1, a2)
        th2 = kine.Th2(pos, th0, th1, d1, a2)
        deg = lambda x: str(int(math.degrees(x)))
        duino.move(ardu, deg(th0), deg(th1), deg(th2))
        cont = 0
    else:
        cont += 1
duino.stop(ardu)