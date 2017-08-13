from lib import duino, kine

ardu = duino.start()

t = 'a'

if t == 'a':
    tray = [[0, 0, 180],
            [0, 180, 0],
            [180, 0, 180],
            [180, 180, 0]]

else:
    tray = [[0, 90, 0],
            [0, 90, 180]]

trunc = lambda x: x.round(3)

for row in tray:
    th0, th1, th2 = row
    duino.move(ardu, str(th0), str(th1), str(th2))
    pos = kine.PosMat(th0, th1, th2, d1=46.25, a1=82.25, a2=34)
    pos = pos.applyfunc(trunc)
    kine.pprint(pos)

duino.stop(ardu)