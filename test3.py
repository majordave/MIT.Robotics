from lib import duino

while True:
    th0, th1, th2 = input("Enter angles: ").split(" ")
    duino.move(th0, th1, th2)