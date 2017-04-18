import serial
import time
import sys

# Move robot via arduino
def move(th0, th1, th2):

    if sys.platform == 'win32':
        usbport = 'com3'
    elif sys.platform == 'darwin':
        usbport = '/dev/cu.usbmodem1411'

    ardu = serial.Serial(usbport, 9600, timeout=1)
    time.sleep(2)
    angles = " ".join([th0, th1, th2])
    angles = bytearray(bytes(angles, "ascii"))
    ardu.write(angles)