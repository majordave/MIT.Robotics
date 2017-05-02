import serial
import time
import sys

# Initialize arduino
def start():
    if sys.platform == 'darwin':
        usbport = '/dev/cu.usbmodem1421'
    else:
        usbport = 'com3'
    return serial.Serial(usbport, 9600, timeout=1)

# When finalized moving, set to default position
def stop(robot):
    move(robot, '90', '90', '90')

# Move robot via arduino
def move(robot, th0, th1, th2):
    time.sleep(2)
    angles = " ".join([th0, th1, th2])
    angles = bytearray(bytes(angles, "ascii"))
    robot.write(angles)