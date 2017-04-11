#include <Servo.h>

int pos = 90;
int servoPin = 9;
int servoDelay = 15;

Servo servo1;

void setup() {
  Serial.begin(9600);
  servo1.attach(servoPin);
}

void loop() {
  while (Serial.available() == 0){};
  pos = Serial.read();
  servo1.write(pos);
  delay(servoDelay);
}
