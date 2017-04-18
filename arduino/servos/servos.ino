#include <Servo.h>
 
Servo servo1;
Servo servo2;
Servo servo3;
int input[2];

void setup()
{
   Serial.begin(9600);
   servo1.attach(9);
   servo2.attach(10);
   servo3.attach(11);
}
 
void loop(){
  if (Serial.available() >= 3){
    for (int i=0; i < 3; i++){
      int inp = Serial.parseInt();
      if  (inp <= 0){
        inp = 1;
      }else if(inp >= 180){
        inp = 179;
      }
      input[i] = inp;
    }   
    servo1.write(input[0]);
    servo2.write(input[1]);
    servo3.write(input[2]);
  }
}
