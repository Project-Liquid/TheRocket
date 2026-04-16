#include <Arduino.h>

class Transducer {
private:
  int pin;
  float V_REF = 5.0;
  float P_MIN;
  float P_MAX;
  int R_SHUNT = 250;
  int I_MIN = 4;
  int I_MAX = 20;
  float V_MIN = I_MIN * R_SHUNT / 1000;
  float V_MAX = I_MAX * R_SHUNT / 1000;

  int raw;
  float v;
  float i_mA;
  float pressure;

public:
  Transducer(int pin, float P_MIN, float P_MAX);
  static float barToPSI(float bar);
  static float PSIToBar(float psi);
  float readPressure();
  void status();
  void value();
};