#include <Arduino.h>

class Transducer {
private:
  int pin;
  float V_REF = 1.1;
  float P_MIN;
  float P_MAX;
  float R_SHUNT = 150;

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