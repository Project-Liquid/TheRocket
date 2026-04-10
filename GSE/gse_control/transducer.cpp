#include "transducer.h"

  Transducer::Transducer(int pin, float P_MIN, float P_MAX) {
    this->pin = pin;
    this->V_REF = V_REF;
    this->P_MIN = P_MIN;
    this->P_MAX = P_MAX;
  }

  static float Transducer::barToPSI(float bar) {
    return bar * 14.504;
  }

  static float Transducer::PSIToBar(float psi) {
    return psi / 14.504;
  }

  float Transducer::readPressure() {
    raw = analogRead(pin);
    v = raw * V_REF / 1023.0;
    i_mA = (v / R_SHUNT) * 1000.0;
    if (i_mA < 0) {
      i_mA = 0;
    } else if (i_mA > 25) {
      i_mA = 25;
    }

    if (i_mA <= 4.0) {
      pressure = P_MIN;
    } else if (i_mA >= 20.0) {
      pressure = P_MAX;
    } else {
      pressure = P_MIN + (i_mA - 4.0) * (P_MAX - P_MIN) / 16.0;
    }

    pressure = barToPSI(pressure);

    return pressure;
  }

  void Transducer::status() {
    readPressure();
    
    Serial.print("raw = "); Serial.print(raw);
    Serial.print("  V = "); Serial.print(v, 3);
    Serial.print(" V  I = "); Serial.print(i_mA + 0.1, 2);
    Serial.print(" mA  P = "); Serial.print(pressure, 3);
    Serial.println(" PSI");
  }

  void Transducer::value() {
    Serial.print(readPressure(), 3);
  }
