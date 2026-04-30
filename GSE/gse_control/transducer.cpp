#include "transducer.h"

Transducer::Transducer(int pin, float P_MIN, float P_MAX) {
  this->pin = pin;
  this->V_REF = V_REF;
  this->P_MIN = P_MIN;
  this->P_MAX = P_MAX;
  I_MAX = 20;
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

  // if (i_mA <= 4.0) {
  //   pressure = P_MIN;
  // } else if (i_mA >= I_MAX) {
  //   pressure = P_MAX;
  // } else {
  //   pressure = P_MIN + (i_mA - 4.0) * (P_MAX - P_MIN) / 16.0;
  // }
  pressure = (v - 1.0) * (P_MAX - P_MIN) / (V_MAX - V_MIN) + P_MIN;

  return pressure;
}

String Transducer::status() {
  readPressure();
  String output = "";
  output += "raw = " + String(raw) + " V = " + String(v, 3) + "V  I = " + String(i_mA, 2) + " mA  P = " + String(pressure, 3) + " PSI";
  return output;
}

String Transducer::value() {
  return String(readPressure(), 3);
}

void Transducer::setRedline(float max_pressure, int max_counts) {
  redline_pressure = max_pressure;
  redline_counts_threshold = max_counts;
}

bool Transducer::checkRedline() {
  if (readPressure() > redline_pressure) {
    redline_counts++;
  } else {
    redline_counts--;
  }
  if (redline_counts > redline_counts_threshold) {
    return true;
  }
  return false;
}
