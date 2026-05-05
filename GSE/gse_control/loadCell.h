#pragma once

#include <Arduino.h>
#include <HX711.h>

class LoadCell {
private:
  HX711 scale;
  float cal = 0;
  int gain;
  LoadCell* LC2 = nullptr;
  LoadCell* LC3 = nullptr;
  bool joint = false;
  float raw = 0;

public:
  LoadCell(int DT_PIN, int SCK_PIN, int gain = 128);
  float getAverage(int samples = 100);
  float calibrateCell(float knownWeight);
  void setCalFactor(float calFactor);
  void tare();
  float read(int samples = 1);
  static void calibrateCells(LoadCell &scale1, LoadCell &scale2, LoadCell &scale3);
  void join(LoadCell* LC2, LoadCell* LC3);
  float readJoint(int samples = 1);
};
