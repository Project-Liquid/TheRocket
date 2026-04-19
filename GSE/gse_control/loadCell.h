#include <Arduino.h>
#include <HX711.h>

class LoadCell {
private:
  HX711 scale;
  float cal = 0;
  int gain;

public:
  LoadCell(int DT_PIN, int SCK_PIN, int gain = 128);
  float getAverage(int samples = 100);
  float calibrateCell(float knownWeight);
  void setCalFactor(float calFactor);
  float read(int samples = 50);
  static void LoadCell::calibrateCells(LoadCell &scale1, LoadCell &scale2, LoadCell &scale3);
};
