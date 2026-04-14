#include <Wire.h>
#include <SparkFun_Qwiic_Scale_NAU7802_Arduino_Library.h>

class ThrustCell {
private:
  NAU7802 scale;
  float calibrationFactor = 1.0;
  bool calibrated = false;

public:
  ThrustCell();
  long getAverageReading(int samples = 10);
  void calibrate();
  void setCalFactor(float calFactor);
};