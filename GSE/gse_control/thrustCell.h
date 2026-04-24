#include <Wire.h>
#include <SparkFun_Qwiic_Scale_NAU7802_Arduino_Library.h>

// I2C Address: 0x2A (required)

class ThrustCell {
private:
  NAU7802 scale;
  float calibrationFactor = 1.0;
  bool calibrated = false;
  int offset = 0;

public:
  ThrustCell();
  long getAverageReading(int samples = 10);
  double read(int samples = 1);
  void calibrate();
  void tare();
  void setCalFactor(float calFactor);
};