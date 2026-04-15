#include <Wire.h>
#include <Adafruit_MCP9601.h>

//#define I2C_ADDRESS 0x67

class Heater {
private:
  int relay_pin;
  int I2C_ADDRESS;
  bool blanketOn = false;
  Adafruit_MCP9601 mcp;
  const float ON_THRESHOLD_F = 120.5;
  const float OFF_THRESHOLD_F = 121.5;

public:
  Heater(int relay_pin, int I2C_ADDRESS = 0x67);
  float cToF(float c);
  void on();
  void off();
  void update();
};