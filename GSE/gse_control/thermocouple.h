#include <Wire.h>
// #include <Adafruit_I2CDevice.h>
// #include <Adafruit_I2CRegister.h>
#include <Adafruit_MCP9601.h>

class Thermocouple {
private:
  int I2C_ADDRESS;
  Adafruit_MCP9601 mcp;
  Ambient_Resolution ambientRes = RES_ZERO_POINT_0625;

public:
  Thermocouple(int I2C_ADDRESS = 0x67);
  static float cToF(float c);
  bool checkConnection();
  float readHot();
  float readCold();
  float readADC();
};