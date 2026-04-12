#include <Wire.h>
#include <Adafruit_I2CDevice.h>
#include <Adafruit_I2CRegister.h>
#include "Adafruit_MCP9600.h"

class Thermocouple {
private:
  #define I2C_ADDRESS (0x67)
  Adafruit_MCP9600 mcp;
  Ambient_Resolution ambientRes = RES_ZERO_POINT_0625;

public:
  Thermocouple();
  bool checkConnection();
  void setAmbientResolution(Ambient_Resolution res);
  void setThermocoupleType();
  float readHot();
  float readCold();
  float readADC();
};