#include "thermocouple.h"

Thermocouple::Thermocouple(int I2C_ADDRESS) {
  this->I2C_ADDRESS = I2C_ADDRESS;
  checkConnection();
  mcp.setThermocoupleType(MCP9600_TYPE_K);
}

static float Thermocouple::cToF(float c) {
  return c * 9.0 / 5.0 + 32.0;
}

bool Thermocouple::checkConnection() {
  if (!mcp.begin(I2C_ADDRESS)) {
    return false;
  }
  return true;
}

float Thermocouple::readHot() {
  return mcp.readThermocouple();
}

float Thermocouple::readCold() {
  return mcp.readAmbient();
}

float Thermocouple::readADC() {
  return mcp.readADC() * 2;
}