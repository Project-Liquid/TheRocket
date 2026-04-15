#include "thermocouple.h"

Thermocouple::Thermocouple(int I2C_ADDRESS) {
  this->I2C_ADDRESS = I2C_ADDRESS;
  checkConnection();
  //setAmbientResolution(ambientRes);
  setThermocoupleType();
  //mcp.setFilterCoefficient(3);
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

void Thermocouple::setAmbientResolution(Ambient_Resolution res) {
  mcp.setAmbientResolution(res);
  // Serial.print("Ambient Resolution set to: ");
  // switch (ambientRes) {
  //   case RES_ZERO_POINT_25:    Serial.println("0.25°C"); break;
  //   case RES_ZERO_POINT_125:   Serial.println("0.125°C"); break;
  //   case RES_ZERO_POINT_0625:  Serial.println("0.0625°C"); break;
  //   case RES_ZERO_POINT_03125: Serial.println("0.03125°C"); break;
  // }
  // Serial.println(" bits");
}

void Thermocouple::setThermocoupleType() {
  mcp.setThermocoupleType(MCP9600_TYPE_K);
  // Serial.print("Thermocouple type set to ");
  // switch (mcp.getThermocoupleType()) {
  //   case MCP9600_TYPE_K:  Serial.print("K"); break;
  //   case MCP9600_TYPE_J:  Serial.print("J"); break;
  //   case MCP9600_TYPE_T:  Serial.print("T"); break;
  //   case MCP9600_TYPE_N:  Serial.print("N"); break;
  //   case MCP9600_TYPE_S:  Serial.print("S"); break;
  //   case MCP9600_TYPE_E:  Serial.print("E"); break;
  //   case MCP9600_TYPE_B:  Serial.print("B"); break;
  //   case MCP9600_TYPE_R:  Serial.print("R"); break;
  // }
  // Serial.println(" type");
}

float Thermocouple::readHot() {
  return cToF(mcp.readThermocouple());
}

float Thermocouple::readCold() {
  return cToF(mcp.readAmbient());
}

float Thermocouple::readADC() {
  return mcp.readADC() * 2;
}