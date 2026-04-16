#include "heater.h"

Heater::Heater(int relay_pin, int I2C_ADDRESS) {
  this->relay_pin = relay_pin;
  this->I2C_ADDRESS = I2C_ADDRESS;
  pinMode(relay_pin, OUTPUT);
  off();

  Wire.begin();

  if (!mcp.begin(I2C_ADDRESS)) {
    Serial.println("ERROR: MCP9601 not found.");
    while (1) delay(100);
  }
  mcp.setThermocoupleType(MCP9600_TYPE_K);
}

float Heater::cToF(float c) {
  return c * 9.0 / 5.0 + 32.0;
}

bool Heater::isOn() {
  return blanketOn;
}

float Heater::getTemp() {
  float hotC = mcp.readThermocouple();
  //float coldC = mcp.readAmbient();
  float hotF = cToF(hotC);
  return hotF;
}

void Heater::on() {
  digitalWrite(relay_pin, HIGH);
  blanketOn = true;
}

void Heater::off() {
  digitalWrite(relay_pin, LOW);
  blanketOn = false;
}

void Heater::update() {
  float hotF = getTemp();
  // Hysteresis control
  if (!blanketOn && hotF < ON_THRESHOLD_F) {
    on();
  } else if (blanketOn && hotF > OFF_THRESHOLD_F) {
    off();
  }
  // Serial.print(hotC, 4);
  // Serial.print(",");
  // Serial.print(hotF, 4);
  // Serial.print(",");
  // Serial.print(coldC, 4);
  // Serial.print(",");
  // Serial.println(blanketOn ? "ON" : "OFF");
}