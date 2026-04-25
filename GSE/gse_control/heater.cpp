#include "heater.h"

Heater::Heater(int relay_pin, Transducer* PT) 
  : relay_pin(relay_pin), PT(PT) {
  pinMode(relay_pin, OUTPUT);
  off();
}

void Heater::setTarget(float target_pressure, float threshold_error) {
  this->target_pressure = target_pressure; 
  this->threshold_error = threshold_error;
}

bool Heater::isOn() {
  return blanketOn;
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
  float pressure = PT->readPressure();
  // Hysteresis control
  if (!blanketOn && pressure < (target_pressure - threshold_error)) {
    on();
  } else if (blanketOn && pressure > (target_pressure - threshold_error)) {
    off();
  }
  // Serial.print(pressure, 4);
  // Serial.print(",");
  // Serial.println(blanketOn ? "ON" : "OFF");
}