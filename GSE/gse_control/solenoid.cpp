#include "solenoid.h"

Solenoid::Solenoid(int pin) {
  this->pin = pin;
  pinMode(pin, OUTPUT);
  digitalWrite(pin, LOW);
  is_open = false;
}

Solenoid::open() {
  is_open = true;
  digitalWrite(pin, HIGH);
}

Solenoid::close() {
  is_open = false;
  digitalWrite(pin, LOW);
}

Solenoid::toggle() {
  is_open = !is_open;
  digitalWrite(pin, (is_open ? HIGH : LOW));
}