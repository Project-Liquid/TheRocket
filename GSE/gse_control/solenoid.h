#pragma once

#include <Arduino.h>

class Solenoid {
private:
  int pin;
  int is_open;

public:
  Solenoid(int pin);
  open();
  close();
  toggle();
};