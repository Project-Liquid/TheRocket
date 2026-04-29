#include "solenoid.h"

Solenoid::Solenoid(int pin) {
  this->pin = pin;
  pinMode(pin, OUTPUT);
  digitalWrite(pin, LOW);
  is_open = false;
}

void Solenoid::open() {
  is_open = true;
  digitalWrite(pin, HIGH);
}

void Solenoid::close() {
  is_open = false;
  digitalWrite(pin, LOW);
}

void Solenoid::toggle() {
  is_open = !is_open;
  digitalWrite(pin, (is_open ? HIGH : LOW));
}

bool Solenoid::state() {
  return is_open;
}

void Solenoid::setNextActuation(int delay, bool open) {
  scheduled_actuations.push_back({ millis() + (unsigned long)delay, open });
}

void Solenoid::checkScheduledActuation() {
  if (!scheduled_actuations.empty() && millis() >= scheduled_actuations[0].trigger_ms) {
    if (scheduled_actuations[0].open) {
      open();
    } else {
      close();
    }
    scheduled_actuations.erase(scheduled_actuations.begin());
  }
}

void Solenoid::clearSchedule() {
  scheduled_actuations.clear();
}