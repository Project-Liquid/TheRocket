#include "relay.h"

Relay::Relay(int pin) {
  this->pin = pin;
  pinMode(pin, OUTPUT);
  digitalWrite(pin, LOW);
  is_open = false;
}

void Relay::open() {
  is_open = true;
  digitalWrite(pin, HIGH);
}

void Relay::close() {
  is_open = false;
  digitalWrite(pin, LOW);
}

void Relay::toggle() {
  is_open = !is_open;
  digitalWrite(pin, (is_open ? HIGH : LOW));
}

bool Relay::state() {
  return is_open;
}

void Relay::setNextActuation(int delay, bool open) {
  scheduled_actuations.push_back({ millis() + (unsigned long)delay, open });
}

void Relay::checkScheduledActuation() {
  if (!scheduled_actuations.empty() && millis() >= scheduled_actuations[0].trigger_ms) {
    if (scheduled_actuations[0].open) {
      open();
    } else {
      close();
    }
    scheduled_actuations.erase(scheduled_actuations.begin());
  }
}

void Relay::clearSchedule() {
  scheduled_actuations.clear();
}