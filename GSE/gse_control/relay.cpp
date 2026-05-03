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
  // Check if we have room in the buffer
  if (count < MAX_SCHEDULE) {
    schedule[tail] = { millis() + (unsigned long)delay, open };
    tail = (tail + 1) % MAX_SCHEDULE; // Wrap around
    count++;
  } else {
    Serial.println("Warning: Relay actuation queue is full!");
  }
}

void Relay::checkScheduledActuation() {
  // If there are tasks and the oldest one's time has arrived
  if (count > 0 && millis() >= schedule[head].trigger_ms) {
    
    if (schedule[head].open) {
      open();
    } else {
      close();
    }
    
    head = (head + 1) % MAX_SCHEDULE; // Move head to "erase"
    count--;
  }
}

void Relay::clearSchedule() {
  // Instantly reset the queue pointers
  head = 0;
  tail = 0;
  count = 0;
}

void Relay::neutralize() {
  clearSchedule();
  close();
}
