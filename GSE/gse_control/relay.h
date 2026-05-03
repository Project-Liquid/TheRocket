#pragma once

#include <Arduino.h>
#include <StandardCplusplus.h>
#include <vector>

class Relay {
private:
  int pin;
  int is_open;
  // --- Circular Buffer Implementation ---
  static const int MAX_SCHEDULE = 50; // Adjust this limit as needed

  struct ScheduledActuation {
    unsigned long trigger_ms;
    bool open;

  ScheduledActuation schedule[MAX_SCHEDULE];
  int head = 0;  // Index of the oldest task
  int tail = 0;  // Index of the newest task
  int count = 0; // Current number of tasks in the queue


public:
  Relay(int pin);
  void open();
  void close();
  void toggle();
  bool state();
  void setNextActuation(int delay, bool open);
  void checkScheduledActuation();
  void clearSchedule();
  void neutralize();
};
