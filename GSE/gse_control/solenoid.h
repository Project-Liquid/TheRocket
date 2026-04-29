#pragma once

#include <Arduino.h>
#include <StandardCplusplus.h>
#include <vector>

class Solenoid {
private:
  int pin;
  int is_open;
  struct ScheduledActuation {
    unsigned long trigger_ms;
    bool open;
  };
  std::vector<ScheduledActuation> scheduled_actuations;


public:
  Solenoid(int pin);
  void open();
  void close();
  void toggle();
  bool state();
  void setNextActuation(int delay, bool open);
  void checkScheduledActuation();
  void clearSchedule();
};