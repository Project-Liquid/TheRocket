#pragma once

#include <Arduino.h>
#include <Encoder.h>

class MBV {
private:
  int pwm_pin;
  long target_position;
  long current_position;
  bool moving;
  long move_start_ms;
  Encoder enc;

  const int MOTOR_FAST = 255;
  const int MOTOR_SLOW = 190;
  const int SLOW_ZONE = 150;
  const int STOP_ZONE = 5;
  const int TIMEOUT_MS = 3000;
  const int SKIP_ZONE = 50;

  const float counts_per_motor_rev = 500;
  const float gear_ratio = (30.0 / 14.0);
  const float ff = 0.965;
  const float counts_per_degree = (ff * counts_per_motor_rev * gear_ratio) / 360.0;
  const float counts_per_90 = (long)(90 * counts_per_degree);

  // --- Circular Buffer Implementation ---
  static const int MAX_SCHEDULE = 10; // Adjust this limit as needed

  struct ScheduledActuation {
    unsigned long trigger_ms;
    float degrees;
  };
  
  ScheduledActuation schedule[MAX_SCHEDULE];
  int head = 0;  // Index of the oldest task (read point)
  int tail = 0;  // Index of the newest task (write point)
  int count = 0; // Current number of tasks in the queue

public:
  MBV(int pwm_pin, int encoder_pin_1, int encoder_pin_2);
  bool next_90();
  bool move_degrees(float degrees);
  void reset();
  void update();
  void status();
  long getCurrentPosition();
  float getCurrentDegrees();
  bool isOpen();
  void setNextActuation(int delay, float degrees = 90);
  void checkScheduledActuation();
  void clearSchedule();
  void neutralize();
};

