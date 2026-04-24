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
  const int MOTOR_SLOW = 215;
  const int SLOW_ZONE = 150;
  const int STOP_ZONE = 5;
  const int TIMEOUT_MS = 3000;
  const int SKIP_ZONE = 50;

  const float counts_per_motor_rev = 751.8;
  const float gear_ratio = (30.0 / 14.0);
  const float counts_per_degree = (counts_per_motor_rev * gear_ratio) / 360.0;
  const float counts_per_90 = (long)(90 * counts_per_degree);

public:
  MBV(int pwm_pin, int encoder_pin_1, int encoder_pin_2);
  bool next_90();
  bool move_degrees(float degrees);
  void reset();
  void update();
  void status();
};

