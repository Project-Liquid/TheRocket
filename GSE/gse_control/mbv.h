#include <Arduino.h>
#include <Encoder.h>

class MBV {
private:
  int pwm_pin;
  long target_position;
  long current_position;
  bool moving;
  long move_start_ms;
  Encoder* enc = nullptr;

  const int MOTOR_FAST = 255;
  const int MOTOR_SLOW = 215;
  const int SLOW_ZONE = 150;
  const int STOP_ZONE = 5;
  const int TIMEOUT_MS = 3000;
  const int SKIP_ZONE = 50;

  static const long counts_per_motor_rev = 751.8;
  static const long gear_ratio = (30.0 / 14.0);
  static const long counts_per_degree = (counts_per_motor_rev * gear_ratio) / 360.0;
  static const long counts_per_90 = (long)(90 * counts_per_degree);

public:
  MBV(int relay_pin, int encoder_pin_1, int encoder_pin_2);
  bool next_90();
  bool move_10();
  bool move_360();
  void reset();
  void update();
  void status();
};

// class MBV {
// private:
//   int relay_pin;
//   long target_position;
//   long current_position;
//   bool moving;
//   Encoder* enc = nullptr;

//   static const long counts_per_motor_rev = 751.8;
//   static const long gear_ratio = (31.0 / 14.0);
//   static const long counts_per_valve_rev = counts_per_motor_rev * gear_ratio;
//   static const long counts_per_degree = counts_per_valve_rev / 360.0;
//   static const long counts_per_90 = (long)(90 * counts_per_degree);
//   static const long counts_per_5 = (long)(5 * counts_per_degree);
//   static const long drift = 50;

// public:
//   MBV(int relay_pin, int encoder_pin_1, int encoder_pin_2);
//   bool move_90();
//   bool move_small();
//   bool move_360();
//   void reset();
//   void update();
//   void status();
// };
