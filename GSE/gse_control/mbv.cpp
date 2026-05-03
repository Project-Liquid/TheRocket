#include "mbv.h"

MBV::MBV(int pwm_pin, int encoder_pin_1, int encoder_pin_2)
  : pwm_pin(pwm_pin), enc(encoder_pin_1, encoder_pin_2) {
  current_position = 0;
  target_position = 0;
  moving = false;

  pinMode(pwm_pin, OUTPUT);
  analogWrite(pwm_pin, 0);
  enc.write(0);
}

bool MBV::next_90() {
  // current_position = enc.read();
  if (moving) return false;
  float current_deg = current_position / counts_per_degree;
  float next_deg = (floor(current_deg / 90.0) + 1.0) * 90.0;
  if ((next_deg - current_deg) * counts_per_degree <= SKIP_ZONE) {
      next_deg += 90.0;
  }
  target_position = (long)(next_deg * counts_per_degree);
  moving = true;
  move_start_ms = millis();
  // Serial.print("Moving to: "); Serial.println(target_position / counts_per_degree);
  return true;
}

bool MBV::move_degrees(float degrees) {
  if(moving) return false;
  // current_position = enc.read();
  target_position = current_position + (long)(degrees * counts_per_degree);
  moving = true;
  move_start_ms = millis();
  return true;
}

void MBV::reset() {
  enc.write(0);
  target_position = 0;
}

/**
 * Checks if at target position
 */
void MBV::update() {
  current_position = enc.read();
  if (moving) {
    long error = target_position - current_position;
    //Serial.print("Error: "); Serial.println(error);

    if (abs(error) <= STOP_ZONE || millis() - move_start_ms > TIMEOUT_MS) {
      analogWrite(pwm_pin, 0);
      moving = false;
      //delay(250);
      // Serial.print("Stopped at: ");
      // Serial.print(current_position / counts_per_degree);
      // Serial.println(" degrees");
    } else if (abs(error) < SLOW_ZONE) {
      analogWrite(pwm_pin, MOTOR_SLOW);
    } else {
      analogWrite(pwm_pin, MOTOR_FAST);
    }
  }
}

void MBV::status() {
  Serial.print("Target: ");
  Serial.print(target_position);
  Serial.print(" (degrees: ");
  Serial.print(target_position / counts_per_degree);
  Serial.println(")");

  Serial.print("Position: ");
  Serial.print(current_position);
  Serial.print(" (degrees: ");
  Serial.print(current_position / counts_per_degree);
  Serial.println(")");
}

long MBV::getCurrentPosition() {
  return current_position;
}

float MBV::getCurrentDegrees() {
  return getCurrentPosition() / counts_per_degree;
} 

bool MBV::isOpen() {
  if (((int)getCurrentDegrees() % 180) < 45 || ((int)getCurrentDegrees() % 180) > 135) {
    return false;
  }
  return true;
}

void MBV::setNextActuation(int delay, float degrees) {
  if (count < MAX_SCHEDULE) {
    schedule[tail] = { millis() + (unsigned long)delay, degrees };
    tail = (tail + 1) % MAX_SCHEDULE; // Wrap around if we hit the end of the array
    count++;
  } else {
    Serial.println("Warning: Actuation queue is full!");
  }
}

void MBV::checkScheduledActuation() {
  if (count > 0 && millis() >= schedule[head].trigger_ms) {
    next_90();
    head = (head + 1) % MAX_SCHEDULE; // Move the head forward to "erase" the task
    count--;
  }
}

void MBV::clearSchedule() {
  head = 0;
  tail = 0;
  count = 0;
}

void MBV::neutralize() {
  clearSchedule();
}
