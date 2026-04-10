#include "mbv.h"

MBV::MBV(int relay_pin, int encoder_pin_1, int encoder_pin_2) {
  this->relay_pin = relay_pin;
  enc = &Encoder(encoder_pin_1, encoder_pin_2);
  current_position = 0;
  target_position = 0;
  moving = false;

  pinMode(relay_pin, OUTPUT);
  digitalWrite(relay_pin, LOW);
  (*enc).write(0);
}

bool MBV::move_90() {
  if (moving) return false;
  current_position = (*enc).read();
  target_position = ((current_position / counts_per_90) + 1) * counts_per_90 - drift;
  moving = true;
  return true;
}

bool MBV::move_small() {
  if (moving) return false;
  current_position = (*enc).read();
  target_position = ((current_position / counts_per_5) + 1) * counts_per_5;
  moving = true;
  return true;
}

bool MBV::move_360() {
  if(moving) return false;
  current_position = (*enc).read();
  target_position = current_position + (long)(360 * counts_per_degree) - drift;
  moving = true;
  return true;
}

void MBV::reset() {
  (*enc).write(0);
  target_position = 0;
}

void MBV::update() {
  if (moving) {
    digitalWrite(relay_pin, HIGH);
    current_position = (*enc).read();
    if (abs(current_position) >= target_position) {
      digitalWrite(relay_pin, LOW);
      moving = false;
      delay(250);
      /*long final_pos = (*enc).read();
      if (full_output) {
        Serial.print("Stopped at: ");
        Serial.print(final_pos);
        Serial.print(" (degrees: ");
        Serial.print(final_pos / counts_per_degree);
        Serial.println(")");
      }*/
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


