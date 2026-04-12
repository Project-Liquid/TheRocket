#include "mbv.h"

MBV::MBV(int pwn_pin, int encoder_pin_1, int encoder_pin_2) {
  this->pwm_pin = pwm_pin;
  enc = &Encoder(encoder_pin_1, encoder_pin_2);
  current_position = 0;
  target_position = 0;
  moving = false;

  pinMode(pwm_pin, OUTPUT);
  analogWrite(pwm_pin, 0);
  (*enc).write(0);
}

bool MBV::next_90() {
  current_position = (*enc).read();
  if (moving) return false;
  float current_deg = current_position / counts_per_degree;
  float next_deg = (floor(current_deg / 90.0) + 1.0) * 90.0;
  if ((next_deg - current_deg) * counts_per_degree <= SKIP_ZONE) {
      next_deg += 90.0;
  }
  target_position = (long)(next_deg * counts_per_degree);
  moving = true;
  return true;
}

bool MBV::move_10() {
  if (moving) return false;
  current_position = (*enc).read();
  target_position = current_position + (long)(10 * counts_per_degree);
  moving = true;
  return true;
}

bool MBV::move_360() {
  if(moving) return false;
  current_position = (*enc).read();
  target_position = current_position + (long)(360 * counts_per_degree);
  moving = true;
  return true;
}

void MBV::reset() {
  (*enc).write(0);
  target_position = 0;
}

void MBV::update() {
  if (moving) {
    long error = target_position - (*enc).read();

    if (abs(error) <= STOP_ZONE || millis() - move_start_ms > TIMEOUT_MS) {
      analogWrite(pwm_pin, 0);
      moving = false;
      delay(250);
      Serial.print("Stopped at: ");
      Serial.print((*enc).read() / counts_per_degree);
      Serial.println(" degrees");
    } else if (error < SLOW_ZONE) {
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



// #include "mbv.h"

// MBV::MBV(int relay_pin, int encoder_pin_1, int encoder_pin_2) {
//   this->relay_pin = relay_pin;
//   enc = &Encoder(encoder_pin_1, encoder_pin_2);
//   current_position = 0;
//   target_position = 0;
//   moving = false;

//   pinMode(relay_pin, OUTPUT);
//   digitalWrite(relay_pin, LOW);
//   (*enc).write(0);
// }

// bool MBV::move_90() {
//   if (moving) return false;
//   current_position = (*enc).read();
//   target_position = ((current_position / counts_per_90) + 1) * counts_per_90 - drift;
//   moving = true;
//   return true;
// }

// bool MBV::move_small() {
//   if (moving) return false;
//   current_position = (*enc).read();
//   target_position = ((current_position / counts_per_5) + 1) * counts_per_5;
//   moving = true;
//   return true;
// }

// bool MBV::move_360() {
//   if(moving) return false;
//   current_position = (*enc).read();
//   target_position = current_position + (long)(360 * counts_per_degree) - drift;
//   moving = true;
//   return true;
// }

// void MBV::reset() {
//   (*enc).write(0);
//   target_position = 0;
// }

// void MBV::update() {
//   if (moving) {
//     digitalWrite(relay_pin, HIGH);
//     current_position = (*enc).read();
//     if (abs(current_position) >= target_position) {
//       digitalWrite(relay_pin, LOW);
//       moving = false;
//       delay(250);
//       /*long final_pos = (*enc).read();
//       if (full_output) {
//         Serial.print("Stopped at: ");
//         Serial.print(final_pos);
//         Serial.print(" (degrees: ");
//         Serial.print(final_pos / counts_per_degree);
//         Serial.println(")");
//       }*/
//     }
//   }
// }

// void MBV::status() {
//   Serial.print("Target: ");
//   Serial.print(target_position);
//   Serial.print(" (degrees: ");
//   Serial.print(target_position / counts_per_degree);
//   Serial.println(")");

//   Serial.print("Position: ");
//   Serial.print(current_position);
//   Serial.print(" (degrees: ");
//   Serial.print(current_position / counts_per_degree);
//   Serial.println(")");
// }



