#include "thrustCell.h"

ThrustCell::ThrustCell() {
  Wire.begin();

  if (!scale.begin()) {
    Serial.println("NAU7802 not found!");
    while (1);
  }
}

long ThrustCell::getAverageReading(int samples) {
  long sum = 0;

  for (int i = 0; i < samples; i++) {
    sum += scale.getReading();
    delay(10); // stable ADC timing
  }

  return sum / samples;
}

void ThrustCell::calibrate() {
  if (Serial.available() && !calibrated) {

    float knownWeight = Serial.parseFloat();

    // 🔥 IMPORTANT: clear leftover newline / garbage
    while (Serial.available()) {
      Serial.read();
    }

    Serial.print("You entered: ");
    Serial.println(knownWeight);

    // safety check (prevents inf)
    if (knownWeight <= 0) {
      Serial.println("ERROR: weight must be > 0");
      return;
    }

    long raw = getAverageReading(20);
    long zero = scale.getZeroOffset();

    calibrationFactor = (raw - zero) / knownWeight;

    Serial.print("Calibration Factor: ");
    Serial.println(calibrationFactor);

    Serial.println("Use this factor in final code.");

    calibrated = true; // 🔒 lock so it never runs again
  }
}

void ThrustCell::setCalFactor(float calFactor) {
  calibrationFactor = calFactor;
}