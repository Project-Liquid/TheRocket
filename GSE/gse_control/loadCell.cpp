#include "loadCell.h"

LoadCell::LoadCell(int DT_PIN, int SCK_PIN, int gain) {
  this->gain = gain;
  scale.begin(DT_PIN, SCK_PIN);
  scale.tare();
} 

float LoadCell::getAverage(int samples = 100) {
  float sum = 0;
  int count = 0;

  for (int i = 0; i < samples; i++) {
    if (scale.is_ready()) {
      sum += scale.get_units();
      count++;
    }
    delay(20);  // longer spacing for stability
  }

  if (count == 0) return 0;
  return sum / count;
}

float LoadCell::calibrateCell(float knownWeight) {
  Serial.println("Stabilizing...");
  delay(5000);

  float m1 = getAverage();
  Serial.print("Measurement 1: ");
  Serial.println(m1);

  delay(3000);

  float m2 = getAverage();
  Serial.print("Measurement 2: ");
  Serial.println(m2);

  delay(3000);

  float m3 = getAverage();
  Serial.print("Measurement 3: ");
  Serial.println(m3);

  float avgReading = (m1 + m2 + m3) / 3.0;

  Serial.print("Average reading: ");
  Serial.println(avgReading);

  float calFactor = avgReading / knownWeight;

  Serial.print("Calibration Factor: ");
  Serial.println(calFactor);
  cal = calFactor;

  return calFactor;
}

float LoadCell::read() {
  float raw = getAverage(50);
  float weight_grams = raw / cal;
  float weight_lbs = weight_grams / 453.592;
  return weight_lbs;
}
