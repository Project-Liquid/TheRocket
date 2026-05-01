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

void LoadCell::setCalFactor(float calFactor) {
  cal = calFactor;
}

float LoadCell::read(int samples) {
  float raw = getAverage(samples);
  float weight_grams = raw / cal;
  float weight_lbs = weight_grams / 453.592;
  return weight_lbs;
}

static void LoadCell::calibrateCells(LoadCell &scale1, LoadCell &scale2, LoadCell &scale3) {
  static int step = 1;
  if (Serial.available() > 0) {

    float knownWeight = Serial.parseFloat();

    while (Serial.available()) Serial.read();

    if (knownWeight <= 0) {
      Serial.println("Invalid weight. Try again.");
      return;
    }

    if (step == 1) {
      float cal1 = scale1.calibrateCell(knownWeight);

      Serial.println("\nMove SAME weight to Load Cell 2.");
      Serial.println("Enter weight again:");
      step = 2;
    }

    else if (step == 2) {
      float cal2 = scale2.calibrateCell(knownWeight);

      Serial.println("\nMove SAME weight to Load Cell 3.");
      Serial.println("Enter weight again:");
      step = 3;
    }

    else if (step == 3) {
      float cal3 = scale3.calibrateCell(knownWeight);

      Serial.println("\n=== CALIBRATION COMPLETE ===");
      Serial.println("Record these 3 calibration factors.");

      //while (1); // stop forever
    }
  }
}

void LoadCell::join(LoadCell* LC2, LoadCell* LC3) {
  this->LC2 = LC2;
  this->LC3 = LC3;
  joint = true;
}

float LoadCell::readJoint(int samples) {
  if (joint) {
    return this->read(samples) + LC2->read(samples) + LC3->read(samples);
  }
  return -1;
}

void LoadCell::setRedline(float min_weight, int max_counts) {
  this->redline_weight = min_weight;
  this->redline_counts_threshold = max_counts;
}

bool LoadCell::checkRedline() {
  if (!joint) return false;
  float weight = readJoint(1);

  if (weight < redline_weight) {
    redline_counts++;
    Serial.print("Extreme weight: "); Serial.println(weight);
  } else if (redline_counts > 0) {
    redline_counts--;
  }

  if (redline_counts > redline_counts_threshold) {
    redline_counts = 0;
    return true;
  }
  return false;
}

