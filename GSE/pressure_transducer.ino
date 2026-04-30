const int SENSOR_PIN = A0;

// ADC resolution (Teensy 4.1 supports 8, 10, 12, or 16 bit)
const int ADC_BITS = 12;
const int ADC_MAX  = (1 << ADC_BITS) - 1;  // 4095 for 12-bit

// Teensy 4.1 ADC reference voltage
const float VREF = 3.3f;

// Voltage divider: R1=15kΩ, R2=27kΩ
const float DIVIDER_RATIO = 27.0f / (15.0f + 27.0f);  // 0.6429

// Sensor electrical range (at the sensor output, before divider)
const float V_SENSOR_MIN = 1.0f;   // volts -> 0 bar
const float V_SENSOR_MAX = 5.0f;   // volts -> 350 bar

// Sensor pressure range (bar)
const float P_MIN_BAR = 0.0f;
const float P_MAX_BAR = 350.0f;

// Bar to PSI conversion factor
const float BAR_TO_PSI = 14.5038f;

// Number of software samples to average per reading
// Higher = smoother but slightly slower response
const int NUM_SAMPLES = 128;

// How often to print a reading (milliseconds)
const unsigned long PRINT_INTERVAL_MS = 100;

// ── Calibration ───────────────────────────────────────────────────────────────
// pressureOffset_bar: set so that atmospheric reads as 1.0 bar
// Formula: offset = -(raw_reading_at_atmosphere - 1.0)
// Adjust this value if your zero drifts:
float pressureOffset_bar = -182.6f;
float pressureScale      = 1.0f;

// ── Globals ───────────────────────────────────────────────────────────────────
unsigned long lastPrintTime = 0;

// ── Function prototypes ───────────────────────────────────────────────────────
float readAveragedVoltage(int pin, int numSamples);
float sensorVoltageToPressure(float vSensor_V);
bool  isSensorFault(float vSensor_V);

// ─────────────────────────────────────────────────────────────────────────────

void setup() {
  Serial.begin(115200);
  while (!Serial && millis() < 3000);

  analogReadResolution(ADC_BITS);

  Serial.println("=== Pressure Transducer Reader ===");
  Serial.print("Resistors      : R1=15k, R2=27k  ratio=");
  Serial.println(DIVIDER_RATIO, 4);
  Serial.print("ADC resolution : "); Serial.print(ADC_BITS); Serial.println(" bit");
  Serial.print("VREF           : "); Serial.print(VREF, 2);  Serial.println(" V");
  Serial.print("Sensor range   : "); Serial.print(V_SENSOR_MIN);
  Serial.print(" - ");               Serial.print(V_SENSOR_MAX); Serial.println(" V");
  Serial.print("Pressure range : "); Serial.print(P_MIN_BAR);
  Serial.print(" - ");               Serial.print(P_MAX_BAR);    Serial.println(" bar");
  Serial.print("               : "); Serial.print(P_MIN_BAR * BAR_TO_PSI, 1);
  Serial.print(" - ");               Serial.print(P_MAX_BAR * BAR_TO_PSI, 1); Serial.println(" psi");
  Serial.println("==================================");
  Serial.println("ADC_raw, V_pin(V), V_sensor(V), Pressure(bar), Pressure(psi), Status");
}

void loop() {
  unsigned long now = millis();

  if (now - lastPrintTime >= PRINT_INTERVAL_MS) {
    lastPrintTime = now;

    // 1. Read averaged voltage at the ADC pin (after divider)
    float vPin_V = readAveragedVoltage(SENSOR_PIN, NUM_SAMPLES);

    // 2. Recover actual sensor output voltage (reverse the divider)
    float vSensor_V = vPin_V / DIVIDER_RATIO;

    // 3. Raw ADC count for diagnostics
    int rawADC = (int)((vPin_V / VREF) * ADC_MAX + 0.5f);

    // 4. Fault detection
    const char* status;
    float pressure_bar = 0.0f;
    float pressure_psi = 0.0f;

    if (isSensorFault(vSensor_V)) {
      status = (vSensor_V < V_SENSOR_MIN) ? "FAULT:OPEN" : "FAULT:SHORT";
      pressure_bar = -1.0f;
      pressure_psi = -1.0f;
    } else {
      // 5. Convert to pressure in bar
      pressure_bar = sensorVoltageToPressure(vSensor_V);

      // 6. Apply calibration offset (zeroed so atmosphere = 1.0 bar)
      pressure_bar = (pressure_bar + pressureOffset_bar) * pressureScale;

      // 7. Clamp to valid range
      pressure_bar = constrain(pressure_bar, 0.0f, P_MAX_BAR);

      // 8. Convert to PSI
      pressure_psi = pressure_bar * BAR_TO_PSI;

      status = "OK";
    }

    // 9. Print CSV line
    Serial.print(rawADC);          Serial.print(", ");
    Serial.print(vPin_V,    4);    Serial.print(", ");
    Serial.print(vSensor_V, 4);    Serial.print(", ");
    if (pressure_bar >= 0.0f) {
      Serial.print(pressure_bar, 2);
      Serial.print(", ");
      Serial.print(pressure_psi, 2);
    } else {
      Serial.print("---, ---");
    }
    Serial.print(", ");
    Serial.println(status);
  }
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/*
 * Read numSamples ADC values and return the averaged voltage at the pin.
 */
float readAveragedVoltage(int pin, int numSamples) {
  long sum = 0;
  for (int i = 0; i < numSamples; i++) {
    sum += analogRead(pin);
  }
  float avgRaw = (float)sum / numSamples;
  return (avgRaw / ADC_MAX) * VREF;
}

/*
 * Linear interpolation: maps sensor voltage (1-5V) to pressure (0-350 bar).
 */
float sensorVoltageToPressure(float vSensor_V) {
  return P_MIN_BAR +
         (P_MAX_BAR - P_MIN_BAR) *
         (vSensor_V - V_SENSOR_MIN) /
         (V_SENSOR_MAX - V_SENSOR_MIN);
}

/*
 * Returns true if the sensor voltage is outside the valid 1-5V window.
 * Below 0.8V = open wire. Above 5.2V = short circuit.
 */
bool isSensorFault(float vSensor_V) {
  const float FAULT_LOW  = 0.80f;
  const float FAULT_HIGH = 5.20f;
  return (vSensor_V < FAULT_LOW || vSensor_V > FAULT_HIGH);
}
