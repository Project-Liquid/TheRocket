// rfd900x_led_control.ino
// Arduino Uno — RFD900x on pins 0/1
// Send "1" to turn LED on  → replies "LED on"
// Send "0" to turn LED off → replies "LED off"
//
// Wiring:
//   RFD900x TX  ->  Arduino pin 0 (RX)
//   RFD900x RX  ->  Arduino pin 1 (TX)
//   RFD900x GND ->  Arduino GND
//   RFD900x 5V  ->  Arduino 5V
//
// NOTE: Disconnect RFD900x from pins 0/1 before uploading, reconnect after.

#define BAUD_RATE 57600
#define LED_PIN   13

String inputBuffer = "";

void setup() {
  Serial.begin(BAUD_RATE);
  pinMode(LED_PIN, OUTPUT);
  digitalWrite(LED_PIN, LOW);
  inputBuffer.reserve(32);
}

void loop() {
  while (Serial.available()) {
    char c = (char)Serial.read();

    if (c == '\n') {
      inputBuffer.trim();

      if (inputBuffer == "1") {
        digitalWrite(LED_PIN, HIGH);
        Serial.println("LED on");
      } else if (inputBuffer == "0") {
        digitalWrite(LED_PIN, LOW);
        Serial.println("LED off");
      }

      inputBuffer = "";
    } else {
      inputBuffer += c;
    }
  }
}
