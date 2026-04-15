#include "transducer.h"
#include "solenoid.h"
#include "mbv.h"
#include "loadCell.h"
#include "thermocouple.h"
#include "heater.h"
#include "thrustCell.h"

// Run Control
bool print_data = false;
long start_time = 1410065408; // max value placeholder
static unsigned long lastPressureMs = 0;
const int runtime = (10*60) + 8;
const int log_interval_ms = 500;
const bool full_output = false;

// Pressure Transducers
const int ETHANE_UPSTREAM_PIN = A12;
const int ETHANE_DOWNSTREAM_PIN = A11;
const int NITROUS_UPSTREAM_PIN = A10;
const int NITROUS_DOWNSTREAM_PIN = A8;
const int REROUTE_PT_PIN = A8;
const float V_REF = 1.1;
const float R_SHUNT = 150.0;

const float P_MIN = 0.0;
const float P_MAX_ETHANE = 68.9;
const float P_MAX_NITROUS = 103.4214;

Transducer EthaneUpstreamPT(ETHANE_UPSTREAM_PIN, P_MIN, P_MAX_ETHANE);
Transducer EthaneDownstreamPT(ETHANE_DOWNSTREAM_PIN, P_MIN, P_MAX_ETHANE);
Transducer NitrousUpstreamPT(NITROUS_UPSTREAM_PIN, P_MIN, P_MAX_NITROUS);
Transducer NitrousDownstreamPT(NITROUS_DOWNSTREAM_PIN, P_MIN, P_MAX_NITROUS);
Transducer ReroutePT(REROUTE_PT_PIN, P_MIN, P_MAX_NITROUS);

// Solenoids
const int ETHANE_RUN_PIN = 47;
const int ETHANE_VENT_PIN = 46;
const int NITROUS_RUN_PIN = 53;
const int NITROUS_VENT_PIN = 52;

Solenoid EthaneRunValve(ETHANE_RUN_PIN);
Solenoid EthaneVent(ETHANE_VENT_PIN);
Solenoid NitrousRunValve(NITROUS_RUN_PIN);
Solenoid NitrousVent(NITROUS_VENT_PIN);

// Motorized Ball Valves
const int ETHANE_MBV_PIN = 6;
const int NITROUS_MBV_PIN = 7;

MBV EthaneMBV(ETHANE_MBV_PIN, 44, 42);
MBV NitrousMBV(NITROUS_MBV_PIN, 50, 48);

// Load Cells
LoadCell EthaneLC1(22, 21);
LoadCell EthaneLC2(24, 23);
LoadCell EthaneLC3(26, 25);
LoadCell NitrousLC1(34, 35);
LoadCell NitrousLC2(32, 33);
LoadCell NitrousLC3(30, 31);

ThrustCell ThrustLC;

// Tank Heaters
bool heaters_active = false;
Heater Heater1(49, 0x66);
Heater Heater2(51, 0x65);

// Thermocouple
Thermocouple RerouteTC(0x67);

//===========================FUNCTIONS============================//
void calibrateCells(LoadCell &scale1, LoadCell &scale2, LoadCell &scale3) {
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

float readEthaneLC() {
  return EthaneLC1.read() + EthaneLC2.read() + EthaneLC3.read();
}

float readNitrousLC() {
  return NitrousLC1.read() + NitrousLC2.read() + NitrousLC3.read();
}

//===========================EXECUTION============================//
void setup()
{
  Serial.begin(9600);
  analogReference(INTERNAL1V1);

  EthaneLC1.setCalFactor(1.0);
  EthaneLC2.setCalFactor(1.0);
  EthaneLC3.setCalFactor(1.0);
  NitrousLC1.setCalFactor(1.0);
  NitrousLC2.setCalFactor(1.0);
  NitrousLC3.setCalFactor(1.0);

  ThrustLC.setCalFactor(5.83);

  delay(1200);

  Serial.println("Commands:");
  Serial.println("  ETHANE_VENT_ON / ETHANE_VENT_OFF");
  Serial.println("  NITROUS_VENT_ON / NITROUS_VENT_OFF");
  Serial.println("  ETHANE_RUN_ON / ETHANE_RUN_OFF");
  Serial.println("  NITROUS_RUN_ON / NITROUS_RUN_OFF");
}

String cmd = "";

void loop()
{
  long elapsed = millis()-start_time;

  if (millis() - lastPressureMs >= log_interval_ms && print_data)
  {
    lastPressureMs += log_interval_ms;
    if(full_output) {
      Serial.println("====================");
      Serial.print(elapsed/1000.0); Serial.println("s");
      Serial.print("Ethane Upstream: \t");
      EthaneUpstreamPT.status();
      Serial.print("Ethane Downstream: \t");
      EthaneDownstreamPT.status();
      Serial.print("Nitrous Upstream: \t");
      NitrousUpstreamPT.status();
      Serial.print("Nitrous Downstream: \t");
      NitrousDownstreamPT.status();
    } else {
      Serial.print(elapsed/1000.0);
      Serial.print(", ");
      EthaneUpstreamPT.value();
      Serial.print(", ");
      EthaneDownstreamPT.value();
      Serial.print(", ");
      NitrousUpstreamPT.value();
      Serial.print(", ");
      NitrousDownstreamPT.value();
      Serial.println();
    }
  }

  if (elapsed/1000.0 > runtime && runtime >= 0) {
    print_data = false;
  }

  // COMMANDS
  if (Serial.available())
  {
    cmd = Serial.readStringUntil('\n');
    cmd.trim();

    // vent
    if (cmd.equalsIgnoreCase("ETHANE_VENT_ON")) {
      EthaneVent.open();
      if(full_output) Serial.println("ETHANE VENT OPEN");
    }
    else if (cmd.equalsIgnoreCase("ETHANE_VENT_OFF")) {
      EthaneVent.close();
      if(full_output) Serial.println("ETHANE VENT CLOSED");
    } 
    else if (cmd.equalsIgnoreCase("NITROUS_VENT_ON")) {
      NitrousVent.open();
      if(full_output) Serial.println("NITROUS VENT OPEN");
    }
    else if (cmd.equalsIgnoreCase("NITROUS_VENT_OFF")) {
      NitrousVent.close();
      if(full_output) Serial.println("NITROUS VENT CLOSED");
    }

    // solenoid
    else if (cmd.equalsIgnoreCase("ETHANE_RUN_ON")) {
      EthaneRunValve.open();
      if(full_output) Serial.println("ETHANE RUN VALVE OPEN");
    }
    else if (cmd.equalsIgnoreCase("ETHANE_RUN_OFF")) {
      EthaneRunValve.close();
      if(full_output) Serial.println("ETHANE RUN VALVE CLOSED");
    } 
    else if (cmd.equalsIgnoreCase("NITROUS_RUN_ON")) {
      NitrousRunValve.open();
      if(full_output) Serial.println("NITROUS RUN VALVE OPEN");
    } 
    else if (cmd.equalsIgnoreCase("NITROUS_RUN_OFF")) {
      NitrousRunValve.close();
      if(full_output) Serial.println("NITROUS RUN VALVE CLOSED");
    } 
    
    // run control
    else if (cmd.equalsIgnoreCase("START")) {
      print_data = true;
      start_time = millis();
      lastPressureMs = start_time;
    } 
    else if (cmd.equalsIgnoreCase("STOP")) {
      print_data = false;
    }
    
    // ball
    else if (cmd.equalsIgnoreCase("ETHANE_A")) {


      EthaneMBV.next_90();
    } 
    else if (cmd.equalsIgnoreCase("ETHANE_S")) {
      EthaneMBV.move_10();
    } 
    else if (cmd.equalsIgnoreCase("ETHANE_D")) {
      EthaneMBV.move_360();
    } 
    else if (cmd.equalsIgnoreCase("ETHANE_R")) {
      EthaneMBV.reset();
      if(full_output) Serial.println("Position reset to zero");
    }

    else {
      if(full_output) Serial.println("Unknown command.");
    }
  }

  EthaneMBV.update();
  NitrousMBV.update();

  if(heaters_active) {
    Heater1.update();
    Heater2.update();
  }

  delay(10);
}