#include "transducer.h"
#include "solenoid.h"
#include "mbv.h"
#include "loadCell.h"
#include "thermocouple.h"
#include "heater.h"
#include "thrustCell.h"
//#include <string>

// Run Control
bool print_data = false;
long start_time = 1410065408; // max value placeholder
static unsigned long lastPressureMs = 0;
const int runtime = (10*60) + 8;
const int log_interval_ms = 500;
const bool full_output = true;

// Pressure Transducers
const int ETHANE_UPSTREAM_PIN = A12;
const int ETHANE_DOWNSTREAM_PIN = A11;
const int NITROUS_UPSTREAM_PIN = A10;
const int NITROUS_DOWNSTREAM_PIN = A9;
const int REROUTE_PT_PIN = A8;

const float P_MIN = 0.0;
const float P_MAX_ETHANE = 1000;//68.9;
const float P_MAX_NITROUS = 1500;//103.4214;

Transducer EthaneUpstreamPT(ETHANE_UPSTREAM_PIN, P_MIN, P_MAX_ETHANE);
Transducer EthaneDownstreamPT(ETHANE_DOWNSTREAM_PIN, P_MIN, P_MAX_ETHANE);
Transducer NitrousUpstreamPT(NITROUS_UPSTREAM_PIN, P_MIN, P_MAX_NITROUS);
Transducer NitrousDownstreamPT(NITROUS_DOWNSTREAM_PIN, P_MIN, P_MAX_NITROUS);
Transducer ReroutePT(REROUTE_PT_PIN, P_MIN, P_MAX_NITROUS);

// Solenoids
const int ETHANE_RUN_PIN = 47;
const int ETHANE_VENT_PIN = 46;
const int NITROUS_RUN_PIN = 52;
const int NITROUS_VENT_PIN = 41;

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

// ThrustCell ThrustLC;

// Tank Heaters
// bool heaters_active = false;
// Heater Heater1(49, 0x66);
// Heater Heater2(51, 0x65);

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
  Serial.flush();
  delay(2000);
  //analogReference(INTERNAL1V1);

  EthaneLC1.setCalFactor(1.0);
  EthaneLC2.setCalFactor(1.0);
  EthaneLC3.setCalFactor(1.0);
  NitrousLC1.setCalFactor(1.0);
  NitrousLC2.setCalFactor(1.0);
  NitrousLC3.setCalFactor(1.0);

  //ThrustLC.setCalFactor(5.83);
  
  delay(1200);

  if(full_output) {
    Serial.println("Commands:");
    Serial.println("  ETHANE_VENT_ON / ETHANE_VENT_OFF");
    Serial.println("  NITROUS_VENT_ON / NITROUS_VENT_OFF");
    Serial.println("  ETHANE_RUN_ON / ETHANE_RUN_OFF");
    Serial.println("  NITROUS_RUN_ON / NITROUS_RUN_OFF");
  }
}

String cmd = "";

void loop()
{
  //Serial.flush();
  long elapsed = millis()-start_time;

  // LOG
  if (millis() - lastPressureMs >= log_interval_ms && print_data)
  {
    lastPressureMs += log_interval_ms;
    if(full_output) {
      Serial.println("====================");
      // Time
      Serial.print(elapsed/1000.0); Serial.println("s");
      // PTs
      Serial.print("Ethane Upstream: \t");
      EthaneUpstreamPT.status();
      Serial.print("Ethane Downstream: \t");
      EthaneDownstreamPT.status();
      Serial.print("Nitrous Upstream: \t");
      NitrousUpstreamPT.status();
      Serial.print("Nitrous Downstream: \t");
      NitrousDownstreamPT.status();
      // Load Cells
      Serial.print("Ethane: \t");
      Serial.print("LC1: "); Serial.print(EthaneLC1.read(1)); 
      Serial.print("\tLC2: "); Serial.print(EthaneLC2.read(1)); 
      Serial.print("\tLC3: "); Serial.print(EthaneLC3.read(1)); 
      // Serial.print("\tTotal: "); Serial.print(readEthaneLC());
      Serial.println();
      // Serial.print("Nitrous: \t");
      // Serial.print("LC1: "); Serial.print(NitrousLC1.read(1)); 
      // Serial.print("\tLC2: "); Serial.print(NitrousLC2.read(1)); 
      // Serial.print("\tLC3: "); Serial.print(NitrousLC3.read(1)); 
      // Serial.print("\tTotal: "); Serial.print(readNitrousLC());
      Serial.println();
      // Tank Heaters
      // Serial.print("Tank 1: "); 
      // Serial.print(Heater1.getTemp()); Serial.print("F -- "); 
      // Serial.println(Heater1.isOn() ? "ON" : "OFF");
      // Serial.print("Tank 2: "); 
      // Serial.print(Heater2.getTemp()); Serial.print("F -- "); 
      // Serial.println(Heater2.isOn() ? "ON" : "OFF");
      // Thermocouple
      //Serial.print("Reroute TC: ");
      //Serial.print(RerouteTC.readHot()); Serial.print("F");
      Serial.println();
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
      Serial.print(", ");
      Serial.print(EthaneLC1.read(1));
      Serial.print(", ");
      Serial.print(EthaneLC2.read(1));
      Serial.print(", ");
      Serial.print(EthaneLC3.read(1));
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

  // if(heaters_active) {
  //   Heater1.update();
  //   Heater2.update();
  // }

  delay(20);
}