#include "transducer.h"
#include "solenoid.h"
#include "mbv.h"
#include "loadCell.h"
#include "thermocouple.h"
#include "heater.h"
#include "thrustCell.h"

// Run Control
bool print_data = true;
long start_time = 1410065408; // max value placeholder
static unsigned long lastPressureMs = 0;
const int runtime = (10*60) + 8;
const int log_interval_ms = 500;
const bool full_output = false;
long elapsed = 0;

// Pressure Transducers
const int ETHANE_UPSTREAM_PIN = A12;
const int ETHANE_DOWNSTREAM_PIN = A11;
const int NITROUS_UPSTREAM_PIN = A10;
const int NITROUS_DOWNSTREAM_PIN = A9;
const int REROUTE_PT_PIN = A8;

const float P_MIN = 0.0;
const float P_MAX_ETHANE = 1000;
const float P_MAX_NITROUS = 1500;

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

MBV* EthaneMBV;
MBV* NitrousMBV;

// Load Cells
LoadCell* EthaneLC2;
LoadCell* EthaneLC1;
LoadCell* EthaneLC3;
LoadCell* NitrousLC1;
LoadCell* NitrousLC2;
LoadCell* NitrousLC3;

ThrustCell* ThrustLC;

// Tank Heaters
bool heaters_active = false;
Heater* EthaneHeater1;
Heater* EthaneHeater2;
Heater* NitrousHeater1;
Heater* NitrousHeater2;

// Thermocouple
Thermocouple* RerouteTC;
Thermocouple* RerouteTC2;
Thermocouple* RerouteTC3;

//===========================FUNCTIONS============================//

float readEthaneLC() {
  return EthaneLC1->read(1) + EthaneLC2->read(1) + EthaneLC3->read(1);
}

float readNitrousLC() {
  return NitrousLC1->read(1) + NitrousLC2->read(1) + NitrousLC3->read(1);
}

void status(bool verbose = true) {
  elapsed = millis()-start_time;

  if (verbose) {  // Full diagnostic output
    Serial.println("====================");
    Serial2.println("====================");
    // Time
    Serial.print(elapsed/1000.0); Serial.println("s");
    Serial2.print(elapsed/1000.0); Serial2.println("s");
    // PTs
    Serial.print("Ethane Upstream: \t");
    Serial2.print("Ethane Upstream: \t");
    EthaneUpstreamPT.status();
    Serial.print("Ethane Downstream: \t");
    Serial2.print("Ethane Downstream: \t");
    EthaneDownstreamPT.status();
    Serial.print("Nitrous Upstream: \t");
    Serial2.print("Nitrous Upstream: \t");
    NitrousUpstreamPT.status();
    Serial.print("Nitrous Downstream: \t");
    Serial2.print("Nitrous Downstream: \t");
    NitrousDownstreamPT.status();
    Serial.print("Reroute: \t");
    Serial2.print("Reroute: \t");
    ReroutePT.status();
    // Load Cells
    Serial.print("Ethane: \t");
    Serial2.print("Ethane: \t");
    Serial.print("LC1: "); Serial.print(EthaneLC1->read(1)); 
    Serial2.print("LC1: "); Serial2.print(EthaneLC1->read(1)); 
    Serial.print("\tLC2: "); Serial.print(EthaneLC2->read(1));
    Serial2.print("\tLC2: "); Serial2.print(EthaneLC2->read(1));
    Serial.print("\tLC3: "); Serial.print(EthaneLC3->read(1)); 
    Serial2.print("\tLC3: "); Serial2.print(EthaneLC3->read(1)); 
    Serial.print("\tTotal: "); Serial.print(readEthaneLC());
    Serial2.print("\tTotal: "); Serial2.print(readEthaneLC());
    Serial.println();
    Serial2.println();
    Serial.print("Nitrous: \t");
    Serial2.print("Nitrous: \t");
    Serial.print("LC1: "); Serial.print(NitrousLC1->read(1)); 
    Serial2.print("LC1: "); Serial2.print(NitrousLC1->read(1)); 
    Serial.print("\tLC2: "); Serial.print(NitrousLC2->read(1)); 
    Serial2.print("\tLC2: "); Serial2.print(NitrousLC2->read(1)); 
    Serial.print("\tLC3: "); Serial.print(NitrousLC3->read(1)); 
    Serial2.print("\tLC3: "); Serial2.print(NitrousLC3->read(1)); 
    Serial.print("\tTotal: "); Serial.print(readNitrousLC());
    Serial2.print("\tTotal: "); Serial2.print(readNitrousLC());
    Serial.println();
    Serial2.println();
    Serial.print("Thrust: "); Serial.print(ThrustLC->read());
    Serial2.print("Thrust: "); Serial2.print(ThrustLC->read());
    Serial.println();
    Serial2.println();
    // Serial.println();
    // Tank Heaters
    Serial.print("Ethane Tank 1: "); 
    Serial2.print("Ethane Tank 1: "); 
    Serial.println(EthaneHeater1->isOn() ? "ON" : "OFF");
    Serial2.println(EthaneHeater1->isOn() ? "ON" : "OFF");
    Serial.print("Nitrous Tank 1: "); 
    Serial2.print("Nitrous Tank 1: "); 
    Serial.println(NitrousHeater1->isOn() ? "ON" : "OFF");
    Serial2.println(NitrousHeater1->isOn() ? "ON" : "OFF");
    // Thermocouple
    Serial.print("Reroute TC: ");
    Serial2.print("Reroute TC: ");
    Serial.print(RerouteTC->readHot()); Serial.print("F");
    Serial2.print(RerouteTC->readHot()); Serial2.print("F");
    Serial.println();
    Serial2.println();
    Serial.print("Reroute TC2: ");
    Serial2.print("Reroute TC2: ");
    Serial.print(RerouteTC2->readHot()); Serial.print("F");
    Serial2.print(RerouteTC2->readHot()); Serial2.print("F");
    Serial.println();
    Serial2.println();
    Serial.print("Reroute TC3: ");
    Serial2.print("Reroute TC3: ");
    Serial.print(RerouteTC3->readHot()); Serial.print("F");
    Serial2.print(RerouteTC3->readHot()); Serial2.print("F");
    Serial.println();
    Serial2.println();
  } else { // Simplified output for log
    // FORMAT: DATA|millis|KEY:VAL|KEY:VAL|...
    Serial.print("DATA|");
    Serial2.print("DATA|");
    Serial.print(elapsed/1000.0);
    Serial2.print(elapsed/1000.0);
    Serial.print("|PT_EU:");
    Serial2.print("|PT_EU:");
    EthaneUpstreamPT.value();
    Serial.print("|PT_ED:");
    Serial2.print("|PT_ED:");
    EthaneDownstreamPT.value();
    Serial.print("|PT_NU:");
    Serial2.print("|PT_NU:");
    NitrousUpstreamPT.value();
    Serial.print("|PT_ND:");
    Serial2.print("|PT_ND:");
    NitrousDownstreamPT.value();
    Serial.print("|LC_E1:");
    Serial2.print("|LC_E1:");
    Serial.print(EthaneLC1->read(1));
    Serial2.print(EthaneLC1->read(1));
    Serial.print("|LC_E2:");
    Serial2.print("|LC_E2:");
    Serial.print(EthaneLC2->read(1));
    Serial2.print(EthaneLC2->read(1));
    Serial.print("|LC_E3:");
    Serial2.print("|LC_E3:");
    Serial.print(EthaneLC3->read(1));
    Serial2.print(EthaneLC3->read(1));
    Serial.print("|LC_N1:");
    Serial2.print("|LC_N1:");
    Serial.print(NitrousLC1->read(1));
    Serial2.print(NitrousLC1->read(1));
    Serial.print("|LC_N2:");
    Serial2.print("|LC_N2:");
    Serial.print(NitrousLC2->read(1));
    Serial2.print(NitrousLC2->read(1));
    Serial.print("|LC_N3:");
    Serial2.print("|LC_N3:");
    Serial.print(NitrousLC3->read(1));
    Serial2.print(NitrousLC3->read(1));
    Serial.print("|LC_T:");
    Serial2.print("|LC_T:");
    Serial.print(ThrustLC->read());
    Serial2.print(ThrustLC->read());
    //Serial.print("|");
    Serial.println();
    Serial2.println();
  }
}

//===========================EXECUTION============================//
void setup()
{
  Serial.begin(57600);
  Serial2.begin(57600);
  Serial.flush();
  Serial2.flush();
  Serial.println("START");
  Serial2.println("START");
  delay(2000);


  EthaneMBV = new MBV(ETHANE_MBV_PIN, 42, 44);
  NitrousMBV = new MBV(NITROUS_MBV_PIN, 50, 48);

  EthaneLC1 = new LoadCell(22, 23);
  EthaneLC2 = new LoadCell(24, 25);
  EthaneLC3 = new LoadCell(26, 27);
  NitrousLC1 = new LoadCell(34, 35);
  NitrousLC2 = new LoadCell(32, 33);
  NitrousLC3 = new LoadCell(30, 31);

  EthaneLC1->setCalFactor(43.27);
  EthaneLC2->setCalFactor(43.92);
  EthaneLC3->setCalFactor(42.46);
  NitrousLC1->setCalFactor(40.0);
  NitrousLC2->setCalFactor(40.0);
  NitrousLC3->setCalFactor(40.0);

  ThrustLC = new ThrustCell();
  ThrustLC->setCalFactor(-5.83);

  RerouteTC = new Thermocouple(0x67);
  RerouteTC2 = new Thermocouple(0x65);
  RerouteTC3 = new Thermocouple(0x66);

  EthaneHeater1 = new Heater(49, &EthaneUpstreamPT);
  EthaneHeater2 = new Heater(51, &EthaneUpstreamPT);
  NitrousHeater1 = new Heater(49, &NitrousUpstreamPT);
  NitrousHeater2 = new Heater(51, &NitrousUpstreamPT);
  
  delay(1200);
  if (print_data) {
    start_time = millis();
  }
  if(full_output) {
    Serial.println("Commands:");
    Serial2.println("Commands:");
    Serial.println("  ETHANE_VENT_ON / ETHANE_VENT_OFF");
    Serial2.println("  ETHANE_VENT_ON / ETHANE_VENT_OFF");
    Serial.println("  NITROUS_VENT_ON / NITROUS_VENT_OFF");
    Serial2.println("  NITROUS_VENT_ON / NITROUS_VENT_OFF");
    Serial.println("  ETHANE_RUN_ON / ETHANE_RUN_OFF");
    Serial2.println("  ETHANE_RUN_ON / ETHANE_RUN_OFF");
    Serial.println("  NITROUS_RUN_ON / NITROUS_RUN_OFF");
    Serial2.println("  NITROUS_RUN_ON / NITROUS_RUN_OFF");
  }
}

String cmd = "";

void loop()
{
  //Serial.flush();
  elapsed = millis()-start_time;

  // LOG
  if (millis() - lastPressureMs >= log_interval_ms && print_data)
  {
    lastPressureMs += log_interval_ms;
    status(full_output);
  }

  if (elapsed/1000.0 > runtime && runtime >= 0) {
    print_data = false;
  }

  // COMMANDS
  if (Serial.available() || Serial2.available())
  {
    if (Serial.available()) cmd = Serial.readStringUntil('\n');
    else cmd = Serial2.readStringUntil('\n');
    cmd.trim();

    // vent
    if (cmd.equalsIgnoreCase("ETHANE_VENT_ON")) {
      EthaneVent.open();
      if(full_output) { Serial.println("ETHANE VENT OPEN"); Serial2.println("ETHANE VENT OPEN"); }
    }
    else if (cmd.equalsIgnoreCase("ETHANE_VENT_OFF")) {
      EthaneVent.close();
      if(full_output) { Serial.println("ETHANE VENT CLOSED"); Serial2.println("ETHANE VENT CLOSED"); }
    } 
    else if (cmd.equalsIgnoreCase("NITROUS_VENT_ON")) {
      NitrousVent.open();
      if(full_output) { Serial.println("NITROUS VENT OPEN"); Serial2.println("NITROUS VENT OPEN"); }
    }
    else if (cmd.equalsIgnoreCase("NITROUS_VENT_OFF")) {
      NitrousVent.close();
      if(full_output) { Serial.println("NITROUS VENT CLOSED"); Serial2.println("NITROUS VENT CLOSED"); }
    }

    // solenoid
    else if (cmd.equalsIgnoreCase("ETHANE_RUN_ON")) {
      EthaneRunValve.open();
      if(full_output) { Serial.println("ETHANE RUN VALVE OPEN"); Serial2.println("ETHANE RUN VALVE OPEN"); }
    }
    else if (cmd.equalsIgnoreCase("ETHANE_RUN_OFF")) {
      EthaneRunValve.close();
      if(full_output) { Serial.println("ETHANE RUN VALVE CLOSED"); Serial2.println("ETHANE RUN VALVE CLOSED"); }
    } 
    else if (cmd.equalsIgnoreCase("NITROUS_RUN_ON")) {
      NitrousRunValve.open();
      if(full_output) { Serial.println("NITROUS RUN VALVE OPEN"); Serial2.println("NITROUS RUN VALVE OPEN"); }
    } 
    else if (cmd.equalsIgnoreCase("NITROUS_RUN_OFF")) {
      NitrousRunValve.close();
      if(full_output) { Serial.println("NITROUS RUN VALVE CLOSED"); Serial2.println("NITROUS RUN VALVE CLOSED"); }
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
    else if (cmd.equalsIgnoreCase("STATUS")) {
      status();
    }
    
    // ball
    else if (cmd.equalsIgnoreCase("ETHANE_A")) {
      EthaneMBV->next_90();
    } 
    else if (cmd.equalsIgnoreCase("ETHANE_S")) {
      EthaneMBV->move_degrees(10);
    } 
    else if (cmd.equalsIgnoreCase("ETHANE_D")) {
      EthaneMBV->move_degrees(360);
    } 
    else if (cmd.equalsIgnoreCase("ETHANE_R")) {
      EthaneMBV->reset();
      if(full_output) { Serial.println("Position reset to zero"); Serial2.println("Position reset to zero"); }
    }
    else if (cmd.equalsIgnoreCase("ETHANE_STATUS")) {
      EthaneMBV->status();
    }

    else {
      if(full_output) { Serial.println("Unknown command."); Serial2.println("Unknown command."); }
    }
  }

  EthaneMBV->update();
  NitrousMBV->update();

  if(heaters_active) {
    EthaneHeater1->update();
    EthaneHeater2->update();
    NitrousHeater1->update();
    NitrousHeater2->update();
  }

  delay(20);
}