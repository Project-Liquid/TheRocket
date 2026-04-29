#include "transducer.h"
#include "solenoid.h"
#include "mbv.h"
#include "loadCell.h"
#include "thermocouple.h"
#include "heater.h"
#include "thrustCell.h"
#include "serial.h"

// Run Control
bool print_data = true;
long start_time = 1410065408; // max value placeholder
static unsigned long lastPressureMs = 0;
const int runtime = -1;
const int log_interval_ms = 50;
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

// Serial Communication for Radios
SerialDualClass SerialDual(Serial, Serial2);

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
    SerialDual.println("====================");
    // Time
    SerialDual.print(elapsed/1000.0); SerialDual.println("s");
    // PTs
    SerialDual.print("Ethane Upstream: \t");
    EthaneUpstreamPT.status();
    SerialDual.print("Ethane Downstream: \t");
    EthaneDownstreamPT.status();
    SerialDual.print("Nitrous Upstream: \t");
    NitrousUpstreamPT.status();
    SerialDual.print("Nitrous Downstream: \t");
    NitrousDownstreamPT.status();
    SerialDual.print("Reroute: \t");
    ReroutePT.status();
    // Load Cells
    SerialDual.print("Ethane: \t");
    SerialDual.print("LC1: "); SerialDual.print(EthaneLC1->read(1)); 
    SerialDual.print("\tLC2: "); SerialDual.print(EthaneLC2->read(1));
    SerialDual.print("\tLC3: "); SerialDual.print(EthaneLC3->read(1)); 
    // SerialDual.print("\tTotal: "); SerialDual.print(readEthaneLC());
    SerialDual.println();
    SerialDual.print("Nitrous: \t");
    SerialDual.print("LC1: "); SerialDual.print(NitrousLC1->read(1)); 
    SerialDual.print("\tLC2: "); SerialDual.print(NitrousLC2->read(1)); 
    SerialDual.print("\tLC3: "); SerialDual.print(NitrousLC3->read(1)); 
    // SerialDual.print("\tTotal: "); SerialDual.print(readNitrousLC());
    SerialDual.println();
    // SerialDual.print("Thrust: "); SerialDual.print(ThrustLC->read());
    // SerialDual.println();
    // SerialDual.println();
    // Tank Heaters
    // SerialDual.print("Ethane Tank 1: "); 
    // SerialDual.println(EthaneHeater1->isOn() ? "ON" : "OFF");
    // SerialDual.print("Nitrous Tank 1: "); 
    // SerialDual.println(NitrousHeater1->isOn() ? "ON" : "OFF");
    // Thermocouple
    // SerialDual.print("Reroute TC: ");
    // SerialDual.print(RerouteTC->readHot()); SerialDual.print("F");
    // SerialDual.println();
    // SerialDual.print("Reroute TC2: ");
    // SerialDual.print(RerouteTC2->readHot()); SerialDual.print("F");
    // SerialDual.println();
    // SerialDual.print("Reroute TC3: ");
    // SerialDual.print(RerouteTC3->readHot()); SerialDual.print("F");
    // SerialDual.println();
    // MBVs
    //SerialDual.print("Ethane MBV: "); SerialDual.println(EthaneMBV->getCurrentDegrees());
  } else { // Simplified output for log
    // FORMAT: DATA|millis|KEY:VAL|KEY:VAL|...
    SerialDual.print("DATA|");
    SerialDual.print(elapsed/1000.0);
    SerialDual.print("|PT_EU:");
    EthaneUpstreamPT.value();
    SerialDual.print("|PT_ED:");
    EthaneDownstreamPT.value();
    // SerialDual.print("|PT_NU:");
    // NitrousUpstreamPT.value();
    // SerialDual.print("|PT_ND:");
    // NitrousDownstreamPT.value();
    // SerialDual.print("|LC_E1:");
    // SerialDual.print(EthaneLC1->read(1));
    // SerialDual.print("|LC_E2:");
    // SerialDual.print(EthaneLC2->read(1));
    // SerialDual.print("|LC_E3:");
    // SerialDual.print(EthaneLC3->read(1));
    // SerialDual.print("|LC_N1:");
    // SerialDual.print(NitrousLC1->read(1));
    // SerialDual.print("|LC_N2:");
    // SerialDual.print(NitrousLC2->read(1));
    // SerialDual.print("|LC_N3:");
    // SerialDual.print(NitrousLC3->read(1));
    // SerialDual.print("|LC_T:");
    // SerialDual.print(ThrustLC->read());
    // SerialDual.print("|ERV:");
    // SerialDual.print(EthaneRunValve.state());
    // SerialDual.print("|EV:");
    // SerialDual.print(EthaneVent.state());
    // SerialDual.print("|NRV:");
    // SerialDual.print(NitrousRunValve.state());
    // SerialDual.print("|NV:");
    // SerialDual.print(NitrousVent.state());
    // SerialDual.print("|MBV_E:");
    // SerialDual.print(EthaneMBV->getCurrentDegrees());
    // SerialDual.print("|MBV_N:");
    // SerialDual.print(NitrousMBV->getCurrentDegrees());
    // //SerialDual.print("|");
    SerialDual.println();
  }
}

void EMERGENCY_STOP() {
  ventEthane();
  ventNitrous();
}

void ventEthane() {
  EthaneRunValve.close();
  EthaneVent.setNextActuation(1000, true);
  EthaneVent.setNextActuation(4000, false);
}

void ventNitrous() {
  NitrousRunValve.close();
  NitrousVent.setNextActuation(1000, true);
  NitrousVent.setNextActuation(4000, false);
}

void coldFlowEthane(long duration_ms) {
  EthaneRunValve.clearSchedule();
  EthaneMBV->clearSchedule();

  EthaneRunValve.open();

  EthaneMBV->setNextActuation(5000);
  EthaneMBV->setNextActuation(5000 + duration_ms);

  EthaneRunValve.setNextActuation(10000 + duration_ms, false);
}

void coldFlowNitrous(long duration_ms) {
  NitrousRunValve.open();

  NitrousMBV->setNextActuation(5000);
  NitrousMBV->setNextActuation(5000 + duration_ms);

  NitrousRunValve.setNextActuation(10000 + duration_ms, false);
}

//===========================EXECUTION============================//
void setup()
{
  SerialDual.begin(57600);
  //SerialDual.flush();
  SerialDual.println("START");
  delay(2000);


  EthaneMBV = new MBV(ETHANE_MBV_PIN, 44, 42);
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
    SerialDual.println("Commands:");
    SerialDual.println("  ETHANE_VENT_ON / ETHANE_VENT_OFF");
    SerialDual.println("  NITROUS_VENT_ON / NITROUS_VENT_OFF");
    SerialDual.println("  ETHANE_RUN_ON / ETHANE_RUN_OFF");
    SerialDual.println("  NITROUS_RUN_ON / NITROUS_RUN_OFF");
  }
}

String cmd = "";

void loop()
{
  //SerialDual.flush();
  elapsed = millis()-start_time;

  EthaneMBV->update();
  NitrousMBV->update();

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
      if(full_output) { SerialDual.println("ETHANE VENT OPEN"); }
    }
    else if (cmd.equalsIgnoreCase("ETHANE_VENT_OFF")) {
      EthaneVent.close();
      if(full_output) { SerialDual.println("ETHANE VENT CLOSED"); }
    } 
    else if (cmd.equalsIgnoreCase("NITROUS_VENT_ON")) {
      NitrousVent.open();
      if(full_output) { SerialDual.println("NITROUS VENT OPEN"); }
    }
    else if (cmd.equalsIgnoreCase("NITROUS_VENT_OFF")) {
      NitrousVent.close();
      if(full_output) { SerialDual.println("NITROUS VENT CLOSED"); }
    }

    // solenoid
    else if (cmd.equalsIgnoreCase("ETHANE_RUN_ON")) {
      EthaneRunValve.open();
      if(full_output) { SerialDual.println("ETHANE RUN VALVE OPEN"); }
    }
    else if (cmd.equalsIgnoreCase("ETHANE_RUN_OFF")) {
      EthaneRunValve.close();
      if(full_output) { SerialDual.println("ETHANE RUN VALVE CLOSED"); }
    } 
    else if (cmd.equalsIgnoreCase("NITROUS_RUN_ON")) {
      NitrousRunValve.open();
      if(full_output) { SerialDual.println("NITROUS RUN VALVE OPEN"); }
    } 
    else if (cmd.equalsIgnoreCase("NITROUS_RUN_OFF")) {
      NitrousRunValve.close();
      if(full_output) { SerialDual.println("NITROUS RUN VALVE CLOSED"); }
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
    else if (cmd.equalsIgnoreCase("ETHANE_90")) {
      EthaneMBV->next_90();
    } 
    else if (cmd.equalsIgnoreCase("ETHANE_10")) {
      EthaneMBV->move_degrees(10);
    } 
    else if (cmd.equalsIgnoreCase("ETHANE_360")) {
      EthaneMBV->move_degrees(360);
    } 
    else if (cmd.equalsIgnoreCase("ETHANE_RESET")) {
      EthaneMBV->reset();
      if(full_output) { SerialDual.println("Ethane MBV position reset to zero"); }
    }
    else if (cmd.equalsIgnoreCase("ETHANE_MBV_STATUS")) {
      EthaneMBV->status();
    }
    else if (cmd.equalsIgnoreCase("NITROUS_90")) {
      NitrousMBV->next_90();
    } 
    else if (cmd.equalsIgnoreCase("NITROUS_10")) {
      NitrousMBV->move_degrees(10);
    } 
    else if (cmd.equalsIgnoreCase("NITROUS_360")) {
      NitrousMBV->move_degrees(360);
    } 
    else if (cmd.equalsIgnoreCase("NITROUS_RESET")) {
      NitrousMBV->reset();
      if(full_output) { SerialDual.println("Nitrous MBV position reset to zero"); }
    }
    else if (cmd.equalsIgnoreCase("NITROUS_MBV_STATUS")) {
      NitrousMBV->status();
    }
    else if (cmd.equalsIgnoreCase("E_STOP")) {
      EMERGENCY_STOP();
    }
    else if (cmd.equalsIgnoreCase("COLD_FLOW")) {
      coldFlowEthane(5000);
    }

    else {
      if(full_output) { SerialDual.println("Unknown command."); }
    }
  }

  // EthaneMBV->update();
  // NitrousMBV->update();
  EthaneVent.checkScheduledActuation();
  NitrousVent.checkScheduledActuation();
  EthaneRunValve.checkScheduledActuation();
  NitrousRunValve.checkScheduledActuation();
  EthaneMBV->checkScheduledActuation();
  NitrousMBV->checkScheduledActuation();

  if(heaters_active) {
    EthaneHeater1->update();
    EthaneHeater2->update();
    NitrousHeater1->update();
    NitrousHeater2->update();
  }

  //delay(10);
}