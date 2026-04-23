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
long elapsed = 0;

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

MBV* EthaneMBV;
MBV* NitrousMBV;

// Load Cells
LoadCell* EthaneLC2;
LoadCell* EthaneLC1;
LoadCell* EthaneLC3;
LoadCell* NitrousLC1;
LoadCell* NitrousLC2;
LoadCell* NitrousLC3;

// ThrustCell* ThrustLC;

// Tank Heaters
// bool heaters_active = false;
// Heater* Heater1;
// Heater* Heater2;

// Thermocouple
// Thermocouple* RerouteTC;


//===========================FUNCTIONS============================//

float readEthaneLC() {
  return EthaneLC1->read(1) + EthaneLC2->read(1) + EthaneLC3->read(1);
}

float readNitrousLC() {
  return NitrousLC1->read(1) + NitrousLC2->read(1) + NitrousLC3->read(1);
}

void status() {
  elapsed = millis()-start_time;

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
  Serial.print("LC1: "); Serial.print(EthaneLC1->read(1)); 
  Serial.print("\tLC2: "); Serial.print(EthaneLC2->read(1));
  Serial.print("\tLC3: "); Serial.print(EthaneLC3->read(1)); 
  // Serial.print("\tTotal: "); Serial.print(readEthaneLC());
  Serial.println();
  Serial.print("Nitrous: \t");
  Serial.print("LC1: "); Serial.print(NitrousLC1->read(1)); 
  Serial.print("\tLC2: "); Serial.print(NitrousLC2->read(1)); 
  Serial.print("\tLC3: "); Serial.print(NitrousLC3->read(1)); 
  Serial.print("\tTotal: "); Serial.print(readNitrousLC());
  Serial.println();
  // Serial.print("Thrust: "); Serial.print(ThrustLC->getAverageReading(1));
  // Serial.println();
  // Tank Heaters
  // Serial.print("Tank 1: "); 
  // Serial.print(Heater1->getTemp()); Serial.print("F -- "); 
  // Serial.println(Heater1->isOn() ? "ON" : "OFF");
  // Serial.print("Tank 2: "); 
  // Serial.print(Heater2->getTemp()); Serial.print("F -- "); 
  // Serial.println(Heater2->isOn() ? "ON" : "OFF");
  // Thermocouple
  //Serial.print("Reroute TC: ");
  //Serial.print(RerouteTC->readHot()); Serial.print("F");
  Serial.println();
}

//===========================EXECUTION============================//
void setup()
{
  Serial.begin(9600);
  Serial.flush();
  Serial.println("START");
  delay(2000);


  EthaneMBV = new MBV(ETHANE_MBV_PIN, 42, 44);
  NitrousMBV = new MBV(NITROUS_MBV_PIN, 50, 48);

  EthaneLC1 = new LoadCell(24, 23);
  EthaneLC2 = new LoadCell(22, 21);
  EthaneLC3 = new LoadCell(26, 25);
  NitrousLC1 = new LoadCell(34, 35);
  NitrousLC2 = new LoadCell(32, 33);
  NitrousLC3 = new LoadCell(30, 31);

  EthaneLC1->setCalFactor(39.32);
  EthaneLC2->setCalFactor(43.89);
  EthaneLC3->setCalFactor(47.28);
  NitrousLC1->setCalFactor(40.0);
  NitrousLC2->setCalFactor(60.0);
  NitrousLC3->setCalFactor(60.0);

  // ThrustLC = new ThrustCell();
  // ThrustLC->setCalFactor(5.83);

  // RerouteTC = new Thermocouple(0x67);

  // Heater1 = new Heater(49, 0x66);
  // Heater2 = new Heater(51, 0x65);
  
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
  elapsed = millis()-start_time;

  // LOG
  if (millis() - lastPressureMs >= log_interval_ms && print_data)
  {
    lastPressureMs += log_interval_ms;
    if(full_output) {
      status();
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
      Serial.print(EthaneLC1->read(1));
      Serial.print(", ");
      Serial.print(EthaneLC2->read(1));
      Serial.print(", ");
      Serial.print(EthaneLC3->read(1));
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
      if(full_output) Serial.println("Position reset to zero");
    }
    else if (cmd.equalsIgnoreCase("ETHANE_STATUS")) {
      EthaneMBV->status();
    }

    else {
      if(full_output) Serial.println("Unknown command.");
    }
  }

  EthaneMBV->update();
  NitrousMBV->update();

  // if(heaters_active) {
  //   Heater1->update();
  //   Heater2->update();
  // }

  delay(20);
}