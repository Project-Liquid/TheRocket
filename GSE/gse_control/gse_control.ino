#include "transducer.h"
#include "solenoid.h"
#include "mbv.h"
#include "loadCell.h"
#include "thermocouple.h"
#include "heater.h"
#include "thrustCell.h"
#include "serial.h"
#include <StandardCplusplus.h>
#include <string>
#include "ADS1118.h"
#include <SPI.h>

// Run Control
bool print_data = true;
long start_time = 0; //1410065408; // max value placeholder
const int runtime = -1;
//static unsigned long lastPressureMs = 0;
//const int log_interval_ms = 50;
const bool full_output = false;
long elapsed = 0;
bool static_fire_initializing = false;
long static_fire_duration_ms = 0;

struct LogInterval {
  const unsigned int interval_ms;
  unsigned long last_trigger_ms;
};

LogInterval PTLog{50, 0};
LogInterval LCLog{200, 0};
LogInterval ValveLog{200, 0};
LogInterval MBVLog{200, 0};
LogInterval TCLog{200, 0};

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
const int NITROUS_RUN_PIN = 36;
const int NITROUS_VENT_PIN = 40;

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

const float ETHANE_TARGET_PRESSURE = 90;
const float NITROUS_TARGET_PRESSURE = 90;

// Thermocouple
Thermocouple* RerouteTC;
Thermocouple* RerouteTC2;
Thermocouple* RerouteTC3;
ADS1118* ChamberTC;

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
    SerialDual.print(String(elapsed/1000.0) + "s");
    // PTs
    SerialDual.print("Ethane Upstream: \t" + EthaneUpstreamPT.status());
    SerialDual.print("Ethane Downstream: \t" + EthaneDownstreamPT.status());
    SerialDual.print("Nitrous Upstream: \t" + NitrousUpstreamPT.status());
    SerialDual.print("Nitrous Downstream: \t" + NitrousDownstreamPT.status());
    SerialDual.print("Reroute: \t" + ReroutePT.status());
    // Load Cells
    SerialDual.print("Ethane: \t");
    SerialDual.print("LC1: " + String(EthaneLC1->read(1))); 
    SerialDual.print("\tLC2: " + String(EthaneLC2->read(1)));
    SerialDual.print("\tLC3: " + String(EthaneLC3->read(1))); 
    // SerialDual.print("\tTotal: " + String(SerialDual.print(readEthaneLC()));
    SerialDual.println();
    SerialDual.print("Nitrous: \t");
    SerialDual.print("LC1: " + String(NitrousLC1->read(1))); 
    SerialDual.print("\tLC2: " + String(NitrousLC2->read(1))); 
    SerialDual.print("\tLC3: " + String(NitrousLC3->read(1))); 
    // SerialDual.print("\tTotal: " + String(readNitrousLC()));
    SerialDual.println();
    // SerialDual.print("Thrust: " + String(ThrustLC->read()));
    // SerialDual.println();
    // Tank Heaters
    // SerialDual.print("Ethane Tank 1: " + String(EthaneHeater1->isOn() ? "ON" : "OFF"));
    // SerialDual.print("Nitrous Tank 1: " + String(NitrousHeater1->isOn() ? "ON" : "OFF"));
    // Thermocouple
    // SerialDual.println("Reroute TC: " + String(RerouteTC->readHot()) + "F");
    // SerialDual.println("Reroute TC2: " + String(RerouteTC2->readHot()) + "F");
    // SerialDual.println("Reroute TC3: " + String(RerouteTC3->readHot()) + "F");
    // MBVs
    //SerialDual.println("Ethane MBV: " + String(EthaneMBV->getCurrentDegrees()));
  } else { // Simplified output for log
    String datastream = "";
    long time = millis();
    // FORMAT: DATA|millis|KEY:VAL|KEY:VAL|...
    datastream += "DATA|" + String((elapsed/1000.0), 3);
    if(time - PTLog.last_trigger_ms >= PTLog.interval_ms) {
      datastream += "|PT_EU:" + String(EthaneUpstreamPT.readPressure(), 3);
      datastream += "|PT_ED:" + String(EthaneDownstreamPT.readPressure(), 3);
      datastream += "|PT_NU:" + String(NitrousUpstreamPT.readPressure(), 3);
      datastream += "|PT_ND:" + String(NitrousDownstreamPT.readPressure(), 3);
      PTLog.last_trigger_ms = time - (time % PTLog.interval_ms);
    }
    
    if(time - LCLog.last_trigger_ms >= LCLog.interval_ms) {
      datastream += "|LC_E1:" + String(EthaneLC1->read(1));
      datastream += "|LC_E2:" + String(EthaneLC2->read(1));
      datastream += "|LC_E3:" + String(EthaneLC3->read(1));
      datastream += "|LC_N1:" + String(NitrousLC1->read(1));
      datastream += "|LC_N2:" + String(NitrousLC2->read(1));
      datastream += "|LC_N3:" + String(NitrousLC3->read(1));
      datastream += "|LC_T:" + String(ThrustLC->read());
      LCLog.last_trigger_ms = time - (time % LCLog.interval_ms);
    }

    if(time - ValveLog.last_trigger_ms >= ValveLog.interval_ms) {
      datastream += "|ERV:" + String(EthaneRunValve.state());
      datastream += "|EV:" + String(EthaneVent.state());
      datastream += "|NRV:" + String(NitrousRunValve.state());
      datastream += "|NV:" + String(NitrousVent.state());
      datastream += "|MBV_E:" + String(EthaneMBV->getCurrentDegrees());
      datastream += "|MBV_N:" + String(NitrousMBV->getCurrentDegrees());
      ValveLog.last_trigger_ms = time - (time % ValveLog.interval_ms);
    }

    if(time - TCLog.last_trigger_ms >= TCLog.interval_ms) {
      // datastream += "|TC_C:" + String(ChamberTC->getTemperature());
      datastream += "|TC_R:" + String(RerouteTC->readHot());
      TCLog.last_trigger_ms = time - (time % TCLog.interval_ms);
    }

    if (datastream.lastIndexOf("|") > 4) SerialDual.println(datastream);
  }
}

void EMERGENCY_STOP() {
  if (EthaneMBV->isOpen()) {
    EthaneMBV->next_90();
  }
  if (NitrousMBV->isOpen()) {
    NitrousMBV->next_90();
  }
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

void staticFire() {
  //NitrousRunValve.open();
  EthaneRunValve.open();

  //Ignitor.setNextActutation(5000);
  delay(1000); // REMOVE REMOVE REMOVE
  //if (ChamberTC->getTemperature() > 100) {
    //NitrousMBV->next_90();
    EthaneMBV->setNextActuation(200);
    //NitrousMBV->setNextActuation(200 + static_fire_duration_ms);
    EthaneMBV->setNextActuation(200 + static_fire_duration_ms + 500);
    //NitrousRunValve.setNextActuation(200 + static_fire_duration_ms+600, false);
    EthaneRunValve.setNextActuation(200 + static_fire_duration_ms + 600, false);
    
    // Vent
    EthaneVent.setNextActuation(800 + static_fire_duration_ms + 1000, true);
    EthaneVent.setNextActuation(800 + static_fire_duration_ms + 4000, false);
    //NitrousVent.setNextActuation(800 + static_fire_duration_ms + 5000, true);
    //NitrousVent.setNextActuation(800 + static_fire_duration_ms + 9000, false);

    static_fire_initializing = false;
  //}
}

//===========================EXECUTION============================//
void setup()
{
  // TURN SERIAL OFF; TURN SERIAL2 ON
  SerialDual.setActive(false, true);
  SerialDual.begin(57600);
  //SerialDual.flush();
  SerialDual.println("START");
  delay(2000);

  EthaneMBV = new MBV(ETHANE_MBV_PIN, 44, 42);
  NitrousMBV = new MBV(NITROUS_MBV_PIN, 38, 48);

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

  pinMode(53, OUTPUT);      // Force SS high to lock Mega in master mode
  digitalWrite(53, HIGH);
  ChamberTC = new ADS1118(53);
  delay(100);

  ChamberTC->setSamplingRate(ChamberTC->RATE_860SPS);
  ChamberTC->setInputSelected(ChamberTC->DIFF_0_1);
  ChamberTC->setFullScaleRange(ChamberTC->FSR_0256);

  EthaneHeater1 = new Heater(49, &EthaneUpstreamPT);
  EthaneHeater2 = new Heater(41, &EthaneUpstreamPT);
  NitrousHeater1 = new Heater(43, &NitrousUpstreamPT);
  NitrousHeater2 = new Heater(45, &NitrousUpstreamPT);

  EthaneHeater1->setTarget(ETHANE_TARGET_PRESSURE);
  EthaneHeater2->setTarget(ETHANE_TARGET_PRESSURE);
  NitrousHeater1->setTarget(NITROUS_TARGET_PRESSURE);
  NitrousHeater2->setTarget(NITROUS_TARGET_PRESSURE);

  EthaneUpstreamPT.setRedline(1100, 25);
  NitrousUpstreamPT.setRedline(1100, 25);
  
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
  elapsed = millis()-start_time;

  EthaneMBV->update();
  NitrousMBV->update();

  // LOG
  if (print_data)
  {
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
    cmd.toUpperCase();

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
    } 
    else if (cmd.equalsIgnoreCase("STOP")) {
      print_data = false;
    }
    else if (cmd.equalsIgnoreCase("STATUS")) {
      status();
    }
    
    // ball
    else if (cmd.indexOf("ETHANE_MBV_") >= 0) {
      String degrees = cmd.substring(11);
      if(degrees.indexOf("N") >= 0) {
        EthaneMBV->next_90();
      } else if(degrees.toInt()) {
        EthaneMBV->move_degrees(degrees.toInt());
      }
    }
    else if (cmd.equalsIgnoreCase("ETHANE_MBV_RESET")) {
      EthaneMBV->reset();
      if(full_output) { SerialDual.println("Ethane MBV position reset to zero"); }
    }
    else if (cmd.equalsIgnoreCase("ETHANE_MBV_STATUS")) {
      EthaneMBV->status();
    }
    else if (cmd.indexOf("NITROUS_MBV_") >= 0) {
      String degrees = cmd.substring(12);
      if(degrees.indexOf("N") >= 0) {
        NitrousMBV->next_90();
      } else if(degrees.toInt()) {
        NitrousMBV->move_degrees(degrees.toInt());
      }
    }
    else if (cmd.equalsIgnoreCase("NITROUS_MBV_RESET")) {
      NitrousMBV->reset();
      if(full_output) { SerialDual.println("Nitrous MBV position reset to zero"); }
    }
    else if (cmd.equalsIgnoreCase("NITROUS_MBV_STATUS")) {
      NitrousMBV->status();
    }
    else if (cmd.equalsIgnoreCase("E_STOP")) {
      EMERGENCY_STOP();
    }
    else if (cmd.indexOf("COLD_FLOW_") >= 0) {
      String prop = cmd.substring(10, cmd.lastIndexOf("_"));
      String duration = cmd.substring(cmd.lastIndexOf("_") + 1);
      if (duration.toInt()) {
        if (prop.equalsIgnoreCase("ETHANE")) {
          coldFlowEthane(duration.toInt());
        } else if (prop.equalsIgnoreCase("NITROUS")) {
          coldFlowNitrous(duration.toInt());
        }
      }
    }
    else if (cmd.indexOf("STATIC_FIRE_") >= 0) {
      String duration = cmd.substring(cmd.lastIndexOf("_") + 1);
      if (duration.toInt()) {
        static_fire_duration_ms = duration.toInt();
        static_fire_initializing = true;
      }
    }

    else {
      if(full_output) { SerialDual.println("Unknown command."); }
    }
  }

  // Scheduling Updates
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

  // Redlines
  if (EthaneUpstreamPT.checkRedline()) {
    SerialDual.println("CRITICAL ERROR: ETHANE PRESSURE REDLINE EXCEEDED");
    if (EthaneMBV->isOpen()) {
      EthaneMBV->next_90();
    }
    EthaneVent.open();
    EthaneVent.setNextActuation(10000, false);
  }

  if (NitrousUpstreamPT.checkRedline()) {
    SerialDual.println("CRITICAL ERROR: NITROUS PRESSURE REDLINE EXCEEDED");
    if (NitrousMBV->isOpen()) {
      NitrousMBV->next_90();
    }
    NitrousVent.open();
    NitrousVent.setNextActuation(10000, false);
  }

  if(static_fire_initializing) {
    staticFire();
  }

  //delay(10);
}