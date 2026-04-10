#include "transducer.h"
#include "solenoid.h"
#include "mbv.h"

const int ETHANE_UPSTREAM_PIN = A8;
const int ETHANE_DOWNSTREAM_PIN = A9;
const int NITROUS_UPSTREAM_PIN = A6;
const int NITROUS_DOWNSTREAM_PIN = A7;
const float V_REF = 1.1;
const float R_SHUNT = 150.0;

const float P_MIN = 0.0;
const float P_MAX_ETHANE = 68.9;
const float P_MAX_NITROUS = 103.4214;

const int ETHANE_RUN_PIN = 48;
const int ETHANE_VENT_PIN = 50;
const int NITROUS_RUN_PIN = 22;
const int NITROUS_VENT_PIN = 24;

bool print_data = false;
long start_time = 1410065408; // max value placeholder
static unsigned long lastPressureMs = 0;
const int runtime = (10*60) + 8;
const int log_interval_ms = 500;

const bool full_output = false;

Transducer EthaneUpstreamPT = Transducer(ETHANE_UPSTREAM_PIN, P_MIN, P_MAX_ETHANE);
Transducer EthaneDownstreamPT = Transducer(ETHANE_DOWNSTREAM_PIN, P_MIN, P_MAX_ETHANE);
Transducer NitrousUpstreamPT = Transducer(NITROUS_UPSTREAM_PIN, P_MIN, P_MAX_NITROUS);
Transducer NitrousDownstreamPT = Transducer(NITROUS_DOWNSTREAM_PIN, P_MIN, P_MAX_NITROUS);

Solenoid EthaneRunValve = Solenoid(ETHANE_RUN_PIN);
Solenoid EthaneVent = Solenoid(ETHANE_VENT_PIN);
Solenoid NitrousRunValve = Solenoid(NITROUS_RUN_PIN);
Solenoid NitrousVent = Solenoid(NITROUS_VENT_PIN);

const int ETHANE_MBV_PIN = 52;
const int NITROUS_MBV_PIN = 26;

MBV EthaneMBV = MBV(ETHANE_MBV_PIN, 38, 36);
MBV NitrousMBV = MBV(NITROUS_MBV_PIN, 30, 32);

void setup()
{
  Serial.begin(9600);
  analogReference(INTERNAL1V1);

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
      EthaneMBV.move_90();
    } 
    else if (cmd.equalsIgnoreCase("ETHANE_S")) {
      EthaneMBV.move_small();
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

  delay(10);
}