#include "main.h"

bool print_data = true;
unsigned long start_time = 0; 
const int runtime = -1;
TransmissionType transmissionFormat = COMPRESSED;
unsigned long time_absolute = 0;
unsigned long time_elapsed = 0;
bool static_fire_initializing = false;
unsigned long static_fire_duration_ms = 0;
bool heaters_active = false;
float ETHANE_TARGET_PRESSURE = 0;
float NITROUS_TARGET_PRESSURE = 0;
PollInterval PTLog{50, 0};
PollInterval LCLog{200, 0};
PollInterval ValveLog{200, 0};
PollInterval MBVLog{200, 0};
PollInterval TCLog{200, 0};
PollInterval RedlinePoll{50, 0};
PollInterval ValveSchedulePoll{50, 0};
PollInterval HeaterPoll{500, 0};
float ETHANE_WEIGHT_REDLINE = -27;
float NITROUS_WEIGHT_REDLINE = 19.5;
int current_highest_redline = 0;
String cmd = "";
//======================OBJECT DEFNINTIONS=======================//

// Pressure Transducers
Transducer EthaneUpstreamPT(ETHANE_UPSTREAM_PIN, P_MIN, P_MAX_ETHANE);
Transducer EthaneDownstreamPT(ETHANE_DOWNSTREAM_PIN, P_MIN, P_MAX_ETHANE);
Transducer NitrousUpstreamPT(NITROUS_UPSTREAM_PIN, P_MIN, P_MAX_NITROUS);
Transducer NitrousDownstreamPT(NITROUS_DOWNSTREAM_PIN, P_MIN, P_MAX_NITROUS);
Transducer ReroutePT(REROUTE_PT_PIN, P_MIN, P_MAX_NITROUS);
Transducer ChamberPT(CHAMBER_PT_PIN, P_MIN, P_MAX_CHAMBER);
Relay EthaneRunValve(ETHANE_RUN_PIN);
Relay EthaneVent(ETHANE_VENT_PIN);
Relay NitrousRunValve(NITROUS_RUN_PIN);
Relay NitrousVent(NITROUS_VENT_PIN);
Relay Ignitor(IGNITOR_PIN);
SerialDualClass SerialDual(Serial, Serial2);
CmdBuffer serialBuf;
CmdBuffer serial2Buf;

MBV* EthaneMBV = nullptr;
MBV* NitrousMBV = nullptr;
LoadCell* EthaneLC2 = nullptr;
LoadCell* EthaneLC1 = nullptr;
LoadCell* EthaneLC3 = nullptr;
LoadCell* NitrousLC1 = nullptr;
LoadCell* NitrousLC2 = nullptr;
LoadCell* NitrousLC3 = nullptr;
ThrustCell* ThrustLC = nullptr;
Heater* EthaneHeater1 = nullptr;
Heater* EthaneHeater2 = nullptr;
Heater* NitrousHeater1 = nullptr;
Heater* NitrousHeater2 = nullptr;
// Thermocouple* RerouteTC = nullptr;
// Thermocouple* RerouteTC2 = nullptr;
// Thermocouple* RerouteTC3 = nullptr;
// ADS1118* ChamberTC = nullptr;

Redline EthaneOverpressure(EthaneOverpressureCondition, EthaneOverpressureResponse, ETHANE_OVERPRESSURE_PRIORITY, OVERPRESSURE_COUNTS_THRESHOLD);
Redline NitrousOverpressure(NitrousOverpressureCondition, NitrousOverpressureResponse, NITROUS_OVERPRESSURE_PRIORITY, OVERPRESSURE_COUNTS_THRESHOLD);
Redline EthaneMBVCloseFailure(EthaneMBVCloseFailureCondition, EthaneMBVCloseFailureResponse, ETHANE_MBV_CLOSE_FAILURE_PRIORITY, MBV_FAILURE_COUNTS_THRESHOLD);
Redline NitrousMBVCloseFailure(NitrousMBVCloseFailureCondition, NitrousMBVCloseFailureResponse, NITROUS_MBV_CLOSE_FAILURE_PRIORITY, MBV_FAILURE_COUNTS_THRESHOLD);
Redline CombustionPropogation(CombustionPropogationCondition, CombustionPropogationResponse, COMBUSTION_PROPOGATION_PRIORITY, COMBUSTION_PROPOGATION_COUNTS_THRESHOLD);
Redline InlineThermalDecomp(InlineThermalDecompCondition, InlineThermalDecompResponse, INLINE_THERMAL_DECOMP_PRIORITY, INLINE_THERMAL_DECOMP_COUNTS_THRESHOLD);
Redline LostLoadCell(LostLoadCellCondition, LostLoadCellResponse, LOST_LOAD_CELL_PRIORITY, LOST_LOAD_CELL_COUNTS_THRESHOLD);
Redline EthaneUnderweight(EthaneUnderweightCondition, EthaneUnderweightResponse, ETHANE_UNDERWEIGHT_PRIORITY, UNDERWEIGHT_COUNTS_THRESHOLD);
Redline NitrousUnderweight(NitrousUnderweightCondition, NitrousUnderweightResponse, NITROUS_UNDERWEIGHT_PRIORITY, UNDERWEIGHT_COUNTS_THRESHOLD);
Redline EthaneOverweight(EthaneOverweightCondition, EthaneOverweightResponse, ETHANE_OVERWEIGHT_PRIORITY, OVERWEIGHT_COUNTS_THRESHOLD);
Redline NitrousOverweight(NitrousOverweightCondition, NitrousOverweightResponse, NITROUS_OVERWEIGHT_PRIORITY, OVERWEIGHT_COUNTS_THRESHOLD);
Redline Redlines[11] = {EthaneOverpressure, NitrousOverpressure, EthaneMBVCloseFailure, NitrousMBVCloseFailure, 
  CombustionPropogation, InlineThermalDecomp, LostLoadCell, EthaneUnderweight, NitrousUnderweight, EthaneOverweight, NitrousOverweight};

//===========================FUNCTIONS============================//
/** MARK: Helpers */

/**
 * Our output function.
 * Note: Some features commented out to speed running.
 */
void status(TransmissionType format = COMPRESSED) {
  time_elapsed = millis()-start_time;
  bool print_current_poll = false;

  switch (format) {
    case DIAGNOSTIC: {
      SerialDual.println("====================");
      // Time
      SerialDual.print(String(time_elapsed/1000.0) + "s");
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
      SerialDual.println();
      SerialDual.print("Nitrous: \t");
      SerialDual.print("LC1: " + String(NitrousLC1->read(1))); 
      SerialDual.print("\tLC2: " + String(NitrousLC2->read(1))); 
      SerialDual.print("\tLC3: " + String(NitrousLC3->read(1))); 
      SerialDual.println();
      SerialDual.print("Thrust: " + String(ThrustLC->read()));
      SerialDual.println();
      // Tank Heaters
      SerialDual.print("Ethane Tank 1: " + String(EthaneHeater1->isOn() ? "ON" : "OFF"));
      SerialDual.print("Nitrous Tank 1: " + String(NitrousHeater1->isOn() ? "ON" : "OFF"));
      // Thermocouple
      // SerialDual.println("Reroute TC: " + String(RerouteTC->readHot()) + "F");
      // SerialDual.println("Reroute TC2: " + String(RerouteTC2->readHot()) + "F");
      // SerialDual.println("Reroute TC3: " + String(RerouteTC3->readHot()) + "F");
      // MBVs
      SerialDual.println("Ethane MBV: " + String(EthaneMBV->getCurrentDegrees()));
      break;
    }

    case COMPRESSED: {
      // FORMAT: DATA|millis|KEY:VAL|KEY:VAL|...
      if(time_absolute - PTLog.last_trigger_ms >= PTLog.interval_ms) {
        SerialDual.print("DATA|"); SerialDual.print((time_elapsed/1000.0), 3);
        SerialDual.print("|PT_EU:");  SerialDual.print(EthaneUpstreamPT.readPressure(), 3);
        SerialDual.print("|PT_ED:");  SerialDual.print(EthaneDownstreamPT.readPressure(), 3);
        SerialDual.print("|PT_NU:");  SerialDual.print(NitrousUpstreamPT.readPressure(), 3);
        SerialDual.print("|PT_ND:");  SerialDual.print(NitrousDownstreamPT.readPressure(), 3);
        SerialDual.print("|PT_RR:");  SerialDual.print(ReroutePT.readPressure(), 3);
        SerialDual.print("|PT_CH:");  SerialDual.print(ChamberPT.readPressure(), 3);
        PTLog.last_trigger_ms = time_absolute - (time_absolute % PTLog.interval_ms);
        print_current_poll = true;
      }
      
      if(time_absolute - LCLog.last_trigger_ms >= LCLog.interval_ms) {
        if (!print_current_poll) { SerialDual.print("DATA|"); SerialDual.print((time_elapsed/1000.0), 3); }
        SerialDual.print("|LC_E1:");  SerialDual.print(EthaneLC1->read(1));
        SerialDual.print("|LC_E2:");  SerialDual.print(EthaneLC2->read(1));
        SerialDual.print("|LC_E3:");  SerialDual.print(EthaneLC3->read(1));
        SerialDual.print("|LC_N1:");  SerialDual.print(NitrousLC1->read(1));
        SerialDual.print("|LC_N2:");  SerialDual.print(NitrousLC2->read(1));
        SerialDual.print("|LC_N3:");  SerialDual.print(NitrousLC3->read(1));
        SerialDual.print("|LC_T:");   SerialDual.print(ThrustLC->read());
        LCLog.last_trigger_ms = time_absolute - (time_absolute % LCLog.interval_ms);
        print_current_poll = true;
      }

      if(time_absolute - ValveLog.last_trigger_ms >= ValveLog.interval_ms) {
        if (!print_current_poll) { SerialDual.print("DATA|"); SerialDual.print((time_elapsed/1000.0), 3); }
        SerialDual.print("|ERV:");    SerialDual.print(EthaneRunValve.state());
        SerialDual.print("|EV:");     SerialDual.print(EthaneVent.state());
        SerialDual.print("|NRV:");    SerialDual.print(NitrousRunValve.state());
        SerialDual.print("|NV:");     SerialDual.print(NitrousVent.state());
        SerialDual.print("|MBV_E:");  SerialDual.print(EthaneMBV->getCurrentDegrees());
        SerialDual.print("|MBV_N:");  SerialDual.print(NitrousMBV->getCurrentDegrees());
        ValveLog.last_trigger_ms = time_absolute - (time_absolute % ValveLog.interval_ms);
        print_current_poll = true;
      }

      if(time_absolute - TCLog.last_trigger_ms >= TCLog.interval_ms) {
        if (!print_current_poll) { SerialDual.print("DATA|"); SerialDual.print((time_elapsed/1000.0), 3); }
        // SerialDual.print("|TC_C:"); SerialDual.print(Thermocouple::cToF(ChamberTC->getTemperature()));
        // SerialDual.print("|TC_R:"); SerialDual.print(RerouteTC->readHot());
        TCLog.last_trigger_ms = time_absolute - (time_absolute % TCLog.interval_ms);
        print_current_poll = true;
      }

      if (print_current_poll) { SerialDual.println(); }
      break;
    }

    case RAW: {
      // FORMAT: DATA|millis|KEY:VAL|KEY:VAL|...
      if(time_absolute - PTLog.last_trigger_ms >= PTLog.interval_ms) {
        SerialDual.print("DATA|"); SerialDual.print((time_elapsed/1000.0), 3);
        SerialDual.print("|PT_EU:");  SerialDual.print(EthaneUpstreamPT.readPressure(), 3);
        SerialDual.print("|PT_ED:");  SerialDual.print(EthaneDownstreamPT.readPressure(), 3);
        SerialDual.print("|PT_NU:");  SerialDual.print(NitrousUpstreamPT.readPressure(), 3);
        SerialDual.print("|PT_ND:");  SerialDual.print(NitrousDownstreamPT.readPressure(), 3);
        PTLog.last_trigger_ms = time_absolute - (time_absolute % PTLog.interval_ms);
        print_current_poll = true;
      }
      
      if(time_absolute - LCLog.last_trigger_ms >= LCLog.interval_ms) {
        if (!print_current_poll) { SerialDual.print("DATA|"); SerialDual.print((time_elapsed/1000.0), 3); }
        SerialDual.print("|LC_E1:");  SerialDual.print(EthaneLC1->read(1));
        SerialDual.print("|LC_E2:");  SerialDual.print(EthaneLC2->read(1));
        SerialDual.print("|LC_E3:");  SerialDual.print(EthaneLC3->read(1));
        SerialDual.print("|LC_N1:");  SerialDual.print(NitrousLC1->read(1));
        SerialDual.print("|LC_N2:");  SerialDual.print(NitrousLC2->read(1));
        SerialDual.print("|LC_N3:");  SerialDual.print(NitrousLC3->read(1));
        SerialDual.print("|LC_T:");   SerialDual.print(ThrustLC->read());
        LCLog.last_trigger_ms = time_absolute - (time_absolute % LCLog.interval_ms);
        print_current_poll = true;
      }

      if(time_absolute - ValveLog.last_trigger_ms >= ValveLog.interval_ms) {
        if (!print_current_poll) { SerialDual.print("DATA|"); SerialDual.print((time_elapsed/1000.0), 3); }
        SerialDual.print("|ERV:");    SerialDual.print(EthaneRunValve.state());
        SerialDual.print("|EV:");     SerialDual.print(EthaneVent.state());
        SerialDual.print("|NRV:");    SerialDual.print(NitrousRunValve.state());
        SerialDual.print("|NV:");     SerialDual.print(NitrousVent.state());
        SerialDual.print("|MBV_E:");  SerialDual.print(EthaneMBV->getCurrentDegrees());
        SerialDual.print("|MBV_N:");  SerialDual.print(NitrousMBV->getCurrentDegrees());
        ValveLog.last_trigger_ms = time_absolute - (time_absolute % ValveLog.interval_ms);
        print_current_poll = true;
      }

      if(time_absolute - TCLog.last_trigger_ms >= TCLog.interval_ms) {
        if (!print_current_poll) { SerialDual.print("DATA|"); SerialDual.print((time_elapsed/1000.0), 3); }
        // // SerialDual.print("|TC_C:"); SerialDual.print((ChamberTC->getMilliVolts()/0.041) + ChamberTC->getTemperature());
        // SerialDual.print("|TC_R:"); SerialDual.print(RerouteTC->readHot());
        TCLog.last_trigger_ms = time_absolute - (time_absolute % TCLog.interval_ms);
        print_current_poll = true;
      }

      if (print_current_poll) { SerialDual.println(); }
      break;
    }
  }
}

void processCommand() {
  bool gotCmd = false;
  cmd = "";
  if (serialBuf.feed(Serial, cmd)) {
    gotCmd = true;
  } else if (serial2Buf.feed(Serial2, cmd)) {
    gotCmd = true;
  }

  if (!gotCmd) return;  // nothing ready yet — non-blocking exit

  SerialDual.println(cmd);
  cmd.trim();
  cmd.toUpperCase();

  if (cmd.length() == 0) {
    Serial.println("Empty command");
    return;
  }

  // vent
  if (cmd.equals("")) {
    Serial.println("Empty command");
  }
  else if (cmd.equalsIgnoreCase("ETHANE_VENT_ON")) {
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
  
  // run control -- this is now deprecated
  else if (cmd.equalsIgnoreCase("START")) {
    print_data = true;
    start_time = millis();
  } 
  else if (cmd.equalsIgnoreCase("STOP")) {
    print_data = false;
  }
  else if (cmd.equalsIgnoreCase("STATUS")) {
    status(transmissionFormat);
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
    EMERGENCY_VENT();
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
  else if (cmd.indexOf("ETHANE_TARGET_PRESSURE_") >= 0) {
    String pressure = cmd.substring(cmd.lastIndexOf("_") + 1);
    if (pressure.toFloat()) {
      ETHANE_TARGET_PRESSURE = pressure.toFloat();
    }
  }
  else if (cmd.indexOf("NITROUS_TARGET_PRESSURE_") >= 0) {
    String pressure = cmd.substring(cmd.lastIndexOf("_") + 1);
    if (pressure.toFloat()) {
      NITROUS_TARGET_PRESSURE = pressure.toFloat();
    }
  }
  else if (cmd.equalsIgnoreCase("ETHANE_LC_TARE")) {
    EthaneLC1->tare();
    EthaneLC2->tare();
    EthaneLC3->tare();
  }
  else if (cmd.equalsIgnoreCase("NITROUS_LC_TARE")) {
    NitrousLC1->tare();
    NitrousLC2->tare();
    NitrousLC3->tare();
  }

  else {
    if(full_output) { SerialDual.println("Unknown command.");}
  }
  
}

//===========================EXECUTION============================//
/**
 * MARK: Execution
 * */ 

void setup()
{
  // TURN SERIAL OFF; TURN SERIAL2 ON
  SerialDual.setActive(true, true);
  SerialDual.begin(BAUD_RATE);
  //SerialDual.flush();
  SerialDual.println("START");
  cmd.reserve(32);
  delay(2000);

  EthaneMBV = new MBV(ETHANE_MBV_PIN, ETHANE_ENCODER_PINS);
  NitrousMBV = new MBV(NITROUS_MBV_PIN, NITROUS_ENCODER_PINS);

  SerialDual.println("START1");
  EthaneLC1 = new LoadCell(ETHANE_LC1_PINS);
  EthaneLC2 = new LoadCell(ETHANE_LC2_PINS);
  EthaneLC3 = new LoadCell(ETHANE_LC3_PINS);
  NitrousLC1 = new LoadCell(NITROUS_LC1_PINS);
  NitrousLC2 = new LoadCell(NITROUS_LC2_PINS);
  NitrousLC3 = new LoadCell(NITROUS_LC3_PINS);

  EthaneLC1->setCalFactor(43.27);
  EthaneLC2->setCalFactor(43.92);
  EthaneLC3->setCalFactor(42.46);
  NitrousLC1->setCalFactor(40.0);
  NitrousLC2->setCalFactor(40.0);
  NitrousLC3->setCalFactor(40.0);
  EthaneLC1->join(EthaneLC2, EthaneLC3);
  NitrousLC1->join(NitrousLC2, NitrousLC3);


  SerialDual.println("START2");
  ThrustLC = new ThrustCell();
  ThrustLC->setCalFactor(-5.83);
  SerialDual.println("START3");

  // RerouteTC = new Thermocouple(0x67);
  // RerouteTC2 = new Thermocouple(0x65);
  // RerouteTC3 = new Thermocouple(0x66);

  // pinMode(CHAMBER_TC_PIN, OUTPUT);      // Force SS high to lock Mega in master mode
  // digitalWrite(CHAMBER_TC_PIN, HIGH);
  // ChamberTC = new ADS1118(CHAMBER_TC_PIN);
  delay(100);
  SerialDual.println("START4");

  // ChamberTC->begin();
  // // ChamberTC->setSamplingRate(ChamberTC->RATE_16SPS);
  // // ChamberTC->setInputSelected(ChamberTC->DIFF_0_1);
  // // ChamberTC->setFullScaleRange(ChamberTC->FSR_0256);

  EthaneHeater1 = new Heater(ETHANE_HEATER_1_PIN, &EthaneUpstreamPT);
  EthaneHeater2 = new Heater(ETHANE_HEATER_2_PIN, &EthaneUpstreamPT);
  NitrousHeater1 = new Heater(NITROUS_HEATER_1_PIN, &NitrousUpstreamPT);
  NitrousHeater2 = new Heater(NITROUS_HEATER_2_PIN, &NitrousUpstreamPT);

  EthaneHeater1->setTarget(ETHANE_TARGET_PRESSURE);
  EthaneHeater2->setTarget(ETHANE_TARGET_PRESSURE);
  NitrousHeater1->setTarget(NITROUS_TARGET_PRESSURE);
  NitrousHeater2->setTarget(NITROUS_TARGET_PRESSURE);
  
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
  }  SerialDual.println("END");

}

void loop()
{
  time_absolute = millis();
  time_elapsed = time_absolute - start_time;

  EthaneMBV->update();
  NitrousMBV->update();

  // LOG
  if (print_data)
  {
    status(transmissionFormat);
  }

  if (time_elapsed/1000.0 > runtime && runtime >= 0) {
    print_data = false;
  }

  // COMMANDS
  processCommand();

  // Scheduling Updates
  if(time_absolute - ValveSchedulePoll.last_trigger_ms >= ValveSchedulePoll.interval_ms) {
    EthaneVent.checkScheduledActuation();
    NitrousVent.checkScheduledActuation();
    EthaneRunValve.checkScheduledActuation();
    NitrousRunValve.checkScheduledActuation();
    EthaneMBV->checkScheduledActuation();
    NitrousMBV->checkScheduledActuation();

    ValveSchedulePoll.last_trigger_ms = time_absolute - (time_absolute % ValveSchedulePoll.interval_ms);
  }

  if(heaters_active && time_absolute - HeaterPoll.last_trigger_ms >= HeaterPoll.interval_ms) {
    EthaneHeater1->update();
    EthaneHeater2->update();
    NitrousHeater1->update();
    NitrousHeater2->update();

    HeaterPoll.last_trigger_ms = time_absolute - (time_absolute % HeaterPoll.interval_ms);
  }

  // Redlines
  if(time_absolute - RedlinePoll.last_trigger_ms >= RedlinePoll.interval_ms) {

    if (current_highest_redline > 0) current_highest_redline--;
    for (Redline r : Redlines) {
      int current_highest_redline = std::max(r.checkTrigger(current_highest_redline), current_highest_redline);
    }

    RedlinePoll.last_trigger_ms = time_absolute - (time_absolute % RedlinePoll.interval_ms);
  }

  if(static_fire_initializing) {
    staticFire();
  }
  

}