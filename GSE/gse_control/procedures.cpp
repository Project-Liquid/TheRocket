#include "procedures.h"

//============================REDLINES============================//

// Trigger Conditions
bool EthaneOverpressureCondition() { return EthaneUpstreamPT.readPressure() > ETHANE_PRESSURE_REDLINE; }
bool NitrousOverpressureCondition() { return NitrousUpstreamPT.readPressure() > NITROUS_PRESSURE_REDLINE; }
bool EthaneMBVCloseFailureCondition() { return; }  //Unessecary, controlled by cold flow/ static fire
bool NitrousMBVCloseFailureCondition() { return; } //Unessecary, controlled by cold flow/ static fire
bool CombustionPropogationCondition() { 
  return static_fire_steady && ((fabs((ReroutePT.readPressure()-ChamberPT.readPressure()))/ReroutePT.readPressure() >0.1)||
   ((fabs(EthaneDownstreamPT.readPressure()-ChamberPT.readPressure()))/EthaneDownstreamPT.readPressure() >0.1)||
    ((0.95*ReroutePT.readPressure() < ChamberPT.readPressure())) ||
    ((0.95*EthaneDownstreamPT.readPressure() < ChamberPT.readPressure())));}
bool InlineThermalDecompCondition() { return RerouteTC->readHot() > REROUTE_TC_REDLINE;}
bool LostLoadCellCondition() { return; } //Unnecessary, operator control
bool EthaneUnderweightCondition() { return EthaneLC1->readJoint(1) < ETHANE_WEIGHT_REDLINE; }
bool NitrousUnderweightCondition() { return;} //Unnecessary, operator control
bool EthaneOverweightCondition() { return; } //Unnecessary, operator control
bool NitrousOverweightCondition() { return NitrousLC1->readJoint(1) > NITROUS_WEIGHT_REDLINE;}

// Redline Responses
void EthaneOverpressureResponse() {
  // Clear schedules to prevent conflicts
  EthaneRunValve.clearSchedule();
  EthaneVent.clearSchedule();
  
  // Solenoid
  EthaneRunValve.open();

  // VENT
  EthaneVent.setNextActuation(0, true);
  EthaneVent.setNextActuation(VENT_TIME, false);
}

void NitrousOverpressureResponse() {
  // Clear schedules to prevent conflicts
  NitrousRunValve.clearSchedule();
  NitrousVent.clearSchedule();

  // Solenoid
  NitrousRunValve.open();

  // VENT
  NitrousVent.setNextActuation(0, true);
  NitrousVent.setNextActuation(VENT_TIME, false);  
}

void EthaneMBVCloseFailureResponse() {

}
void NitrousMBVCloseFailureResponse() {

}
void CombustionPropogationResponse() {
  neutralizeAll();

  NitrousMBV->next_90();
  NitrousRunValve.close();
  EthaneMBV->next_90();
  EthaneRunValve.close();

  NitrousVent.open();
  NitrousVent.setNextActuation(3000, false);
  EthaneVent.open();
  EthaneVent.setNextActuation(3000, false);
}
void InlineThermalDecompResponse() {
  neutralizeAll();

  NitrousMBV->next_90();
  NitrousRunValve.close();
  EthaneMBV->setNextActuation(1000);
  EthaneRunValve.setNextActuation(1000, false);
}
void LostLoadCellResponse() {}
void EthaneUnderweightResponse() {
  EthaneRunValve.neutralize();
  EthaneRunValve.close();
  }
void NitrousUnderweightResponse() {}
void EthaneOverweightResponse() {}
void NitrousOverweightResponse() {
  NitrousRunValve.neutralize();
  NitrousRunValve.close();
  }

//=========================TEST SEQUENCES=========================//
void neutralizeAll() {
  EthaneVent.neutralize();
  NitrousVent.neutralize();
  EthaneRunValve.neutralize();
  NitrousRunValve.neutralize();
  EthaneMBV->neutralize();
  NitrousMBV->neutralize();
  EthaneHeater1->off();
  EthaneHeater2->off();
  NitrousHeater1->off();
  NitrousHeater2->off();
}

void EMERGENCY_VENT() {
  // clear the schedules
  neutralizeAll();

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
  EthaneVent.setNextActuation(VENT_DELAY, true);
  EthaneVent.setNextActuation(VENT_DELAY + VENT_TIME, false);
}

void ventNitrous() {
  NitrousRunValve.close();
  NitrousVent.setNextActuation(VENT_DELAY, true);
  NitrousVent.setNextActuation(VENT_DELAY + VENT_TIME, false);
}

void coldFlowEthane(long duration_ms) {
  EthaneRunValve.clearSchedule();
  EthaneMBV->clearSchedule();

  EthaneRunValve.open();

  EthaneMBV->setNextActuation(RUN_EQUALIZE_TIME);
  EthaneMBV->setNextActuation(RUN_EQUALIZE_TIME + duration_ms);

  EthaneRunValve.setNextActuation(RUN_EQUALIZE_TIME + duration_ms, false);
}

void coldFlowNitrous(long duration_ms) {
  NitrousRunValve.open();

  NitrousMBV->setNextActuation(RUN_EQUALIZE_TIME);
  NitrousMBV->setNextActuation(RUN_EQUALIZE_TIME + duration_ms);

  NitrousRunValve.setNextActuation(RUN_EQUALIZE_TIME + duration_ms, false);
}

void staticFire() {
  static_fire_start_ms = millis();
  //NitrousRunValve.open();
  EthaneRunValve.open();

  Ignitor.setNextActuation(5000, true);

  // if (ChamberTC->readHot() > 100) {
    NitrousMBV->next_90();
    EthaneMBV->setNextActuation(ETHANE_DELAY);
    NitrousMBV->setNextActuation(ETHANE_DELAY + static_fire_duration_ms);
    EthaneMBV->setNextActuation(ETHANE_DELAY + static_fire_duration_ms + BURNOUT_DELAY);
    NitrousRunValve.setNextActuation(ETHANE_DELAY + static_fire_duration_ms, false);
    EthaneRunValve.setNextActuation(ETHANE_DELAY + static_fire_duration_ms + BURNOUT_DELAY, false);
    
    // Vent
    EthaneVent.setNextActuation(ETHANE_DELAY + static_fire_duration_ms + BURNOUT_DELAY + 100 + VENT_DELAY, true);
    EthaneVent.setNextActuation(ETHANE_DELAY + static_fire_duration_ms + BURNOUT_DELAY + 100 + VENT_DELAY + VENT_TIME, false);
    NitrousVent.setNextActuation(ETHANE_DELAY + static_fire_duration_ms + BURNOUT_DELAY + VENT_DELAY + VENT_TIME + 1100, true);
    NitrousVent.setNextActuation(ETHANE_DELAY + static_fire_duration_ms + BURNOUT_DELAY + VENT_DELAY + 2*VENT_TIME + 1100, false);

    static_fire_initializing = false;
    static_fire_steady = false;
  // }
}
