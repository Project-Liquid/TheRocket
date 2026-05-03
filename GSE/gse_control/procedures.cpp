#include "procedures.h"

//============================REDLINES============================//

// Trigger Conditions
bool EthaneOverpressureCondition() { return EthaneUpstreamPT.readPressure() > ETHANE_PRESSURE_REDLINE; }
bool NitrousOverpressureCondition() { return NitrousUpstreamPT.readPressure() > NITROUS_PRESSURE_REDLINE; }
bool EthaneMBVCloseFailureCondition() { return; }
bool NitrousMBVCloseFailureCondition() { return; }
bool CombustionPropogationCondition() { return; }
bool InlineThermalDecompCondition() { return; }
//bool LostLoadCellCondition() { return; }
bool EthaneUnderweightCondition() { return EthaneLC1->readJoint(1) < ETHANE_WEIGHT_REDLINE; }
//bool NitrousUnderweightCondition() { return;}
//bool EthaneOverweightCondition() { return; }
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
void NitrousMBVCloseFailureResponse() {}
void CombustionPropogationResponse() {
  
}
void InlineThermalDecompResponse() {}
//void LostLoadCellResponse() {}
void EthaneUnderweightResponse() {EthaneRunValve.neutralize();}
//void NitrousUnderweightResponse() {}
//void EthaneOverweightResponse() {}
void NitrousOverweightResponse() {NitrousRunValve.neutralize();}

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
  //NitrousRunValve.open();
  EthaneRunValve.open();

  Ignitor.setNextActuation(5000, true);
  
  if (ChamberTC->getTemperature() > 100) {
    NitrousMBV->next_90();
    EthaneMBV->setNextActuation(ETHANE_DELAY);
    NitrousMBV->setNextActuation(ETHANE_DELAY + static_fire_duration_ms);
    EthaneMBV->setNextActuation(ETHANE_DELAY + static_fire_duration_ms + BURNOUT_DELAY);
    NitrousRunValve.setNextActuation(ETHANE_DELAY + static_fire_duration_ms + BURNOUT_DELAY + 100, false);
    EthaneRunValve.setNextActuation(ETHANE_DELAY + static_fire_duration_ms + BURNOUT_DELAY + 100, false);
    
    // Vent
    EthaneVent.setNextActuation(ETHANE_DELAY + static_fire_duration_ms + BURNOUT_DELAY + 100 + VENT_DELAY, true);
    EthaneVent.setNextActuation(ETHANE_DELAY + static_fire_duration_ms + BURNOUT_DELAY + 100 + VENT_DELAY + VENT_TIME, false);
    NitrousVent.setNextActuation(ETHANE_DELAY + static_fire_duration_ms + BURNOUT_DELAY + VENT_DELAY + VENT_TIME + 1100, true);
    NitrousVent.setNextActuation(ETHANE_DELAY + static_fire_duration_ms + BURNOUT_DELAY + VENT_DELAY + 2*VENT_TIME + 1100, false);

    static_fire_initializing = false;
  }
}
