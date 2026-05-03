#pragma once

#include "main.h"

//============================REDLINES============================//
// Trigger Conditions
bool EthaneOverpressureCondition();
bool NitrousOverpressureCondition();
bool EthaneMBVOpenFailureCondition();
bool NitrousMBVOpenFailureCondition();
bool EthaneMBVCloseFailureCondition();
bool NitrousMBVCloseFailureCondition();
bool CombustionPropogationCondition();
bool InlineThermalDecompCondition();
bool LostLoadCellCondition();
bool EthaneUnderweightCondition();
bool NitrousUnderweightCondition();
bool EthaneOverweightCondition();
bool NitrousOverweightCondition();

// Redline Responses
void EthaneOverpressureResponse();
void NitrousOverpressureResponse();
void EthaneMBVOpenFailureResponse();
void NitrousMBVOpenFailureResponse();
void EthaneMBVCloseFailureResponse();
void NitrousMBVCloseFailureResponse();
void CombustionPropogationResponse();
void InlineThermalDecompResponse();
void LostLoadCellResponse();
void EthaneUnderweightResponse();
void NitrousUnderweightResponse();
void EthaneOverweightResponse();
void NitrousOverweightResponse();

//=========================TEST SEQUENCES=========================//
void neutralizeAll();
void EMERGENCY_VENT();
void ventEthane();
void ventNitrous();
void coldFlowEthane(long duration_ms);
void coldFlowNitrous(long duration_ms);
void staticFire();