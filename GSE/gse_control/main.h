#pragma once

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

//===========================CONSTANTS============================//
// Pins
#define ETHANE_UPSTREAM_PIN     A12
#define ETHANE_DOWNSTREAM_PIN   A11
#define NITROUS_UPSTREAM_PIN    A10
#define NITROUS_DOWNSTREAM_PIN  A9
#define REROUTE_PT_PIN          A8

#define ETHANE_RUN_PIN   47
#define ETHANE_VENT_PIN  46
#define NITROUS_RUN_PIN  36
#define NITROUS_VENT_PIN 40

#define ETHANE_MBV_PIN  6
#define NITROUS_MBV_PIN 7
#define ETHANE_ENCODER_PIN_1  44
#define ETHANE_ENCODER_PIN_2  42
#define NITROUS_ENCODER_PIN_1 38
#define NITROUS_ENCODER_PIN_2 48

#define ETHANE_LC1_PINS   22, 23
#define ETHANE_LC2_PINS   24, 25
#define ETHANE_LC3_PINS   26, 27
#define NITROUS_LC1_PINS  34, 35
#define NITROUS_LC2_PINS  32, 33
#define NITROUS_LC3_PINS  30, 31

#define CHAMBER_TC_PIN 53

#define ETHANE_HEATER_1_PIN   49
#define ETHANE_HEATER_2_PIN   41
#define NITROUS_HEATER_1_PIN  43
#define NITROUS_HEATER_2_PIN  45

// Test Procedure Timings
#define VENT_DELAY        1000
#define VENT_TIME         3000
#define RUN_EQUALIZE_TIME 5000
#define ETHANE_DELAY      200
#define BURNOUT_DELAY     500

// Redlines
#define ETHANE_PRESSURE_REDLINE   1100
#define NITROUS_PRESSURE_REDLINE  1100
#define OVERPRESSURE_COUNTS_THRESHOLD 25

// PT Calibration Data (PSIG)
#define P_MIN         0.0
#define P_MAX_ETHANE  1000
#define P_MAX_NITROUS 1500


//========================GLOBAL VARIABLES========================//
// Run Control
bool print_data = true; // Will always be true
long start_time = 0; //1410065408; // max value placeholder

/**
 * Maximum miliseconds of runtime.
 * Note: Does not stop if negative
 */
const int runtime = -1;

/**
 * Gives diagnostic with full details if true, otherwise gives simplified output for logging.
 * Note: Full output is not recommended for long runs due to performance constraints
 * of Serial printing every loop, which can cause lag and missed readings. Use
 * full output for debugging and simplified output for extended
 */
const bool full_output = false;

/**
 * Current Groundtime
 * - Defined as milis() - start_time */
long elapsed = 0;
bool static_fire_initializing = false;
long static_fire_duration_ms = 0;

bool heaters_active = false;
float ETHANE_TARGET_PRESSURE = 0;
float NITROUS_TARGET_PRESSURE = 0;

/**
 * The timing of each thing in data log
 * - Has the interval at which to log, and the last time it was logged,
 * so we can check if it's time to log again
 */
struct LogInterval {
  const unsigned int interval_ms;
  unsigned long last_trigger_ms;
};

LogInterval PTLog{50, 0};
LogInterval LCLog{200, 0};
LogInterval ValveLog{200, 0};
LogInterval MBVLog{200, 0};
LogInterval TCLog{200, 0};


//======================OBJECT DECLARATIONS=======================//

// Pressure Transducers
Transducer EthaneUpstreamPT(ETHANE_UPSTREAM_PIN, P_MIN, P_MAX_ETHANE);
Transducer EthaneDownstreamPT(ETHANE_DOWNSTREAM_PIN, P_MIN, P_MAX_ETHANE);
Transducer NitrousUpstreamPT(NITROUS_UPSTREAM_PIN, P_MIN, P_MAX_NITROUS);
Transducer NitrousDownstreamPT(NITROUS_DOWNSTREAM_PIN, P_MIN, P_MAX_NITROUS);
Transducer ReroutePT(REROUTE_PT_PIN, P_MIN, P_MAX_NITROUS);

// Solenoids
Solenoid EthaneRunValve(ETHANE_RUN_PIN);
Solenoid EthaneVent(ETHANE_VENT_PIN);
Solenoid NitrousRunValve(NITROUS_RUN_PIN);
Solenoid NitrousVent(NITROUS_VENT_PIN);

// Motorized Ball Valves
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
Heater* EthaneHeater1;
Heater* EthaneHeater2;
Heater* NitrousHeater1;
Heater* NitrousHeater2;

// Thermocouple
Thermocouple* RerouteTC;
Thermocouple* RerouteTC2;
Thermocouple* RerouteTC3;
/** TODO: This doesn't really work. */
ADS1118* ChamberTC; 

/** Serial Communication for Radios
 * - Serial is hardline
 * - Serial2 is radio
*/
SerialDualClass SerialDual(Serial, Serial2);
