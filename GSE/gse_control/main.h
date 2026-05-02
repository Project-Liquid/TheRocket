#pragma once

#include "transducer.h"
#include "relay.h"
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
#define BAUD_RATE 57600

// Pins
#define ETHANE_UPSTREAM_PIN     A12
#define ETHANE_DOWNSTREAM_PIN   A11
#define NITROUS_UPSTREAM_PIN    A10
#define NITROUS_DOWNSTREAM_PIN  A9
#define REROUTE_PT_PIN          A8
#define CHAMBER_PT_PIN          37  //CHECK--doesn't seem right

#define ETHANE_RUN_PIN   47
#define ETHANE_VENT_PIN  46
#define NITROUS_RUN_PIN  36
#define NITROUS_VENT_PIN 40

#define ETHANE_MBV_PIN  6
#define NITROUS_MBV_PIN 7
#define ETHANE_ENCODER_PINS  44, 42
#define NITROUS_ENCODER_PINS 38, 48

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

#define IGNITOR_PIN 39

// Test Procedure Timings
#define VENT_DELAY              1000
#define VENT_TIME               3000
#define RUN_EQUALIZE_TIME       5000
#define ETHANE_DELAY            200
#define BURNOUT_DELAY           500
#define OVERPRESSURE_VENT_TIME  10000

// Redlines
#define ETHANE_PRESSURE_REDLINE   1100
#define NITROUS_PRESSURE_REDLINE  1100
#define OVERPRESSURE_COUNTS_THRESHOLD 25

#define UNDERWEIGHT_COUNTS_THRESHOLD  20

// PT Calibration Data (PSIG)
#define P_MIN         0.0
#define P_MAX_ETHANE  1000
#define P_MAX_NITROUS 1500
#define P_MAX_CHAMBER ?


//========================GLOBAL VARIABLES========================//
// Run Control
extern bool print_data; // Will always be true
extern unsigned long start_time; 

/**
 * Maximum miliseconds of runtime.
 * Note: Does not stop if negative
 */
extern const int runtime;

/**
 * Gives diagnostic with full details if true, otherwise gives simplified output for logging.
 * Note: Full output is not recommended for long runs due to performance constraints
 * of Serial printing every loop, which can cause lag and missed readings. Use
 * full output for debugging and simplified output for extended
 */
//const bool full_output = false;

enum TransmissionType {
  DIAGNOSTIC,
  COMPRESSED,
  RAW
};

extern TransmissionType transmissionFormat;
#define full_output (transmissionFormat == DIAGNOSTIC)

extern unsigned long time_absolute;
/**
 * Current Groundtime
 * - Defined as milis() - start_time */
extern unsigned long time_elapsed;
extern bool static_fire_initializing;
extern unsigned long static_fire_duration_ms;

extern bool heaters_active;
extern float ETHANE_TARGET_PRESSURE;
extern float NITROUS_TARGET_PRESSURE;

/**
 * The timing of each thing in data log
 * - Has the interval at which to log, and the last time it was logged,
 * so we can check if it's time to log again
 */
struct PollInterval {
  const unsigned int interval_ms;
  unsigned long last_trigger_ms;
};

extern PollInterval PTLog;
extern PollInterval LCLog;
extern PollInterval ValveLog;
extern PollInterval MBVLog;
extern PollInterval TCLog;
extern PollInterval RedlinePoll;
extern PollInterval ValveSchedulePoll;
extern PollInterval HeaterPoll;


extern float ETHANE_WEIGHT_REDLINE;
extern float NITROUS_WEIGHT_REDLINE;

//======================OBJECT DECLARATIONS=======================//

// Pressure Transducers
Transducer EthaneUpstreamPT(ETHANE_UPSTREAM_PIN, P_MIN, P_MAX_ETHANE);
Transducer EthaneDownstreamPT(ETHANE_DOWNSTREAM_PIN, P_MIN, P_MAX_ETHANE);
Transducer NitrousUpstreamPT(NITROUS_UPSTREAM_PIN, P_MIN, P_MAX_NITROUS);
Transducer NitrousDownstreamPT(NITROUS_DOWNSTREAM_PIN, P_MIN, P_MAX_NITROUS);
Transducer ReroutePT(REROUTE_PT_PIN, P_MIN, P_MAX_NITROUS);
//Transducer ChamberPT(CHAMBER_PT_PIN, P_MIN, P_MAX_CHAMBER);

// Solenoids
Relay EthaneRunValve(ETHANE_RUN_PIN);
Relay EthaneVent(ETHANE_VENT_PIN);
Relay NitrousRunValve(NITROUS_RUN_PIN);
Relay NitrousVent(NITROUS_VENT_PIN);

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
ADS1118* ChamberTC; 

// Ignitor
Relay Ignitor(IGNITOR_PIN);

/** Serial Communication for Radios
 * - Serial is hardline
 * - Serial2 is radio
*/
SerialDualClass SerialDual(Serial, Serial2);
