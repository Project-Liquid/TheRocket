#pragma once

#include "transducer.h"
#include "relay.h"
#include "mbv.h"
#include "loadCell.h"
#include "thermocouple.h"
#include "heater.h"
#include "thrustCell.h"
#include "serial.h"
#include "redline.h"
#include "procedures.h"

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
#define CHAMBER_PT_PIN          A13  

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

#define ETHANE_OVERPRESSURE_PRIORITY        1
#define NITROUS_OVERPRESSURE_PRIORITY       1
#define ETHANE_MBV_OPEN_FAILURE_PRIORITY    1
#define NITROUS_MBV_OPEN_FAILURE_PRIORITY   1
#define ETHANE_MBV_CLOSE_FAILURE_PRIORITY   1
#define NITROUS_MBV_CLOSE_FAILURE_PRIORITY  1
#define COMBUSTION_PROPOGATION_PRIORITY     1
#define INLINE_THERMAL_DECOMP_PRIORITY      1
#define LOST_LOAD_CELL_PRIORITY             1
#define ETHANE_UNDERWEIGHT_PRIORITY         1
#define NITROUS_UNDERWEIGHT_PRIORITY        1
#define ETHANE_OVERWEIGHT_PRIORITY          1
#define NITROUS_OVERWEIGHT_PRIORITY         1

#define OVERPRESSURE_COUNTS_THRESHOLD               20
#define MBV_FAILURE_COUNTS_THRESHOLD                20
#define COMBUSTION_PROPOGATION_COUNTS_THRESHOLD     20
#define INLINE_THERMAL_DECOMP_COUNTS_THRESHOLD      20
#define LOST_LOAD_CELL_COUNTS_THRESHOLD             20
#define UNDERWEIGHT_COUNTS_THRESHOLD                20
#define OVERWEIGHT_COUNTS_THRESHOLD                 20

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
extern Transducer EthaneUpstreamPT;
extern Transducer EthaneDownstreamPT;
extern Transducer NitrousUpstreamPT;
extern Transducer NitrousDownstreamPT;
extern Transducer ReroutePT;
//Transducer ChamberPT;

// Solenoids
extern Relay EthaneRunValve;
extern Relay EthaneVent;
extern Relay NitrousRunValve;
extern Relay NitrousVent;

// Motorized Ball Valves
extern MBV* EthaneMBV;
extern MBV* NitrousMBV;

// Load Cells
extern LoadCell* EthaneLC2;
extern LoadCell* EthaneLC1;
extern LoadCell* EthaneLC3;
extern LoadCell* NitrousLC1;
extern LoadCell* NitrousLC2;
extern LoadCell* NitrousLC3;

extern ThrustCell* ThrustLC;

// Tank Heaters
extern Heater* EthaneHeater1;
extern Heater* EthaneHeater2;
extern Heater* NitrousHeater1;
extern Heater* NitrousHeater2;

// Thermocouple
extern Thermocouple* RerouteTC;
extern Thermocouple* RerouteTC2;
extern Thermocouple* RerouteTC3;
extern ADS1118* ChamberTC; 

// Ignitor
extern Relay Ignitor;

/** Serial Communication for Radios
 * - Serial is hardline
 * - Serial2 is radio
*/
extern SerialDualClass SerialDual;

// Redlines
extern int current_highest_redline;
extern Redline EthaneOverpressure;
extern Redline NitrousOverpressure;
extern Redline EthaneMBVOpenFailure;
extern Redline NitrousMBVOpenFailure;
extern Redline EthaneMBVCloseFailure;
extern Redline NitrousMBVCloseFailure;
extern Redline CombustionPropogation;
extern Redline InlineThermalDecomp;
extern Redline LostLoadCell;
extern Redline EthaneUnderweight;
extern Redline NitrousUnderweight;
extern Redline EthaneOverweight;
extern Redline NitrousOverweight;

