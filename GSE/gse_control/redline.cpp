#include "redline.h"

Redline::Redline(bool (*trigger_condition)(), void (*response)(), int priority, int counts_threshold) 
  : trigger_condition(trigger_condition), response(response), priority(priority), counts_threshold(counts_threshold) {
}

void Redline::setTriggerCondition(bool (*trigger_condition)()) {
  this->trigger_condition = trigger_condition;
}

void Redline::setResponse(void (*response)()) {
  this->response = response;
}

void Redline::setThreshold(int counts_threshold) {
  this->counts_threshold = counts_threshold;
}

int Redline::checkTrigger(int current_priority) {
  if (trigger_condition()) {
    counts++;
  } else if (counts > 0) {
    counts--;
  }

  if (counts > counts_threshold && priority > current_priority) {
    counts = 0;
    return priority;
  }
  return -1;
}
