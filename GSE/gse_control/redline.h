#pragma once

class Redline {
private:
  bool (*trigger_condition)();
  void (*response)();
  int priority;
  int counts;
  int counts_threshold;

public:
  Redline::Redline(bool (*trigger_condition)(), void (*response)(), int priority, int counts_threshold);
  void Redline::setTriggerCondition(bool (*trigger_condition)());
  void Redline::setResponse(void (*response)());
  void Redline::setThreshold(int counts_threshold);
  int Redline::checkTrigger(int current_priority);
};