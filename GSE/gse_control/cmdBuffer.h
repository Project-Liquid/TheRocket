#pragma once
#include <Arduino.h>
#define CMD_BUF_SIZE 32

class CmdBuffer {
private:
  String buf;
  int head = 0;
  bool overflow = false;

public:
  CmdBuffer();
  bool feed(Stream &s, String &out);
};