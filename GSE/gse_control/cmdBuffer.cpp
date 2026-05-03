#include "cmdBuffer.h"
  // Returns true and populates `out` when a full line is ready

CmdBuffer::CmdBuffer() {
  buf.reserve(CMD_BUF_SIZE);
}
bool CmdBuffer::feed(Stream &s, String &out) {
  while (s.available()) {
    char c = s.read();

    if (c == '\r') continue;

    if (c == '\n') {
      out = buf;
      buf = "";
      overflow = false;
      return true;
    }

    // Write the character, THEN check if we're now full
    buf += c;

    if (head >= CMD_BUF_SIZE - 1) {
      // Buffer full with no newline — discard and reset
      head = 0;
      overflow = true;
    }
  }
  return false;
}