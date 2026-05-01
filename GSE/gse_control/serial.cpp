#include "serial.h"

SerialDualClass::SerialDualClass(HardwareSerial& serialA,
                HardwareSerial& serialB)
    : s1(serialA), s2(serialB) {}

void SerialDualClass::setActive(bool s1_active, bool s2_active) {
    this->s1_active = s1_active;
    this->s2_active = s2_active;
}

void SerialDualClass::begin(unsigned long baud) {
    if (s1_active) s1.begin(baud);
    if (s2_active) s2.begin(baud);
}

// Required override from Print
size_t SerialDualClass::write(uint8_t c) {
    if (s1_active) s1.write(c);
    if (s2_active) s2.write(c);
    return 1;
}

// Optional but more efficient
size_t SerialDualClass::write(const uint8_t* buffer, size_t size) {
    s1.write(buffer, size);
    s2.write(buffer, size);
    return size;
}
