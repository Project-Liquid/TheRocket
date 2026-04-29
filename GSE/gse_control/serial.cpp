#include "serial.h"

SerialDualClass::SerialDualClass(HardwareSerial& serialA,
                HardwareSerial& serialB)
    : s1(serialA), s2(serialB) {}

void SerialDualClass::begin(unsigned long baud) {
    s1.begin(baud);
    s2.begin(baud);
}

// Required override from Print
size_t SerialDualClass::write(uint8_t c) {
    s1.write(c);
    s2.write(c);
    return 1;
}

// Optional but more efficient
size_t SerialDualClass::write(const uint8_t* buffer, size_t size) {
    s1.write(buffer, size);
    s2.write(buffer, size);
    return size;
}
