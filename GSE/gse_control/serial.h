#pragma once

#include <Arduino.h>

class SerialDualClass : public Print {
private:
    HardwareSerial& s1;
    HardwareSerial& s2;

public:
    SerialDualClass(HardwareSerial& serialA, HardwareSerial& serialB);

    void begin(unsigned long baud);

    // Required override from Print
    virtual size_t write(uint8_t c) override;

    // Optional but more efficient
    virtual size_t write(const uint8_t* buffer, size_t size) override;

    using Print::write;
};