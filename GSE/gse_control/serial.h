#pragma once

#include <Arduino.h>

class SerialDualClass : public Print {
private:
    HardwareSerial& s1;
    HardwareSerial& s2;
    bool s1_active = true;
    bool s2_active = true;

public:
    SerialDualClass(HardwareSerial& serialA, HardwareSerial& serialB);

    void begin(unsigned long baud);

    void setActive(bool s1_active, bool s2_active);

    // Required override from Print
    virtual size_t write(uint8_t c) override;

    // Optional but more efficient
    virtual size_t write(const uint8_t* buffer, size_t size) override;

    using Print::write;
};