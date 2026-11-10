#pragma once
#include "freyalab/CircuitComponents.h"
#include <string>
#include <cstdint>

class LogicGates {
public:
    LogicGates() = delete;


    static uint8_t evaluateBusGate(GateType gate, uint8_t a, uint8_t b = 0);

    static uint8_t evaluateScalarGate(GateType gate, double vA, double vB = 0.0, double vThresh = 2.5);

    static std::string toBinary4Bit(uint8_t value);
};

