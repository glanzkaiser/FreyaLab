// LogicGates.h
#pragma once
#include "CircuitComponents.h" // Natively leverages the shared GateType enum
#include <string>
#include <cstdint>

class LogicGates {
public:
    // Delete constructor to prevent useless structural heap/stack instances
    LogicGates() = delete;

    // Core generic single-gate evaluator dispatched by your standard SPICE GateType enum
    static uint8_t evaluateGate(GateType gate, uint8_t a, uint8_t b = 0);

    // Formats a lower 4-bit nibble into a binary string for precise column printing
    static std::string toBinary4Bit(uint8_t value);
};

