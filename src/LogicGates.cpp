#include "freyalab/LogicGates.h"
#include <bitset>

// --- BUS MULTI-BIT ---
uint8_t LogicGates::evaluateBusGate(GateType gate, uint8_t a, uint8_t b) {
    uint8_t maskA = a & 0x0F;
    uint8_t maskB = b & 0x0F;

    switch (gate) {
        case GateType::NOT:  return (~maskA) & 0x0F;
        case GateType::AND:  return (maskA & maskB) & 0x0F;
        case GateType::OR:   return (maskA | maskB) & 0x0F;
        case GateType::XOR:  return (maskA ^ maskB) & 0x0F;
        case GateType::NAND: return (~(maskA & maskB)) & 0x0F;
        case GateType::NOR:  return (~(maskA | maskB)) & 0x0F;
        case GateType::XNOR: return (~(maskA ^ maskB)) & 0x0F;
        default:             return 0x00;
    }
}

// --- SPICE MNA ---
uint8_t LogicGates::evaluateScalarGate(GateType gate, double vA, double vB, double vThresh) {
    uint8_t inA = (vA >= vThresh) ? 1 : 0;
    uint8_t inB = (vB >= vThresh) ? 1 : 0;

    switch (gate) {
        case GateType::NOT:  return !inA;
        case GateType::AND:  return inA & inB;
        case GateType::OR:   return inA | inB;
        case GateType::XOR:  return inA ^ inB;
        case GateType::NAND: return !(inA & inB);
        case GateType::NOR:  return !(inA | inB);
        case GateType::XNOR: return !(inA ^ inB);
        default:             return 0;
    }
}

std::string LogicGates::toBinary4Bit(uint8_t value) {
    return std::bitset<4>(value & 0x0F).to_string();
}
