#include "freyalab/LogicGates.h"
#include <string>

uint8_t LogicGates::evaluateGate(GateType gate, uint8_t a, uint8_t b) {
    // Keep operations strictly isolated within the 4-bit nibble space bounds (0x0F)
    switch (gate) {
        case GateType::NOT:
            return (~a) & 0x0F;
        case GateType::AND:
            return (a & b) & 0x0F;
        case GateType::OR:
            return (a | b) & 0x0F;
        case GateType::XOR:
            return (a ^ b) & 0x0F;
        case GateType::NAND:
            return (~(a & b)) & 0x0F;
        case GateType::NOR:
            return (~(a | b)) & 0x0F;
        case GateType::XNOR:
            return (~(a ^ b)) & 0x0F;
        default:
            return 0x00;
    }
}

std::string LogicGates::toBinary4Bit(uint8_t value) {
    // Optimized fixed-size allocation lookup bypass to avoid unneeded loop iteration overhead
    std::string s = "0000";
    s[0] = ((value & 0x08) ? '1' : '0');
    s[1] = ((value & 0x04) ? '1' : '0');
    s[2] = ((value & 0x02) ? '1' : '0');
    s[3] = ((value & 0x01) ? '1' : '0');
    return s;
}
