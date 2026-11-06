// FileParser.h
#pragma once
#include "Netlist.h"
#include "CircuitComponents.h" // Ensures visibility for shared GateType enum
#include <string>
#include <map>
#include <vector>

class FileParser {
private:
    // Maps netlist string labels (e.g., "GND", "VCC", "net1") to clean 1-based matrix node integers
    std::map<std::string, int> nodeMap;
    int nodeCounter = 0;

    int getOrAssignNode(const std::string& nodeName);
    double parseValue(const std::string& valStr);
    std::vector<int> parseNodes(const std::vector<std::string>& tokens, size_t count);
    
    // Safely maps string tokens into matching internal gate enum options
    GateType parseGateType(const std::string& gateStr);

public:
    FileParser();
    
    // Parses the input file line-by-line and loads instantiated components into the target netlist
    void parseNetlistFile(const std::string& filename, Netlist& netlist);
};

