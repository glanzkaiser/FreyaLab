// FileParser.h
#pragma once
#include "Netlist.h"
#include "CircuitComponents.h" // Ensures visibility for shared GateType enum
#include <string>
#include <map>
#include <vector>

struct SubcircuitBlueprint {
    std::string name;
    std::vector<std::string> formalPorts; // E.g., {"IN+", "IN-", "OUT", "VCC", "VEE"}
    std::vector<std::vector<std::string>> internalLines; // Raw tokens of child components
};

class FileParser {
private:
    // Maps netlist string labels (e.g., "GND", "VCC", "net1") to clean 1-based matrix node integers
    std::map<std::string, int> nodeMap;
    int nodeCounter = 0;

    // Tracks registered subcircuit macro structures by uppercase identifiers
    std::map<std::string, SubcircuitBlueprint> subcircuitRegistry;

    int getOrAssignNode(const std::string& nodeName);
    // Modified: Now accepts the active netlist workspace to record string names dynamically
    int getOrAssignNode(const std::string& nodeName, Netlist& netlist);

    double parseValue(const std::string& valStr);
    std::vector<int> parseNodes(const std::vector<std::string>& tokens, size_t count);
    // Modified: Now passes the working netlist reference down to internal assign routines
    std::vector<int> parseNodes(const std::vector<std::string>& tokens, size_t count, Netlist& netlist);

    // Safely maps string tokens into matching internal gate enum options
    GateType parseGateType(const std::string& gateStr);

    // --- SUB-CIRCUIT MACRO INTERNAL EXPANSION IMPLEMENTATIONS ---
    /**
     * Recursively expands subcircuit declarations (X-cards) down to true primitives,
     * rewriting local port text maps securely to absolute parent tracking coordinates.
     */
    void flattenSubcircuitInstance(
        const std::string& instanceName,
        const std::string& subcircuitTypeName,
        const std::vector<std::string>& actualNodes,
        const std::string& hierarchicalPrefix,
        std::vector<std::vector<std::string>>& globalFlattenedTokens
    );
public:
    FileParser();
    
    // Parses the input file line-by-line and loads instantiated components into the target netlist
    void parseNetlistFile(const std::string& filename, Netlist& netlist);
};

