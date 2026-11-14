// FileParser.h
#pragma once
#include "Netlist.h"
#include "CircuitComponents.h" // Visibility for NetNode, MnaIndex, and GateType
#include "BSIM.h"              // Visibility for BsimCmgParameters and MultiGateGeometry
#include <map>
#include <string>
#include <vector>

struct SubcircuitBlueprint {
    std::string name;
    std::vector<std::string> formalPorts; // E.g., {"IN+", "IN-", "OUT", "VCC", "VEE"}
    std::vector<std::vector<std::string>> internalLines; // Raw tokens of child components
};

class FileParser {
private:
    // Maps uppercase netlist string labels (e.g., "GND", "VCC", "NET1") to type-safe 1-based NetNode entities
    std::map<std::string, NetNode> nodeMap;
    int nodeCounter = 0;

    // Tracks registered subcircuit macro structures by uppercase identifiers
    std::map<std::string, SubcircuitBlueprint> subcircuitRegistry;

    NetNode getOrAssignNode(const std::string& nodeName);
    // Records string node names dynamically inside the active netlist workspace
    NetNode getOrAssignNode(const std::string& nodeName, Netlist& netlist);

    double parseValue(const std::string& valStr);
    std::vector<NetNode> parseNodes(const std::vector<std::string>& tokens, size_t count);
    // Passes the working netlist reference down to internal node assignment routines
    std::vector<NetNode> parseNodes(const std::vector<std::string>& tokens, size_t count, Netlist& netlist);

    // Safely maps string tokens into matching internal gate enum options
    GateType parseGateType(const std::string& gateStr);

    // Parses inline key=value BSIM-CMG 112.1.0 parameters (L, NFIN, TFIN, HFIN, VTH0, COREMOD, etc.)
    void parseBsimInlineParams(const std::vector<std::string>& tokens,
                               size_t startIdx,
                               BsimCmgParameters& params);

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

