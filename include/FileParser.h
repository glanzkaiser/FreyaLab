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

struct ModelCardBlueprint {
    std::string modelName;   // Example: "NFET_01V8", "NFIN_CMG"
    std::string modelType;   // Example: "NMOS", "PMOS", "BSIMCMG", "NPN", "PNP", "D"
    std::map<std::string, double> parameters; // Uppercase key=value parameter map
};

class FileParser {
private:
    // Maps uppercase netlist string labels (e.g., "GND", "VCC", "NET1") to type-safe 1-based NetNode entities
    std::map<std::string, NetNode> nodeMap;
    int nodeCounter = 0;

    // Stores global variables from .PARAM directives (e.g., VDD -> 1.8, WMIN -> 0.42e-6)
    std::map<std::string, double> paramRegistry;

    // Stores .MODEL card definitions (e.g., "NFET_01V8" -> ModelCardBlueprint)
    std::map<std::string, ModelCardBlueprint> modelRegistry;

    // Parses .PARAM and .MODEL directives
    void parseParamDirective(const std::vector<std::string>& tokens);
    void parseModelDirective(const std::vector<std::string>& tokens);

    // Evaluates numeric literals OR braced algebraic expressions such as '{2*WMIN + 0.1u}'
    double evaluateExpression(const std::string& exprStr) const;

    // Applies inherited .MODEL parameters prior to overriding with inline instance parameters
    void applyModelCardToBsim(const std::string& modelName, BsimCmgParameters& params) const;

        // Recursively expands .INCLUDE and .LIB <file> <corner> directives while merging '+' continuations
    void collectExpandedLines(const std::string& filepath,
                              const std::string& targetCorner,
                              std::vector<std::string>& outLogicalLines,
                              int depth = 0);

    // Tracks registered subcircuit macro structures by uppercase identifiers
    std::map<std::string, SubcircuitBlueprint> subcircuitRegistry;

    NetNode getOrAssignNode(const std::string& nodeName);
    // Records string node names dynamically inside the active netlist workspace
    NetNode getOrAssignNode(const std::string& nodeName, Netlist& netlist);

    double parseValue(const std::string& valStr) const;
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

