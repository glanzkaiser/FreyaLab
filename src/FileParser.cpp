#include "freyalab/FileParser.h"
#include "freyalab/CircuitComponents.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <cstring>

FileParser::FileParser() {
    nodeMap["0"] = 0;
    nodeMap["GND"] = 0;
    nodeMap["gnd"] = 0;
}

int FileParser::getOrAssignNode(const std::string& nodeName) {
    auto it = nodeMap.find(nodeName);
    if (it != nodeMap.end()) {
        return it->second;
    }
    nodeCounter++;
    nodeMap[nodeName] = nodeCounter;
    return nodeCounter;
}

int FileParser::getOrAssignNode(const std::string& nodeName, Netlist& netlist) {
    std::string upperName = nodeName;
    std::transform(upperName.begin(), upperName.end(), upperName.begin(), ::toupper);

    if (upperName == "0" || upperName == "GND") {
        return 0;
    }

    auto it = nodeMap.find(upperName);
    if (it != nodeMap.end()) {
        return it->second;
    }

    nodeCounter++;
    nodeMap[upperName] = nodeCounter;
    
    // ADD THIS LINE: Pass the clean string label to the netlist workspace
    netlist.setNodeName(nodeCounter, upperName); 
    
    return nodeCounter;
}


double FileParser::parseValue(const std::string& valStr) {
    if (valStr.empty()) return 0.0;
    
    std::string cleanStr = valStr;
    if (cleanStr.back() == ',' || cleanStr.back() == ';') {
        cleanStr.pop_back();
    }

    size_t idx = 0;
    double baseVal = 0.0;
    try {
        baseVal = std::stod(cleanStr, &idx);
    } catch (const std::invalid_argument&) {
        return 0.0;
    } catch (const std::out_of_range&) {
        std::cerr << "Parser Warning: Value out of range [" << cleanStr << "], clamping to 0.0" << std::endl;
        return 0.0;
    }

    if (idx >= cleanStr.size()) return baseVal;

    std::string suffixStr = cleanStr.substr(idx);
    
    // Check for multibyte UTF-8 micro prefixes before running character transforms
    if (suffixStr.find("μ") == 0 || suffixStr.find("µ") == 0) {
        return baseVal * 1e-6;
    }

    std::transform(suffixStr.begin(), suffixStr.end(), suffixStr.begin(), [](unsigned char c){ return std::tolower(c); });

    if (suffixStr.empty()) return baseVal;

    char suffix = suffixStr[0];
    if (suffix == 'f') return baseVal * 1e-15;
    if (suffix == 'p') return baseVal * 1e-12;
    if (suffix == 'n') return baseVal * 1e-9;
    if (suffix == 'u') return baseVal * 1e-6;
    
    if (suffix == 'm') {
        if (suffixStr.find("meg") == 0) {
            return baseVal * 1e6;
        }
        return baseVal * 1e-3;
    }
    if (suffix == 'k') return baseVal * 1e3;
    if (suffix == 'g') return baseVal * 1e9;
    if (suffix == 't') return baseVal * 1e12;

    return baseVal;
}

std::vector<int> FileParser::parseNodes(const std::vector<std::string>& tokens, size_t count) {
    if (tokens.size() < count + 1) {
        throw std::runtime_error("Malformed declaration: missing structural node parameters");
    }
    
    std::vector<int> mappedNodes;
    mappedNodes.reserve(count);
    for (size_t i = 1; i <= count; ++i) {
        mappedNodes.push_back(getOrAssignNode(tokens[i]));
    }
    return mappedNodes;
}

std::vector<int> FileParser::parseNodes(const std::vector<std::string>& tokens, size_t count, Netlist& netlist) {
    std::vector<int> parsed;
    for (size_t i = 1; i <= count && i < tokens.size(); ++i) {
        parsed.push_back(getOrAssignNode(tokens[i], netlist));
    }
    return parsed;
}

GateType FileParser::parseGateType(const std::string& gateStr) {
    std::string cleanStr = gateStr;
    std::transform(cleanStr.begin(), cleanStr.end(), cleanStr.begin(), [](unsigned char c) { return std::toupper(c); });

    if (cleanStr == "AND")  return GateType::AND;
    if (cleanStr == "OR")   return GateType::OR;
    if (cleanStr == "NAND") return GateType::NAND;
    if (cleanStr == "NOR")  return GateType::NOR;
    if (cleanStr == "XOR")  return GateType::XOR;
    if (cleanStr == "XNOR") return GateType::XNOR;
    if (cleanStr == "NOT")  return GateType::NOT;

    throw std::runtime_error("Unrecognized GateType token parameter: " + gateStr);
}

void FileParser::parseNetlistFile(const std::string& filename, Netlist& netlist) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Fatal Error: Could not locate netlist path descriptor: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        size_t commentIdx = line.find_first_of("*#;");
        if (commentIdx != std::string::npos) {
            line = line.substr(0, commentIdx);
        }

        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> tokens;
        while (ss >> token) {
            tokens.push_back(token);
        }

        if (tokens.empty()) continue; 

        std::string compName = tokens[0];
        char compTypeIndicator = std::tolower(static_cast<unsigned char>(compName[0]));

        try {
            // =================================================================
            // RESISTOR PARSING: Rname N1 N2 Value
            // =================================================================
            if (compTypeIndicator == 'r' && tokens.size() >= 4) {
                std::vector<int> nodes = parseNodes(tokens, 2);
                double val = parseValue(tokens[3]);
                netlist.addComponent(std::make_shared<Resistor>(compName, nodes, val));
            }
            // =================================================================
            // VOLTAGE SOURCE PARSING: Vname N1 N2 [DC] Value [AC Mag Phase]
            // =================================================================
            else if (compTypeIndicator == 'v' && tokens.size() >= 4) {
                std::vector<int> nodes = parseNodes(tokens, 2);
                double dcVal = 0.0;
                double acMag = 0.0;
                double acPhase = 0.0;
                bool acFound = false;
                size_t acIdx = 0;

                for (size_t i = 3; i < tokens.size(); ++i) {
                    std::string lowerToken = tokens[i];
                    std::transform(lowerToken.begin(), lowerToken.end(), lowerToken.begin(), [](unsigned char c){ return std::tolower(c); });
                    if (lowerToken == "ac") {
                        acFound = true;
                        acIdx = i;
                        break;
                    }
                }

                if (tokens[3] == "DC" || tokens[3] == "dc") {
                    if (tokens.size() >= 5) {
                        dcVal = parseValue(tokens[4]);
                    }
                } else if (!acFound || acIdx > 3) {
                    dcVal = parseValue(tokens[3]);
                }

                if (acFound) {
                    if (acIdx + 1 < tokens.size()) {
                        acMag = parseValue(tokens[acIdx + 1]);
                    } else {
                        acMag = 1.0; 
                    }
                    if (acIdx + 2 < tokens.size()) {
                        acPhase = parseValue(tokens[acIdx + 2]);
                    } else {
                        acPhase = 0.0; 
                    }
                }

                auto vSource = std::make_shared<VoltageSource>(compName, nodes, dcVal);
                if (acFound) {
                    vSource->setAcParameters(acMag, acPhase);
                }
                netlist.addComponent(vSource);
            }
            // =================================================================
            // CAPACITOR PARSING: Cname N1 N2 Value
            // =================================================================
            else if (compTypeIndicator == 'c' && tokens.size() >= 4) {
                std::vector<int> nodes = parseNodes(tokens, 2);
                double val = parseValue(tokens[3]);
                netlist.addComponent(std::make_shared<Capacitor>(compName, nodes, val));
            }
            // =================================================================
            // INDUCTOR PARSING: Lname N1 N2 Value
            // =================================================================
            else if (compTypeIndicator == 'l' && tokens.size() >= 4) {
                std::vector<int> nodes = parseNodes(tokens, 2);
                double val = parseValue(tokens[3]);
                netlist.addComponent(std::make_shared<Inductor>(compName, nodes, val));
            }
            // =================================================================
            // DIODE PARSING: Dname Anode Cathode
            // =================================================================
            else if (compTypeIndicator == 'd' && tokens.size() >= 3) {
                std::vector<int> nodes = parseNodes(tokens, 2);
                netlist.addComponent(std::make_shared<Diode>(compName, nodes));
            }
	    // =================================================================
            // BEHAVIORAL LOGIC GATE PARSING (Dynamic Token Inspection Layer)
            // Supports arbitrary input counts:
            //   - Mgate InA Out NOT
            //   - Mgate InA InB Out AND [vThresh vHigh vLow rOut]
            //   - Mgate InA InB InC Out OR  [vThresh vHigh vLow rOut]
            // =================================================================
            else if (compTypeIndicator == 'm' && tokens.size() >= 3) {
                // DYNAMICALLY DETECT POSITION OF THE LOGIC IDENTIFIER
                size_t gateTokenIdx = 0;
                for (size_t i = 2; i < tokens.size(); ++i) {
                    std::string t = tokens[i];
                    std::transform(t.begin(), t.end(), t.begin(), ::toupper);
                    if (t == "NOT"  || t == "AND"  || t == "OR"  || 
                        t == "NAND" || t == "NOR"  || t == "XOR" || t == "XNOR") {
                        gateTokenIdx = i;
                        break;
                    }
                }

                if (gateTokenIdx > 0) {
                    // Explicitly counts tokens representing connectivity nodes before the keyword
                    size_t nodeCount = gateTokenIdx - 1; 
                    std::vector<int> nodes = parseNodes(tokens, nodeCount);
                    GateType gateType = parseGateType(tokens[gateTokenIdx]);
                    
                    // Parameters step forward seamlessly relative to the keyword index position
                    double vThresh = (tokens.size() > gateTokenIdx + 1) ? parseValue(tokens[gateTokenIdx + 1]) : 2.5;
                    double vHigh   = (tokens.size() > gateTokenIdx + 2) ? parseValue(tokens[gateTokenIdx + 2]) : 5.0;
                    double vLow    = (tokens.size() > gateTokenIdx + 3) ? parseValue(tokens[gateTokenIdx + 3]) : 0.0;
                    double rOut    = (tokens.size() > gateTokenIdx + 4) ? parseValue(tokens[gateTokenIdx + 4]) : 10.0;

                    netlist.addComponent(std::make_shared<BehaviorLogicGateComponent>(
                        compName, nodes, gateType, vThresh, vHigh, vLow, rOut
                    ));
                }
                else {
                    // Fallback block if you later want to implement actual analog standard MOSFETs (e.g. M1 D G S NMOS)
                    // std::cerr << "Parser Alert: Non-gate 'M' identifier or model card encountered: " << compName << std::endl;
                }
            }

            // =================================================================
            // MOSFET or BEHAVIORAL LOGIC GATE PARSING (Shared 'M' Identifier)
            // Supports short form: Mgate InA InB Out AND
            // Supports long form:  Mgate InA InB Out AND vThresh vHigh vLow rOut
            // =================================================================
            else if (compTypeIndicator == 'm' && tokens.size() >= 4) {
                if (tokens.size() >= 5 && 
                    (tokens[4] == "NOT"  || tokens[4] == "not"  || 
                     tokens[4] == "AND"  || tokens[4] == "and"  || 
                     tokens[4] == "OR"   || tokens[4] == "or"   || 
                     tokens[4] == "NAND" || tokens[4] == "nand" || 
                     tokens[4] == "NOR"  || tokens[4] == "nor"  || 
                     tokens[4] == "XOR"  || tokens[4] == "xor"  || 
                     tokens[4] == "XNOR" || tokens[4] == "xnor")) {
                     
                    std::vector<int> nodes = parseNodes(tokens, 3); // Maps InA, InB, Out
                    GateType gateType = parseGateType(tokens[4]);   
                    
                    // Default logic levels if omitted from the netlist line
                    double vThresh = (tokens.size() >= 6) ? parseValue(tokens[5]) : 2.5;
                    double vHigh   = (tokens.size() >= 7) ? parseValue(tokens[6]) : 5.0;
                    double vLow    = (tokens.size() >= 8) ? parseValue(tokens[7]) : 0.0;
                    double rOut    = (tokens.size() >= 9) ? parseValue(tokens[8]) : 10.0;
                    
                    // FIXED: Allocating via make_shared to match the Netlist container definition
                    netlist.addComponent(std::make_shared<BehaviorLogicGateComponent>(
                        compName, nodes, gateType, vThresh, vHigh, vLow, rOut
                    ));
                } 
                else {
                    std::vector<int> nodes = parseNodes(tokens, 3);
                    bool isNMOS = true;
                    if (tokens.size() >= 5) {
                        std::string typeStr = tokens[4]; 
                        std::transform(typeStr.begin(), typeStr.end(), typeStr.begin(), [](unsigned char c){ return std::tolower(c); });
                        if (typeStr == "pmos") isNMOS = false;
                    }
                    // FIXED: Allocating via make_shared to preserve object vtable references
                    netlist.addComponent(std::make_shared<Mosfet>(compName, nodes, isNMOS));
                }
            }

            // =================================================================
            // BJT PARSING: Qname Collector Base Emitter [npn/pnp]
            // =================================================================
            else if (compTypeIndicator == 'q' && tokens.size() >= 4) {
                std::vector<int> nodes = parseNodes(tokens, 3); // Correctly parses and maps 3 distinct terminal indices (NC, NB, NE)
                bool isNPN = true;
                if (tokens.size() >= 5) {
                    std::string typeStr = tokens[4]; 
                    std::transform(typeStr.begin(), typeStr.end(), typeStr.begin(), [](unsigned char c){ return std::tolower(c); });
                    if (typeStr == "pnp") isNPN = false;
                }
                // FIXED: Replaced make_unique with make_shared to protect runtime vtable checks
                netlist.addComponent(std::make_shared<VbicBjtAdvanced>(compName, nodes, isNPN));
            }
        }
        catch (const std::exception& e) {
            std::cerr << "Parser Warning: Skipping malformed line entry element [" 
                      << compName << "]: " << e.what() << std::endl;
        }
    }
}
