#include "freyalab/FileParser.h"
#include "freyalab/CircuitComponents.h"
#include "freyalab/BSIM.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <cstring>
#include <functional>

FileParser::FileParser() {
    nodeMap["0"]   = NetNode{0};
    nodeMap["GND"] = NetNode{0};
}

NetNode FileParser::getOrAssignNode(const std::string& nodeName) {
    std::string upperName = nodeName;
    std::transform(upperName.begin(), upperName.end(), upperName.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

    if (upperName == "0" || upperName == "GND") {
        return NetNode{0};
    }

    auto it = nodeMap.find(upperName);
    if (it != nodeMap.end()) {
        return it->second;
    }

    nodeCounter++;
    NetNode assigned{nodeCounter};
    nodeMap[upperName] = assigned;
    return assigned;
}

NetNode FileParser::getOrAssignNode(const std::string& nodeName, Netlist& netlist) {
    std::string upperName = nodeName;
    std::transform(upperName.begin(), upperName.end(), upperName.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

    if (upperName == "0" || upperName == "GND") {
        return NetNode{0};
    }

    auto it = nodeMap.find(upperName);
    if (it != nodeMap.end()) {
        return it->second;
    }

    nodeCounter++;
    NetNode assigned{nodeCounter};
    nodeMap[upperName] = assigned;
    netlist.setNodeName(assigned, upperName);

    return assigned;
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

    std::transform(suffixStr.begin(), suffixStr.end(), suffixStr.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

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

std::vector<NetNode> FileParser::parseNodes(const std::vector<std::string>& tokens, size_t count) {
    if (tokens.size() < count + 1) {
        throw std::runtime_error("Malformed declaration: missing structural node parameters");
    }

    std::vector<NetNode> mappedNodes;
    mappedNodes.reserve(count);
    for (size_t i = 1; i <= count; ++i) {
        mappedNodes.push_back(getOrAssignNode(tokens[i]));
    }
    return mappedNodes;
}

std::vector<NetNode> FileParser::parseNodes(const std::vector<std::string>& tokens, size_t count, Netlist& netlist) {
    if (tokens.size() < count + 1) {
        throw std::runtime_error("Malformed declaration: missing structural node parameters");
    }

    std::vector<NetNode> parsed;
    parsed.reserve(count);
    for (size_t i = 1; i <= count; ++i) {
        parsed.push_back(getOrAssignNode(tokens[i], netlist));
    }
    return parsed;
}

GateType FileParser::parseGateType(const std::string& gateStr) {
    std::string cleanStr = gateStr;
    std::transform(cleanStr.begin(), cleanStr.end(), cleanStr.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

    if (cleanStr == "AND")  return GateType::AND;
    if (cleanStr == "OR")   return GateType::OR;
    if (cleanStr == "NAND") return GateType::NAND;
    if (cleanStr == "NOR")  return GateType::NOR;
    if (cleanStr == "XOR")  return GateType::XOR;
    if (cleanStr == "XNOR") return GateType::XNOR;
    if (cleanStr == "NOT")  return GateType::NOT;

    throw std::runtime_error("Unrecognized GateType token parameter: " + gateStr);
}
// Parses all key=value BSIM-CMG 112.1.0 parameters:
// 1. Populates params.rawOverrides[KEY] so any OSDI parameter is forwarded automatically.
// 2. Updates Phase-1 fallback fields for built-in execution when bsimcmg.osdi is absent.
void FileParser::parseBsimInlineParams(const std::vector<std::string>& tokens,
                                       size_t startIdx,
                                       BsimCmgParameters& params) {
    for (size_t i = startIdx; i < tokens.size(); ++i) {
        size_t eqPos = tokens[i].find('=');
        if (eqPos == std::string::npos || eqPos == 0 || eqPos + 1 >= tokens[i].size()) {
            continue;
        }

        std::string key = tokens[i].substr(0, eqPos);
        std::string valStr = tokens[i].substr(eqPos + 1);
        std::transform(key.begin(), key.end(), key.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

        double val = parseValue(valStr);

        // Forward every key=value pair directly to BSIM-CMG 112.1.0 OSDI parameter table
        params.rawOverrides[key] = val;

        // Synchronize Phase-1 built-in fallback struct fields
        if      (key == "L")       params.L       = val;
        else if (key == "TFIN")    params.TFIN    = val;
        else if (key == "HFIN")    params.HFIN    = val;
        else if (key == "WSHEET")  params.WSHEET  = val;
        else if (key == "NFIN")    params.NFIN    = static_cast<int>(std::max(1.0, val));
        else if (key == "EOT")     params.EOT     = val;
        else if (key == "VTH0")    params.VTH0    = val;
        else if (key == "NFACTOR") params.NFACTOR = val;
        else if (key == "ETA0")    params.ETA0    = val;
        else if (key == "DVT0")    params.DVT0    = val;
        else if (key == "U0")      params.U0      = val;
        else if (key == "VSAT")    params.VSAT    = val;
        else if (key == "PCLM")    params.PCLM    = val;
        else if (key == "DELTA")   params.DELTA   = val;
        else if (key == "CGS0")    params.CGS0    = val;
        else if (key == "CGD0")    params.CGD0    = val;
    }
}

// ==========================================================================
// PASS 2: RECURSIVE SUB-CIRCUIT BLOCKS FLATTENING LAYERS
// ==========================================================================
void FileParser::flattenSubcircuitInstance(
    const std::string& instanceName,
    const std::string& subcircuitTypeName,
    const std::vector<std::string>& actualNodes,
    const std::string& hierarchicalPrefix,
    std::vector<std::vector<std::string>>& globalFlattenedTokens
) {
    auto toUpper = [](std::string s) {
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
        return s;
    };

    std::string typeUpper = toUpper(subcircuitTypeName);

    auto it = subcircuitRegistry.find(typeUpper);
    if (it == subcircuitRegistry.end()) {
        throw std::runtime_error("Hierarchical Error: Subcircuit blueprint macro not found: " + subcircuitTypeName);
    }

    const auto& blueprint = it->second;
    if (actualNodes.size() != blueprint.formalPorts.size()) {
        throw std::runtime_error("Terminal mismatch on subcircuit instance '" + instanceName +
                                 "'. Expected " + std::to_string(blueprint.formalPorts.size()) +
                                 " ports, received " + std::to_string(actualNodes.size()));
    }

    std::map<std::string, std::string> portMapping;
    for (size_t i = 0; i < blueprint.formalPorts.size(); ++i) {
        portMapping[toUpper(blueprint.formalPorts.at(i))] = actualNodes.at(i);
    }

    std::string newPrefix = hierarchicalPrefix.empty() ? instanceName : hierarchicalPrefix + "_" + instanceName;

    auto resolveNode = [&](const std::string& rawNode) -> std::string {
        std::string upperNode = toUpper(rawNode);
        auto portIt = portMapping.find(upperNode);
        if (portIt != portMapping.end()) {
            return portIt->second;
        }
        if (upperNode == "0" || upperNode == "GND") {
            return "0";
        }
        return newPrefix + "_" + rawNode;
    };

    for (const auto& lineTokens : blueprint.internalLines) {
        if (lineTokens.empty()) continue;

        char firstChar = static_cast<char>(std::tolower(static_cast<unsigned char>(lineTokens.front().front())));

        if (firstChar == 'x') {
            if (lineTokens.size() < 3) {
                throw std::runtime_error("Hierarchical Error: Malformed nested subcircuit call inside " + subcircuitTypeName);
            }
            std::string nestedInstance = lineTokens.front();
            std::string nestedType = lineTokens.back();

            std::vector<std::string> nestedActualNodes;
            nestedActualNodes.reserve(lineTokens.size() - 2);
            for (size_t i = 1; i < lineTokens.size() - 1; ++i) {
                nestedActualNodes.push_back(resolveNode(lineTokens.at(i)));
            }

            flattenSubcircuitInstance(nestedInstance, nestedType, nestedActualNodes, newPrefix, globalFlattenedTokens);
        } else {
            std::vector<std::string> expandedLine = lineTokens;
            expandedLine.front() = newPrefix + "_" + lineTokens.front();

            size_t nodeCount = 2;

            if (firstChar == 'm') {
                // Dynamically detect if a logic gate keyword appears at index >= 2
                size_t gateIdx = 0;
                for (size_t k = 2; k < lineTokens.size(); ++k) {
                    std::string tokUpper = toUpper(lineTokens.at(k));
                    if (tokUpper == "NOT"  || tokUpper == "AND"  || tokUpper == "OR"   ||
                        tokUpper == "NAND" || tokUpper == "NOR"  || tokUpper == "XOR"  || tokUpper == "XNOR") {
                        gateIdx = k;
                        break;
                    }
                }
                if (gateIdx > 0) {
                    nodeCount = gateIdx - 1; // Supports arbitrary input counts + Out
                } else {
                    // Check if 4th bulk terminal precedes the model name at index 5
                    nodeCount = 3; // Drain Gate Source
                    if (lineTokens.size() >= 6 && lineTokens.at(4).find('=') == std::string::npos) {
                        std::string tok5 = toUpper(lineTokens.at(5));
                        if (tok5 == "NMOS" || tok5 == "PMOS" ||
                            tok5.find("FINFET") != std::string::npos ||
                            tok5.find("RIBBON") != std::string::npos ||
                            tok5.find("GAA")    != std::string::npos) {
                            nodeCount = 4; // Drain Gate Source Bulk
                        }
                    }
                }
            } else if (firstChar == 'q') {
                nodeCount = 3; // Collector Base Emitter
            } else if (firstChar == 'e' || firstChar == 'g') {
                nodeCount = 4; // N+ N- NC+ NC-
            }

            for (size_t i = 1; i <= nodeCount && i < expandedLine.size(); ++i) {
                expandedLine.at(i) = resolveNode(expandedLine.at(i));
            }
            globalFlattenedTokens.push_back(std::move(expandedLine));
        }
    }
}

void FileParser::parseNetlistFile(const std::string& filename, Netlist& netlist) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Fatal Error: Could not locate netlist path descriptor: " + filename);
    }

    std::string line;
    std::vector<std::vector<std::string>> primaryTokens;
    SubcircuitBlueprint* currentSubcircuitScope = nullptr;

    // -----------------------------------------------------------------
    // PASS 1: EXTRACT SUB-CIRCUIT BLUEPRINTS & REGISTERS
    // -----------------------------------------------------------------
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

        std::string cmdName = tokens.front();
        std::transform(cmdName.begin(), cmdName.end(), cmdName.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

        if (cmdName == ".SUBCKT") {
            if (tokens.size() < 3) {
                throw std::runtime_error("Malformed .SUBCKT description format.");
            }
            std::string subName = tokens.at(1);
            std::transform(subName.begin(), subName.end(), subName.begin(),
                           [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

            SubcircuitBlueprint blueprint;
            blueprint.name = subName;
            for (size_t i = 2; i < tokens.size(); ++i) {
                blueprint.formalPorts.push_back(tokens.at(i));
            }

            subcircuitRegistry[subName] = blueprint;
            currentSubcircuitScope = &subcircuitRegistry[subName];
        }
        else if (cmdName == ".ENDS") {
            currentSubcircuitScope = nullptr;
        }
        else {
            if (currentSubcircuitScope != nullptr) {
                currentSubcircuitScope->internalLines.push_back(tokens);
            } else {
                primaryTokens.push_back(tokens);
            }
        }
    }
    file.close();

    // -----------------------------------------------------------------
    // PASS 2: COMPUTE HIERARCHICAL FLATTENING
    // -----------------------------------------------------------------
    std::vector<std::vector<std::string>> fullyFlattenedTokens;

    for (const auto& tokens : primaryTokens) {
        if (tokens.empty()) continue;
        char firstChar = static_cast<char>(std::tolower(static_cast<unsigned char>(tokens.front().front())));

        if (firstChar == 'x') {
            if (tokens.size() < 3) {
                throw std::runtime_error("Malformed subcircuit call instance (X-card) configuration.");
            }
            std::string instanceName = tokens.front();
            std::string subcircuitTypeName = tokens.back();

            std::vector<std::string> actualNodes;
            actualNodes.reserve(tokens.size() - 2);
            for (size_t i = 1; i < tokens.size() - 1; ++i) {
                actualNodes.push_back(tokens.at(i));
            }

            flattenSubcircuitInstance(instanceName, subcircuitTypeName, actualNodes, "", fullyFlattenedTokens);
        } else {
            fullyFlattenedTokens.push_back(tokens);
        }
    }

    // -----------------------------------------------------------------
    // PASS 3: DUALLY MATCHED MOSFET / BSIM-CMG / GATE & BJT CONSTRUCTORS
    // -----------------------------------------------------------------
    for (const auto& tokens : fullyFlattenedTokens) {
        if (tokens.empty()) continue;

        std::string compName = tokens.front();
        char compTypeIndicator = static_cast<char>(std::tolower(static_cast<unsigned char>(compName.front())));

        try {
            // =================================================================
            // RESISTOR PARSING: Rname N1 N2 Value
            // =================================================================
            if (compTypeIndicator == 'r' && tokens.size() >= 4) {
                std::vector<NetNode> nodes = parseNodes(tokens, 2, netlist);
                double val = parseValue(tokens.at(3));
                netlist.addComponent(std::make_shared<Resistor>(compName, nodes, val));
            }
            // =================================================================
            // VOLTAGE SOURCE PARSING: Vname N1 N2 [DC] Value [AC Mag Phase]
            // =================================================================
            else if (compTypeIndicator == 'v' && tokens.size() >= 4) {
                std::vector<NetNode> nodes = parseNodes(tokens, 2, netlist);
                double dcVal = 0.0;
                double acMag = 0.0;
                double acPhase = 0.0;
                bool acFound = false;
                size_t acIdx = 0;

                for (size_t i = 3; i < tokens.size(); ++i) {
                    std::string lowerToken = tokens.at(i);
                    std::transform(lowerToken.begin(), lowerToken.end(), lowerToken.begin(),
                                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                    if (lowerToken == "ac") {
                        acFound = true;
                        acIdx = i;
                        break;
                    }
                }

                if (tokens.at(3) == "DC" || tokens.at(3) == "dc") {
                    if (tokens.size() >= 5) {
                        dcVal = parseValue(tokens.at(4));
                    }
                } else if (!acFound || acIdx > 3) {
                    dcVal = parseValue(tokens.at(3));
                }

                if (acFound) {
                    if (acIdx + 1 < tokens.size()) {
                        acMag = parseValue(tokens.at(acIdx + 1));
                    } else {
                        acMag = 1.0;
                    }
                    if (acIdx + 2 < tokens.size()) {
                        acPhase = parseValue(tokens.at(acIdx + 2));
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
                std::vector<NetNode> nodes = parseNodes(tokens, 2, netlist);
                double val = parseValue(tokens.at(3));
                netlist.addComponent(std::make_shared<Capacitor>(compName, nodes, val));
            }
	    // =================================================================
            // INDUCTOR PARSING: Lname N1 N2 Value
            // =================================================================
            else if (compTypeIndicator == 'l' && tokens.size() >= 4) {
                std::vector<NetNode> nodes = parseNodes(tokens, 2, netlist);
                double val = parseValue(tokens.at(3));
                netlist.addComponent(std::make_shared<Inductor>(compName, nodes, val));
            }
            // =================================================================
            // DIODE PARSING: Dname Anode Cathode
            // =================================================================
            else if (compTypeIndicator == 'd' && tokens.size() >= 3) {
                std::vector<NetNode> nodes = parseNodes(tokens, 2, netlist);
                netlist.addComponent(std::make_shared<Diode>(compName, nodes));
            }
            // =================================================================
            // UNIFIED 'M' PARSER:
            //   1. Behavioral Logic Gates: Mgate InA [InB ...] Out GATE [vThresh vHigh vLow rOut]
            //   2. BSIM-CMG Multi-Gate:    Mname D G S [B] NFINFET|PFINFET|NRIBBON|PRIBBON|NGAA|PGAA [L=... NFIN=...]
            //   3. Standard Level-1 FET:   Mname D G S [B] NMOS|PMOS
            // =================================================================
            else if (compTypeIndicator == 'm' && tokens.size() >= 3) {
                // Step A: Dynamically detect if this 'M' card is a Behavioral Logic Gate
                size_t gateTokenIdx = 0;
                for (size_t i = 2; i < tokens.size(); ++i) {
                    std::string t = tokens.at(i);
                    std::transform(t.begin(), t.end(), t.begin(),
                                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
                    if (t == "NOT"  || t == "AND"  || t == "OR"  ||
                        t == "NAND" || t == "NOR"  || t == "XOR" || t == "XNOR") {
                        gateTokenIdx = i;
                        break;
                    }
                }

                if (gateTokenIdx > 0) {
                    size_t nodeCount = gateTokenIdx - 1;
                    std::vector<NetNode> nodes = parseNodes(tokens, nodeCount, netlist);
                    GateType gateType = parseGateType(tokens.at(gateTokenIdx));

                    double vThresh = (tokens.size() > gateTokenIdx + 1) ? parseValue(tokens.at(gateTokenIdx + 1)) : 2.5;
                    double vHigh   = (tokens.size() > gateTokenIdx + 2) ? parseValue(tokens.at(gateTokenIdx + 2)) : 5.0;
                    double vLow    = (tokens.size() > gateTokenIdx + 3) ? parseValue(tokens.at(gateTokenIdx + 3)) : 0.0;
                    double rOut    = (tokens.size() > gateTokenIdx + 4) ? parseValue(tokens.at(gateTokenIdx + 4)) : 10.0;

                    netlist.addComponent(std::make_shared<BehaviorLogicGateComponent>(
                        compName, nodes, gateType, vThresh, vHigh, vLow, rOut
                    ));
                }
                else if (tokens.size() >= 4) {
                    // Step B: Physical Transistor (3-node {D, G, S} or 4-node {D, G, S, B})
                    std::string modelTokenUpper = "NMOS";
                    size_t modelIdx = 4;

                    // Check if a 4th Bulk node precedes the model keyword at index 5
                    if (tokens.size() >= 6 && tokens.at(4).find('=') == std::string::npos) {
                        std::string tok5 = tokens.at(5);
                        std::transform(tok5.begin(), tok5.end(), tok5.begin(),
                                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
                        if (tok5 == "NMOS" || tok5 == "PMOS" ||
                            tok5.find("FINFET") != std::string::npos ||
                            tok5.find("RIBBON") != std::string::npos ||
                            tok5.find("GAA")    != std::string::npos) {
                            modelIdx = 5;
                        }
                    }

                    // Parse 4 nodes {D, G, S, B} if Bulk is present, otherwise 3 nodes {D, G, S}
                    size_t terminalCount = (modelIdx == 5) ? 4 : 3;
                    std::vector<NetNode> nodes = parseNodes(tokens, terminalCount, netlist);

                    if (tokens.size() > modelIdx && tokens.at(modelIdx).find('=') == std::string::npos) {
                        modelTokenUpper = tokens.at(modelIdx);
                        std::transform(modelTokenUpper.begin(), modelTokenUpper.end(), modelTokenUpper.begin(),
                                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
                    }

                    // Dispatch BSIM-CMG Multi-Gate FET (FinFET, RibbonFET, GAA-FET)
                    if (modelTokenUpper == "NFINFET"    || modelTokenUpper == "PFINFET"    ||
                        modelTokenUpper == "NRIBBON"    || modelTokenUpper == "PRIBBON"    ||
                        modelTokenUpper == "NRIBBONFET" || modelTokenUpper == "PRIBBONFET" ||
                        modelTokenUpper == "NGAA"       || modelTokenUpper == "PGAA"       ||
                        modelTokenUpper == "NGAAFET"    || modelTokenUpper == "PGAAFET") {

                        BsimCmgParameters bsimParams;
                        bsimParams.isNChannel = (modelTokenUpper.front() == 'N');

                        ComponentType compSubType = ComponentType::FINFET;
                        if (modelTokenUpper.find("RIBBON") != std::string::npos) {
                            compSubType = ComponentType::RIBBONFET;
                            bsimParams.geomod = MultiGateGeometry::GAA_NANOSHEET;
                        } else if (modelTokenUpper.find("GAA") != std::string::npos) {
                            compSubType = ComponentType::GAA_FET;
                            bsimParams.geomod = MultiGateGeometry::GAA_NANOSHEET;
                        } else {
                            compSubType = ComponentType::FINFET;
                            bsimParams.geomod = MultiGateGeometry::FINFET_TRIGATE;
                        }

                        parseBsimInlineParams(tokens, modelIdx + 1, bsimParams);

                        netlist.addComponent(std::make_shared<BsimCmgFet>(
                            compName, nodes, compSubType, bsimParams
                        ));
                    }
                    else {
                        // Dispatch Standard Level-1 Planar MOSFET (NMOS / PMOS)
                        bool isNMOS = (modelTokenUpper != "PMOS");
                        netlist.addComponent(std::make_shared<Mosfet>(compName, nodes, isNMOS));
                    }
                }
            }
            // =================================================================
            // BJT PARSING: Qname Collector Base Emitter [npn/pnp]
            // =================================================================
            else if (compTypeIndicator == 'q' && tokens.size() >= 4) {
                std::vector<NetNode> nodes = parseNodes(tokens, 3, netlist);
                bool isNPN = true;
                if (tokens.size() >= 5) {
                    std::string typeStr = tokens.at(4);
                    std::transform(typeStr.begin(), typeStr.end(), typeStr.begin(),
                                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                    if (typeStr == "pnp") isNPN = false;
                }
                netlist.addComponent(std::make_shared<VbicBjtAdvanced>(compName, nodes, isNPN));
            }
        }
        catch (const std::exception& e) {
            std::cerr << "Parser Warning: Skipping malformed line entry element ["
                      << compName << "]: " << e.what() << std::endl;
        }
    }
}

