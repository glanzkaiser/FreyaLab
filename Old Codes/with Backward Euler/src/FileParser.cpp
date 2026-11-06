#include "freyalab/FileParser.h"
#include "freyalab/CircuitComponents.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <set>
#include <map>

double FileParser::parseValue(const std::string& valueStr) {
    if (valueStr.empty()) return 0.0;
    
    size_t i = 0;
    bool processingExponent = false;
    
    while (i < valueStr.size()) {
        char c = valueStr[i];
        if (std::isdigit(c) || c == '.') {
            i++;
        } else if ((c == '-' || c == '+')) {
            if (i == 0 || processingExponent) {
                i++;
                processingExponent = false;
            } else {
                break;
            }
        } else if ((c == 'e' || c == 'E')) {
            processingExponent = true;
            i++;
        } else {
            break;
        }
    }
    
    std::string numPart = valueStr.substr(0, i);
    std::string suffix = valueStr.substr(i);
    double baseVal = 0.0;
    try { 
        baseVal = std::stod(numPart); 
    } catch (...) { 
        return 0.0; 
    }
    if (suffix.empty()) return baseVal;
    
    std::string ls = suffix; 
    std::transform(ls.begin(), ls.end(), ls.begin(), ::tolower);
    
    if (ls.size() >= 3 && ls.substr(0, 3) == "meg") return baseVal * 1e6;
    
    char scale = ls[0];
    switch (scale) {
        case 'f': return baseVal * 1e-15; 
        case 'p': return baseVal * 1e-12; 
        case 'n': return baseVal * 1e-9;
        case 'u': return baseVal * 1e-6; 
        case 'm': return baseVal * 1e-3; 
        case 'k': return baseVal * 1e3;
        case 'g': return baseVal * 1e9; 
        case 't': return baseVal * 1e12; 
        default: return baseVal;
    }
}

std::string FileParser::cleanAndTokenize(const std::string& line) {
    std::string clean = "";
    for (char c : line) {
        if (c == '(' || c == ')' || c == ',' || c == '=') {
            clean += ' ';
        } else {
            clean += c;
        }
    }
    return clean;
}

bool FileParser::parseNetlist(const std::string& filename, std::vector<std::shared_ptr<Component>>& outComponents, int& numNodes) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    
    outComponents.clear();
    
    std::map<std::string, int> nodeMap;
    nodeMap["0"] = 0;
    nodeMap["gnd"] = 0;
    int nextNodeIdx = 1;

    auto getOrCreateNode = [&](const std::string& nodeStr) {
        std::string lowStr = nodeStr;
        std::transform(lowStr.begin(), lowStr.end(), lowStr.begin(), ::tolower);
        if (nodeMap.find(lowStr) == nodeMap.end()) {
            nodeMap[lowStr] = nextNodeIdx++;
        }
        return nodeMap[lowStr];
    };
    
    std::string line; 
    while (std::getline(file, line)) {
        line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](unsigned char ch) {
            return !std::isspace(ch);
        }));

        if (line.empty() || line[0] == '*' || line[0] == '#') continue;
        
        std::stringstream ss(cleanAndTokenize(line));
        std::string name;
        if (!(ss >> name)) continue;
        
        std::vector<std::string> tokens;
        std::string t;
        while (ss >> t) tokens.push_back(t);
        
        char compType = std::toupper(name[0]);
        
        try {
            // --- 3-Terminal Device (BJT Mapping) ---
            if (compType == 'Q') {
                if (tokens.size() < 3) continue;
                int nC = getOrCreateNode(tokens[0]);
                int nB = getOrCreateNode(tokens[1]);
                int nE = getOrCreateNode(tokens[2]);
                
                bool isNpn = true;
                double Bf = 100.0, Br = 1.0, Is = 1e-14;
                
                for (size_t i = 3; i < tokens.size(); ++i) {
                    std::string pKey = tokens[i];
                    std::transform(pKey.begin(), pKey.end(), pKey.begin(), ::toupper);
                    if (pKey == "PNP") isNpn = false;
                    else if (pKey == "NPN") isNpn = true;
                    else {
                        std::transform(pKey.begin(), pKey.end(), pKey.begin(), ::tolower);
                        if (pKey == "bf" && i + 1 < tokens.size()) Bf = parseValue(tokens[++i]);
                        else if (pKey == "br" && i + 1 < tokens.size()) Br = parseValue(tokens[++i]);
                        else if (pKey == "is" && i + 1 < tokens.size()) Is = parseValue(tokens[++i]);
                    }
                }
                outComponents.push_back(std::make_shared<BJT>(name, nC, nB, nE, isNpn, Bf, Br, Is));
                continue;
            }
            
            // --- 2-Terminal Devices ---
            if (tokens.size() < 3) continue; 
            int nodeA = getOrCreateNode(tokens[0]);
            int nodeB = getOrCreateNode(tokens[1]);
            
            if (compType == 'R') {
                outComponents.push_back(std::make_shared<Resistor>(name, nodeA, nodeB, parseValue(tokens[2])));
            } else if (compType == 'C') {
                outComponents.push_back(std::make_shared<Capacitor>(name, nodeA, nodeB, parseValue(tokens[2])));
            } else if (compType == 'L') {
                outComponents.push_back(std::make_shared<Inductor>(name, nodeA, nodeB, parseValue(tokens[2])));
            } else if (compType == 'D') {
                double Is = 1e-14, Vt = 0.026;
                for (size_t i = 2; i < tokens.size(); ++i) {
                    std::string pKey = tokens[i];
                    std::transform(pKey.begin(), pKey.end(), pKey.begin(), ::tolower);
                    if (pKey == "is" && i + 1 < tokens.size()) Is = parseValue(tokens[++i]);
                    else if (pKey == "vt" && i + 1 < tokens.size()) Vt = parseValue(tokens[++i]);
                }
                outComponents.push_back(std::make_shared<Diode>(name, nodeA, nodeB, Is, Vt));
            } else if (compType == 'V') {
                std::string specType = tokens[2];
                std::transform(specType.begin(), specType.end(), specType.begin(), ::toupper);
                
                if (std::isdigit(specType[0]) || specType[0] == '-' || specType[0] == '+') {
                    outComponents.push_back(std::make_shared<VoltageSource>(name, nodeA, nodeB, parseValue(specType)));
                } else if (specType == "DC" && tokens.size() > 3) {
                    outComponents.push_back(std::make_shared<VoltageSource>(name, nodeA, nodeB, parseValue(tokens[3])));
                } else if (specType == "AC") {
                    double mag = (tokens.size() > 3) ? parseValue(tokens[3]) : 1.0;
                    outComponents.push_back(std::make_shared<VoltageSource>(name, nodeA, nodeB, SourceType::AC, 0.0, mag, 0.0));
                } else if (specType == "PULSE") {
                    double v1 = tokens.size() > 3 ? parseValue(tokens[3]) : 0.0;
                    double v2 = tokens.size() > 4 ? parseValue(tokens[4]) : 1.0;
                    double td = tokens.size() > 5 ? parseValue(tokens[5]) : 0.0;
                    double tr = tokens.size() > 6 ? parseValue(tokens[6]) : 1e-9;
                    double tf = tokens.size() > 7 ? parseValue(tokens[7]) : 1e-9;
                    double tw = tokens.size() > 8 ? parseValue(tokens[8]) : 0.5;
                    double tp = tokens.size() > 9 ? parseValue(tokens[9]) : 1.0;
                    outComponents.push_back(std::make_shared<VoltageSource>(name, nodeA, nodeB, SourceType::PULSE, v1, v2, td, tr, tf, tw, tp));
                }
            }
        } catch (...) {
            std::cerr << "[Warning] Netlist line parsing failure encountered on line: " << line << "\n";
            continue;
        }
    }
    
    numNodes = nextNodeIdx; 
    return true;
}
