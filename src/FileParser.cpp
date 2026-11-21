#include "freyalab/FileParser.h"
#include "freyalab/BehavioralSource.h"
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
#include <cmath>
#include <filesystem>

namespace {

class SpiceExprParser {
public:
    SpiceExprParser(const std::string& src,
                    const std::map<std::string, double>& params,
                    std::function<double(const std::string&)> numParser)
        : s(src), pos(0), paramMap(params), parseEng(std::move(numParser)) {}

    double parse() {
        return parseAddSub();
    }

private:
    std::string s;
    size_t pos;
    const std::map<std::string, double>& paramMap;
    std::function<double(const std::string&)> parseEng;

    void skipWs() {
        while (pos < s.size() && (std::isspace(static_cast<unsigned char>(s[pos])) ||
                                  s[pos] == '{' || s[pos] == '}' || s[pos] == '\'')) {
            ++pos;
        }
    }

    double parseAddSub() {
        double left = parseMulDiv();
        while (true) {
            skipWs();
            if (pos < s.size() && s[pos] == '+') { ++pos; left += parseMulDiv(); }
            else if (pos < s.size() && s[pos] == '-') { ++pos; left -= parseMulDiv(); }
            else break;
        }
        return left;
    }

    double parseMulDiv() {
        double left = parsePower();
        while (true) {
            skipWs();
            if (pos < s.size() && s[pos] == '*') { ++pos; left *= parsePower(); }
            else if (pos < s.size() && s[pos] == '/') {
                ++pos;
                double denom = parsePower();
                left = (std::abs(denom) > 1e-30) ? (left / denom) : 0.0;
            } else break;
        }
        return left;
    }

    double parsePower() {
        double base = parseUnary();
        skipWs();
        if (pos + 1 < s.size() && s[pos] == '*' && s[pos + 1] == '*') {
            pos += 2;
            return std::pow(base, parsePower());
        }
        if (pos < s.size() && s[pos] == '^') {
            ++pos;
            return std::pow(base, parsePower());
        }
        return base;
    }

    double parseUnary() {
        skipWs();
        if (pos < s.size() && s[pos] == '+') { ++pos; return parseUnary(); }
        if (pos < s.size() && s[pos] == '-') { ++pos; return -parseUnary(); }
        return parsePrimary();
    }

    double parsePrimary() {
        skipWs();
        if (pos >= s.size()) return 0.0;

        if (s[pos] == '(') {
            ++pos;
            double v = parseAddSub();
            skipWs();
            if (pos < s.size() && s[pos] == ')') ++pos;
            return v;
        }

        // Numeric literal with SPICE engineering suffix (e.g., 14n, 1.8, 2.5e-6)
        if (std::isdigit(static_cast<unsigned char>(s[pos])) || s[pos] == '.') {
            size_t start = pos;
            while (pos < s.size() && (std::isalnum(static_cast<unsigned char>(s[pos])) || s[pos] == '.')) {
                // Handle signed scientific notation exponents such as 1e-9 or 1e+3
                if ((s[pos] == 'e' || s[pos] == 'E') && pos + 1 < s.size() &&
                    (s[pos + 1] == '+' || s[pos + 1] == '-')) {
                    pos += 2;
                    continue;
                }
                ++pos;
            }
            return parseEng(s.substr(start, pos - start));
        }

        // Parameter identifier or math function (SQRT, EXP, LOG, LN, LOG10, ABS, MIN, MAX)
        if (std::isalpha(static_cast<unsigned char>(s[pos])) || s[pos] == '_') {
            size_t start = pos;
            while (pos < s.size() && (std::isalnum(static_cast<unsigned char>(s[pos])) || s[pos] == '_')) {
                ++pos;
            }
            std::string ident = s.substr(start, pos - start);
            std::transform(ident.begin(), ident.end(), ident.begin(),
                           [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

            skipWs();
            if (pos < s.size() && s[pos] == '(') {
                ++pos;
                double arg1 = parseAddSub();
                double arg2 = 0.0;
                bool hasArg2 = false;
                skipWs();
                if (pos < s.size() && s[pos] == ',') {
                    ++pos;
                    arg2 = parseAddSub();
                    hasArg2 = true;
                }
                skipWs();
                if (pos < s.size() && s[pos] == ')') ++pos;

                if (ident == "SQRT") return std::sqrt(std::max(0.0, arg1));
                if (ident == "EXP")  return std::exp(std::min(arg1, 80.0));
                if (ident == "LN" || ident == "LOG") return std::log(std::max(1e-30, arg1));
                if (ident == "LOG10") return std::log10(std::max(1e-30, arg1));
                if (ident == "ABS")  return std::abs(arg1);
                if (ident == "MIN" && hasArg2) return std::min(arg1, arg2);
                if (ident == "MAX" && hasArg2) return std::max(arg1, arg2);
                return arg1;
            }

            auto it = paramMap.find(ident);
            if (it != paramMap.end()) {
                return it->second;
            }
            return 0.0;
        }

        return 0.0;
    }
};
// Menggabungkan token berurutan jika berada di dalam kurung kurawal { ... } yang mengandung spasi,
// lalu memajukan indeks `idx` ke token tepat setelah ekspresi selesai.
std::string extractBracedExpr(const std::vector<std::string>& tokens, size_t& idx) {
    if (idx >= tokens.size()) return "";

    std::string combined = tokens[idx];
    int braceDepth = 0;
    for (char c : combined) {
        if (c == '{') ++braceDepth;
        else if (c == '}') --braceDepth;
    }
    ++idx;

    while (braceDepth > 0 && idx < tokens.size()) {
        combined += " " + tokens[idx];
        for (char c : tokens[idx]) {
            if (c == '{') ++braceDepth;
            else if (c == '}') --braceDepth;
        }
        ++idx;
    }

    return combined;
}

} // namespace



void FileParser::collectExpandedLines(const std::string& filepath,
                                      const std::string& targetCorner,
                                      std::vector<std::string>& outLogicalLines,
                                      int depth) {
    if (depth > 32) {
        throw std::runtime_error("PDK Error: Maximum .INCLUDE/.LIB recursion depth (32) exceeded.");
    }

    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("PDK Error: Unable to open included file: " + filepath);
    }

    namespace fs = std::filesystem;
    fs::path baseDir = fs::path(filepath).parent_path();

    bool insideLibFile = !targetCorner.empty();
    bool activeCornerMatch = !insideLibFile; // True if reading a normal file or .INCLUDE
    std::string upperTargetCorner = targetCorner;
    std::transform(upperTargetCorner.begin(), upperTargetCorner.end(), upperTargetCorner.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

    std::string rawLine;
    while (std::getline(file, rawLine)) {
        // 1. Find the first non-whitespace character
        size_t first = rawLine.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;

        // 2. Treat '*' as a comment ONLY if it is the first non-whitespace character
        if (rawLine[first] == '*') continue;

        // 3. Strip inline comments starting with ';' or '//'
        size_t semiIdx = rawLine.find(';', first);
        size_t slashIdx = rawLine.find("//", first);
        size_t inlineCommentIdx = std::min(semiIdx, slashIdx);

        if (inlineCommentIdx != std::string::npos) {
            rawLine = rawLine.substr(0, inlineCommentIdx);
        }

        // 4. Trim leading and trailing whitespace
        size_t last = rawLine.find_last_not_of(" \t\r\n");
        if (last == std::string::npos || last < first) continue;
        std::string trimmed = rawLine.substr(first, last - first + 1);

        // Merge '+' continuation lines onto the preceding active statement
        if (trimmed.front() == '+') {
            if (activeCornerMatch && !outLogicalLines.empty()) {
                outLogicalLines.back() += " " + trimmed.substr(1);
            }
            continue;
        }

        std::stringstream ss(trimmed);
        std::string directive;
        ss >> directive;
        std::string dirUpper = directive;
        std::transform(dirUpper.begin(), dirUpper.end(), dirUpper.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

        if (dirUpper == ".LIB") {
            std::string arg1, arg2;
            ss >> arg1;
            if (ss >> arg2) {
                // Calling syntax: .LIB "sky130.lib.spice" tt
                arg1.erase(std::remove(arg1.begin(), arg1.end(), '"'), arg1.end());
                arg1.erase(std::remove(arg1.begin(), arg1.end(), '\''), arg1.end());
                fs::path subPath = fs::path(arg1).is_absolute() ? fs::path(arg1) : (baseDir / arg1);
                collectExpandedLines(subPath.string(), arg2, outLogicalLines, depth + 1);
            } else if (insideLibFile) {
                // Section header inside library file: .LIB tt
                std::string secUpper = arg1;
                std::transform(secUpper.begin(), secUpper.end(), secUpper.begin(),
                               [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
                activeCornerMatch = (secUpper == upperTargetCorner);
            }
        } else if (dirUpper == ".ENDL") {
            if (insideLibFile) {
                activeCornerMatch = false;
            }
        } else if (dirUpper == ".INCLUDE" || dirUpper == ".INC") {
            if (activeCornerMatch) {
                std::string incFile;
                ss >> incFile;
                incFile.erase(std::remove(incFile.begin(), incFile.end(), '"'), incFile.end());
                incFile.erase(std::remove(incFile.begin(), incFile.end(), '\''), incFile.end());
                fs::path subPath = fs::path(incFile).is_absolute() ? fs::path(incFile) : (baseDir / incFile);
                collectExpandedLines(subPath.string(), "", outLogicalLines, depth + 1);
            }
        } else if (activeCornerMatch) {
            outLogicalLines.push_back(std::move(trimmed));
        }
    }
}

double FileParser::evaluateExpression(const std::string& exprStr) const {
    if (exprStr.empty()) return 0.0;
    SpiceExprParser evaluator(exprStr, paramRegistry,
                              [this](const std::string& tok) { return this->parseValue(tok); });
    return evaluator.parse();
}

void FileParser::parseParamDirective(const std::vector<std::string>& tokens) {
    // Syntax: .PARAM VDD=1.8 LMIN=14n W_EFF={4 * 67n}
    std::string joined;
    for (size_t i = 1; i < tokens.size(); ++i) {
        joined += tokens[i] + " ";
    }

    // Split each KEY = EXPR pair while preserving braced { ... } and parenthesized ( ... ) blocks
    size_t i = 0;
    while (i < joined.size()) {
        while (i < joined.size() && std::isspace(static_cast<unsigned char>(joined[i]))) ++i;
        if (i >= joined.size()) break;

        size_t eqPos = joined.find('=', i);
        if (eqPos == std::string::npos) break;

        std::string key = joined.substr(i, eqPos - i);
        key.erase(std::remove_if(key.begin(), key.end(),
                                 [](unsigned char c) { return std::isspace(c); }),
                  key.end());
        std::transform(key.begin(), key.end(), key.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

        size_t valStart = eqPos + 1;
        while (valStart < joined.size() && std::isspace(static_cast<unsigned char>(joined[valStart]))) ++valStart;

        int braceDepth = 0;
        int parenDepth = 0;
        size_t valEnd = valStart;
        while (valEnd < joined.size()) {
            char c = joined[valEnd];
            if (c == '{') ++braceDepth;
            else if (c == '}') --braceDepth;
            else if (c == '(') ++parenDepth;
            else if (c == ')') --parenDepth;
            else if (std::isspace(static_cast<unsigned char>(c)) && braceDepth == 0 && parenDepth == 0) {
                // Stop if the next whitespace-delimited token begins a new KEY= assignment
                size_t lookahead = valEnd;
                while (lookahead < joined.size() && std::isspace(static_cast<unsigned char>(joined[lookahead]))) ++lookahead;
                size_t nextEq = joined.find('=', lookahead);
                size_t nextSp = joined.find_first_of(" \t", lookahead);
                if (nextEq != std::string::npos && (nextSp == std::string::npos || nextEq < nextSp)) {
                    break;
                }
            }
            ++valEnd;
        }

        std::string expr = joined.substr(valStart, valEnd - valStart);
        if (!key.empty() && !expr.empty()) {
            paramRegistry[key] = evaluateExpression(expr);
        }
        i = valEnd;
    }
}

void FileParser::parseModelDirective(const std::vector<std::string>& tokens) {
    // SPICE Syntax: .MODEL <ModelName> <ModelType> (PARAM1=VAL1 PARAM2=VAL2 ...)
    if (tokens.size() < 3) return;

    std::string mName = tokens[1];
    std::string mType = tokens[2];
    std::transform(mName.begin(), mName.end(), mName.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    std::transform(mType.begin(), mType.end(), mType.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

    // Strip attached opening parenthesis if written as: .MODEL NFET NMOS(VTH0=0.4 ...)
    size_t openParen = mType.find('(');
    std::string trailingFromType;
    if (openParen != std::string::npos) {
        trailingFromType = mType.substr(openParen + 1);
        mType = mType.substr(0, openParen);
    }

    ModelCardBlueprint card;
    card.modelName = mName;
    card.modelType = mType;

    // Inherit existing parameters if the model name was previously declared (model binning/overrides)
    auto existing = modelRegistry.find(mName);
    if (existing != modelRegistry.end()) {
        card.parameters = existing->second.parameters;
    }

    std::vector<std::string> paramTokens;
    if (!trailingFromType.empty()) paramTokens.push_back(trailingFromType);
    for (size_t i = 3; i < tokens.size(); ++i) {
        paramTokens.push_back(tokens[i]);
    }

    for (std::string tok : paramTokens) {
        tok.erase(std::remove(tok.begin(), tok.end(), '('), tok.end());
        tok.erase(std::remove(tok.begin(), tok.end(), ')'), tok.end());
        size_t eqPos = tok.find('=');
        if (eqPos == std::string::npos || eqPos == 0 || eqPos + 1 >= tok.size()) continue;

        std::string key = tok.substr(0, eqPos);
        std::string valExpr = tok.substr(eqPos + 1);
        std::transform(key.begin(), key.end(), key.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

        card.parameters[key] = evaluateExpression(valExpr);
    }

    modelRegistry[mName] = std::move(card);
}

void FileParser::applyModelCardToBsim(const std::string& modelName, BsimCmgParameters& params) const {
    std::string upperName = modelName;
    std::transform(upperName.begin(), upperName.end(), upperName.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

    auto it = modelRegistry.find(upperName);
    if (it == modelRegistry.end()) return;

    const auto& card = it->second;
    if (card.modelType == "PMOS" || card.modelType == "PFINFET") {
        params.isNChannel = false;
    } else if (card.modelType == "NMOS" || card.modelType == "NFINFET") {
        params.isNChannel = true;
    }

    for (const auto& kv : card.parameters) {
        const std::string& key = kv.first;
        double val = kv.second;
        params.rawOverrides[key] = val;

        if      (key == "TYPE")    params.isNChannel = (val >= 0.0);
        else if (key == "L")       params.L          = val;
        else if (key == "TFIN")    params.TFIN       = val;
        else if (key == "HFIN")    params.HFIN       = val;
        else if (key == "WSHEET")  params.WSHEET     = val;
        else if (key == "NFIN")    params.NFIN       = static_cast<int>(std::max(1.0, val));
        else if (key == "EOT")     params.EOT        = val;
        else if (key == "VTH0")    params.VTH0       = val;
        else if (key == "NFACTOR") params.NFACTOR    = val;
        else if (key == "ETA0")    params.ETA0       = val;
        else if (key == "DVT0")    params.DVT0       = val;
        else if (key == "U0")      params.U0         = val;
        else if (key == "VSAT")    params.VSAT       = val;
        else if (key == "PCLM")    params.PCLM       = val;
        else if (key == "DELTA")   params.DELTA      = val;
        else if (key == "CGS0")    params.CGS0       = val;
        else if (key == "CGD0")    params.CGD0       = val;
    }
}



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

double FileParser::parseValue(const std::string& valStr) const {
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

        if (key == "TYPE") {
            std::string upperVal = valStr;
            std::transform(upperVal.begin(), upperVal.end(), upperVal.begin(),
                           [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
            params.isNChannel = (upperVal != "PMOS" && evaluateExpression(valStr) >= 0.0);
            continue;
        }

        double val = evaluateExpression(valStr);

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
                    auto isKnownModelToken = [&](const std::string& tok) {
                        return tok == "NMOS" || tok == "PMOS" ||
                               tok.find("FINFET") != std::string::npos ||
                               tok.find("RIBBON") != std::string::npos ||
                               tok.find("GAA")    != std::string::npos ||
                               modelRegistry.find(tok) != modelRegistry.end();
                    };

                    nodeCount = 3; // Default: Drain Gate Source
                    if (lineTokens.size() >= 6 && lineTokens.at(4).find('=') == std::string::npos) {
                        std::string tok4 = toUpper(lineTokens.at(4));
                        std::string tok5 = toUpper(lineTokens.at(5));
                        // 4 terminals {D, G, S, B} if tok4 is a node (not a model) and tok5 is the model name
                        if (!isKnownModelToken(tok4) && tok5.find('=') == std::string::npos) {
                            nodeCount = 4;
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
    std::vector<std::vector<std::string>> primaryTokens;
    SubcircuitBlueprint* currentSubcircuitScope = nullptr;

    // -----------------------------------------------------------------
    // PASS 1: EXPAND .INCLUDE / .LIB CORNERS, MERGE '+', & EXTRACT .SUBCKT
    // -----------------------------------------------------------------
    std::vector<std::string> logicalLines;
    collectExpandedLines(filename, "", logicalLines, 0);

    for (const auto& logicalLine : logicalLines) {
        std::stringstream ss(logicalLine);
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
	else if (cmdName == ".PARAM") {
            // Evaluate parameters early so .MODEL expressions and Pass 2 can use them
            parseParamDirective(tokens);
        }
        else if (cmdName == ".MODEL") {
            // Register models early so Pass 2 (flattenSubcircuitInstance) knows all model names
            parseModelDirective(tokens);
        }
        else {
            if (currentSubcircuitScope != nullptr) {
                currentSubcircuitScope->internalLines.push_back(tokens);
            } else {
                primaryTokens.push_back(tokens);
            }
        }
    }

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
            // SIMULATION & PDK DIRECTIVES: .PARAM, .MODEL, .OP, .DC, .AC, .TRAN
            // =================================================================
            if (compTypeIndicator == '.') {
                std::string directive = compName;
                std::transform(directive.begin(), directive.end(), directive.begin(),
                               [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
		if (directive == ".TEMP" && tokens.size() >= 2) {
		    size_t tempIdx = 1;
                    double tempC = evaluateExpression(extractBracedExpr(tokens, tempIdx));
                    netlist.setTemperatureCelsius(tempC);
		} else if (directive == ".NOISE" && tokens.size() >= 6) {
                    // Syntax: .NOISE V(outPos [, outNeg]) inputSource DEC|OCT|LIN pts fstart fstop
                    AnalysisDirective cmd;
                    cmd.kind = AnalysisDirective::Kind::NOISE;

                    // Merge V(...) tokens if spaces appear inside parentheses, e.g., V(out+ , out-)
                    size_t cursor = 1;
                    std::string outSpec = tokens.at(cursor++);
                    while (outSpec.find('(') != std::string::npos &&
                           outSpec.find(')') == std::string::npos &&
                           cursor < tokens.size()) {
                        outSpec += tokens.at(cursor++);
                    }

                    // Extract node name(s) inside V(...)
                    size_t lParen = outSpec.find('(');
                    size_t rParen = outSpec.rfind(')');
                    if (lParen != std::string::npos && rParen != std::string::npos && rParen > lParen) {
                        std::string inside = outSpec.substr(lParen + 1, rParen - lParen - 1);
                        size_t commaPos = inside.find(',');
                        if (commaPos != std::string::npos) {
                            std::string posStr = inside.substr(0, commaPos);
                            std::string negStr = inside.substr(commaPos + 1);
                            posStr.erase(std::remove_if(posStr.begin(), posStr.end(), ::isspace), posStr.end());
                            negStr.erase(std::remove_if(negStr.begin(), negStr.end(), ::isspace), negStr.end());
                            cmd.noiseOutPos = getOrAssignNode(posStr, netlist);
                            cmd.noiseOutNeg = getOrAssignNode(negStr, netlist);
                        } else {
                            inside.erase(std::remove_if(inside.begin(), inside.end(), ::isspace), inside.end());
                            cmd.noiseOutPos = getOrAssignNode(inside, netlist);
                            cmd.noiseOutNeg = NetNode{0}; // Ground reference
                        }
                    }

                    if (cursor + 4 < tokens.size()) {
                        cmd.noiseInputSource = tokens.at(cursor++);

                        std::string scale = tokens.at(cursor++);
                        std::transform(scale.begin(), scale.end(), scale.begin(),
                                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
                        cmd.acScaleType = scale;

                        int pts = std::max(2, static_cast<int>(evaluateExpression(extractBracedExpr(tokens, cursor))));
                        cmd.acStartFreq = std::max(1e-3, evaluateExpression(extractBracedExpr(tokens, cursor)));
                        cmd.acStopFreq  = std::max(cmd.acStartFreq, evaluateExpression(extractBracedExpr(tokens, cursor)));

                        if (scale == "DEC") {
                            double decades = std::max(1.0, std::log10(cmd.acStopFreq / cmd.acStartFreq));
                            cmd.acPoints = static_cast<int>(std::round(pts * decades)) + 1;
                        } else if (scale == "OCT") {
                            double octaves = std::max(1.0, std::log2(cmd.acStopFreq / cmd.acStartFreq));
                            cmd.acPoints = static_cast<int>(std::round(pts * octaves)) + 1;
                        } else {
                            cmd.acPoints = pts;
                        }

                        netlist.addAnalysisDirective(cmd);
                    }
                } else if (directive == ".OP") {
                    AnalysisDirective cmd;
                    cmd.kind = AnalysisDirective::Kind::DC_OP;
                    netlist.addAnalysisDirective(cmd);
                } else if (directive == ".DC" && tokens.size() >= 5) {
                    AnalysisDirective cmd;
                    cmd.kind = AnalysisDirective::Kind::DC_SWEEP;
                    cmd.dcSourceName = tokens.at(1);
                    cmd.dcStart = parseValue(tokens.at(2));
                    cmd.dcStop  = parseValue(tokens.at(3));
                    cmd.dcStep  = parseValue(tokens.at(4));
                    netlist.addAnalysisDirective(cmd);
                } else if (directive == ".AC" && tokens.size() >= 5) {
                    AnalysisDirective cmd;
                    cmd.kind = AnalysisDirective::Kind::AC_SWEEP;
                    std::string scale = tokens.at(1);
                    std::transform(scale.begin(), scale.end(), scale.begin(),
                                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
                    cmd.acScaleType = scale;
                    int pts = std::max(2, static_cast<int>(parseValue(tokens.at(2))));
                    cmd.acStartFreq = std::max(1e-3, parseValue(tokens.at(3)));
                    cmd.acStopFreq  = std::max(cmd.acStartFreq, parseValue(tokens.at(4)));

                    if (scale == "DEC") {
                        double decades = std::max(1.0, std::log10(cmd.acStopFreq / cmd.acStartFreq));
                        cmd.acPoints = static_cast<int>(std::round(pts * decades)) + 1;
                    } else if (scale == "OCT") {
                        double octaves = std::max(1.0, std::log2(cmd.acStopFreq / cmd.acStartFreq));
                        cmd.acPoints = static_cast<int>(std::round(pts * octaves)) + 1;
                    } else {
                        cmd.acPoints = pts;
                    }
                    netlist.addAnalysisDirective(cmd);
                } else if (directive == ".TRAN" && tokens.size() >= 3) {
                    AnalysisDirective cmd;
                    cmd.kind = AnalysisDirective::Kind::TRANSIENT;
                    cmd.tranStep = std::max(1e-15, parseValue(tokens.at(1)));
                    cmd.tranStop = std::max(cmd.tranStep, parseValue(tokens.at(2)));
                    if (tokens.size() >= 4) cmd.tranStart = std::max(0.0, parseValue(tokens.at(3)));
                    if (tokens.size() >= 5) cmd.tranMaxDt = std::max(cmd.tranStep, parseValue(tokens.at(4)));
                    netlist.addAnalysisDirective(cmd);
                }
            }
	    // =================================================================
            // RESISTOR PARSING: Rname N1 N2 Value_or_{Expr} [TC1=val_or_{Expr}]
            // =================================================================
            else if (compTypeIndicator == 'r' && tokens.size() >= 4) {
                std::vector<NetNode> nodes = parseNodes(tokens, 2, netlist);
                size_t valIdx = 3;
                double val = evaluateExpression(extractBracedExpr(tokens, valIdx));

                auto resistor = std::make_shared<Resistor>(compName, nodes, val);

                // Parse optional trailing parameters (e.g., TC1=3.9e-3 or TC1={expr})
                while (valIdx < tokens.size()) {
                    std::string param = tokens.at(valIdx);
                    std::string upperParam = param;
                    std::transform(upperParam.begin(), upperParam.end(), upperParam.begin(),
                                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

                    if (upperParam.rfind("TC1=", 0) == 0) {
                        std::string rhs = param.substr(4);
                        if (!rhs.empty()) {
                            // Handle inline TC1=3.9e-3 or TC1={expr}
                            std::vector<std::string> subTokens = { rhs };
                            for (size_t k = valIdx + 1; k < tokens.size() && rhs.find('{') != std::string::npos && rhs.find('}') == std::string::npos; ++k) {
                                subTokens.push_back(tokens.at(k));
                                valIdx = k;
                            }
                            size_t subIdx = 0;
                            double tc1Val = evaluateExpression(extractBracedExpr(subTokens, subIdx));
                            resistor->setTc1(tc1Val);
                        }
                        valIdx++;
                    } else {
                        valIdx++;
                    }
                }

                netlist.addComponent(resistor);
            }
            // =================================================================
            // VOLTAGE SOURCE PARSING:
            //   Vname N1 N2 [DC] Val [AC Mag Phase] [PULSE(V1 V2 TD TR TF PW PER)]
            // =================================================================
            else if (compTypeIndicator == 'v' && tokens.size() >= 4) {
                std::vector<NetNode> nodes = parseNodes(tokens, 2, netlist);
                double dcVal = 0.0;
                double acMag = 0.0;
                double acPhase = 0.0;
                bool acFound = false;

                bool pulseFound = false;
                double pulseV1 = 0.0, pulseV2 = 5.0;
                double pulseTd = 0.0, pulseTr = 1e-12, pulseTf = 1e-12;
                double pulsePw = 10e-6, pulsePer = 20e-6;

                size_t cursor = 3;
                std::string tok3Lower = tokens.at(cursor);
                std::transform(tok3Lower.begin(), tok3Lower.end(), tok3Lower.begin(),
                               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

                if (tok3Lower == "dc") {
                    ++cursor;
                    if (cursor < tokens.size()) {
                        dcVal = evaluateExpression(extractBracedExpr(tokens, cursor));
                    }
                } else if (tok3Lower != "ac" && tok3Lower.rfind("pulse", 0) != 0) {
                    dcVal = evaluateExpression(extractBracedExpr(tokens, cursor));
                }

                // Scan remaining tokens for AC and/or PULSE(...) specifications
                while (cursor < tokens.size()) {
                    std::string kw = tokens.at(cursor);
                    std::transform(kw.begin(), kw.end(), kw.begin(),
                                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

                    if (kw == "ac") {
                        acFound = true;
                        ++cursor;
                        acMag = (cursor < tokens.size())
                                    ? evaluateExpression(extractBracedExpr(tokens, cursor))
                                    : 1.0;
                        acPhase = (cursor < tokens.size())
                                      ? evaluateExpression(extractBracedExpr(tokens, cursor))
                                      : 0.0;
                    } else if (kw.rfind("pulse", 0) == 0) {
                        pulseFound = true;

                        // Collect all tokens belonging to PULSE(...)
                        std::string pulseBlock = tokens.at(cursor);
                        ++cursor;
                        int parenDepth = 0;
                        for (char c : pulseBlock) {
                            if (c == '(') ++parenDepth;
                            else if (c == ')') --parenDepth;
                        }
                        while (parenDepth > 0 && cursor < tokens.size()) {
                            pulseBlock += " " + tokens.at(cursor);
                            for (char c : tokens.at(cursor)) {
                                if (c == '(') ++parenDepth;
                                else if (c == ')') --parenDepth;
                            }
                            ++cursor;
                        }

                        // Extract interior inside PULSE( ... )
                        size_t lParen = pulseBlock.find('(');
                        size_t rParen = pulseBlock.rfind(')');
                        std::string inner = (lParen != std::string::npos && rParen != std::string::npos && rParen > lParen)
                                                ? pulseBlock.substr(lParen + 1, rParen - lParen - 1)
                                                : "";

                        // Replace commas with spaces and split braced expressions cleanly
                        std::replace(inner.begin(), inner.end(), ',', ' ');
                        std::vector<std::string> rawPulseTokens;
                        std::stringstream pss(inner);
                        std::string pTok;
                        while (pss >> pTok) {
                            rawPulseTokens.push_back(pTok);
                        }

                        std::vector<double> pVals;
                        size_t pIdx = 0;
                        while (pIdx < rawPulseTokens.size()) {
                            pVals.push_back(evaluateExpression(extractBracedExpr(rawPulseTokens, pIdx)));
                        }

                        if (pVals.size() >= 1) pulseV1  = pVals[0];
                        if (pVals.size() >= 2) pulseV2  = pVals[1];
                        if (pVals.size() >= 3) pulseTd  = pVals[2];
                        if (pVals.size() >= 4) pulseTr  = std::max(1e-15, pVals[3]);
                        if (pVals.size() >= 5) pulseTf  = std::max(1e-15, pVals[4]);
                        if (pVals.size() >= 6) pulsePw  = pVals[5];
                        if (pVals.size() >= 7) pulsePer = std::max(1e-15, pVals[6]);

                        dcVal = pulseV1; // Initial t=0 DC operating point matches V1
                    } else {
                        ++cursor;
                    }
                }

                auto vSource = std::make_shared<VoltageSource>(compName, nodes, dcVal);
                if (acFound) {
                    vSource->setAcParameters(acMag, acPhase);
                }
                if (pulseFound) {
                    vSource->setPulseParameters(pulseV1, pulseV2, pulseTd, pulseTr, pulseTf, pulsePw, pulsePer);
                }
                netlist.addComponent(vSource);
                
            }
            // =================================================================
            // CAPACITOR PARSING: Cname N1 N2 Value_or_{Expr}
            // =================================================================
            else if (compTypeIndicator == 'c' && tokens.size() >= 4) {
                std::vector<NetNode> nodes = parseNodes(tokens, 2, netlist);
                size_t valIdx = 3;
                double val = evaluateExpression(extractBracedExpr(tokens, valIdx));
                netlist.addComponent(std::make_shared<Capacitor>(compName, nodes, val));
            }
            // =================================================================
            // INDUCTOR PARSING: Lname N1 N2 Value_or_{Expr}
            // =================================================================
            else if (compTypeIndicator == 'l' && tokens.size() >= 4) {
                std::vector<NetNode> nodes = parseNodes(tokens, 2, netlist);
                size_t valIdx = 3;
                double val = evaluateExpression(extractBracedExpr(tokens, valIdx));
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
            //   2. BSIM-CMG Multi-Gate:    Mname D G S [B] FINFET|NFINFET|PFINFET|NRIBBON|NGAA|ModelName [L=... NFIN=...]
            //   3. Standard Level-1 FET:   Mname D G S [B] NMOS|PMOS|ModelName
            // =================================================================
            else if (compTypeIndicator == 'm' && tokens.size() >= 3) {
                auto toUpper = [](std::string s) {
                    std::transform(s.begin(), s.end(), s.begin(),
                                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
                    return s;
                };

                // Step A: Dynamically detect if this 'M' card is a Behavioral Logic Gate
                size_t gateTokenIdx = 0;
                for (size_t i = 2; i < tokens.size(); ++i) {
                    std::string t = toUpper(tokens.at(i));
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

                    double vThresh = (tokens.size() > gateTokenIdx + 1) ? evaluateExpression(tokens.at(gateTokenIdx + 1)) : 2.5;
                    double vHigh   = (tokens.size() > gateTokenIdx + 2) ? evaluateExpression(tokens.at(gateTokenIdx + 2)) : 5.0;
                    double vLow    = (tokens.size() > gateTokenIdx + 3) ? evaluateExpression(tokens.at(gateTokenIdx + 3)) : 0.0;
                    double rOut    = (tokens.size() > gateTokenIdx + 4) ? evaluateExpression(tokens.at(gateTokenIdx + 4)) : 10.0;

                    netlist.addComponent(std::make_shared<BehaviorLogicGateComponent>(
                        compName, nodes, gateType, vThresh, vHigh, vLow, rOut
                    ));
                }
                else if (tokens.size() >= 4) {
                    // Step B: Physical Transistor (3-node {D, G, S} or 4-node {D, G, S, B})
                    auto isKnownModelToken = [&](const std::string& tok) {
                        return tok == "NMOS" || tok == "PMOS" ||
                               tok.find("FINFET") != std::string::npos ||
                               tok.find("RIBBON") != std::string::npos ||
                               tok.find("GAA")    != std::string::npos ||
                               modelRegistry.find(tok) != modelRegistry.end();
                    };

                    std::string modelTokenUpper = "NMOS";
                    size_t modelIdx = 4;

                    // 4 terminals {D, G, S, B} if tok4 is a Bulk node (not a model) and tok5 is the model name
                    if (tokens.size() >= 6 && tokens.at(4).find('=') == std::string::npos) {
                        std::string tok4 = toUpper(tokens.at(4));
                        std::string tok5 = toUpper(tokens.at(5));
                        if (!isKnownModelToken(tok4) && tok5.find('=') == std::string::npos) {
                            modelIdx = 5;
                        }
                    }

                    // Parse 4 nodes {D, G, S, B} if Bulk is present, otherwise 3 nodes {D, G, S}
                    size_t terminalCount = (modelIdx == 5) ? 4 : 3;
                    std::vector<NetNode> nodes = parseNodes(tokens, terminalCount, netlist);

                    if (tokens.size() > modelIdx && tokens.at(modelIdx).find('=') == std::string::npos) {
                        modelTokenUpper = toUpper(tokens.at(modelIdx));
                    }

                    // Resolve .MODEL card metadata if registered
                    auto modelIt = modelRegistry.find(modelTokenUpper);
                    std::string resolvedModelType = (modelIt != modelRegistry.end())
                                                        ? modelIt->second.modelType
                                                        : modelTokenUpper;

                    bool isMultiGateKeyword =
                        (modelTokenUpper == "FINFET"     ||
                         modelTokenUpper == "NFINFET"    || modelTokenUpper == "PFINFET"    ||
                         modelTokenUpper == "NRIBBON"    || modelTokenUpper == "PRIBBON"    ||
                         modelTokenUpper == "NRIBBONFET" || modelTokenUpper == "PRIBBONFET" ||
                         modelTokenUpper == "NGAA"       || modelTokenUpper == "PGAA"       ||
                         modelTokenUpper == "NGAAFET"    || modelTokenUpper == "PGAAFET"    ||
                         resolvedModelType == "BSIMCMG"  || resolvedModelType == "BSIMCMG_VA" ||
                         resolvedModelType == "NFINFET"  || resolvedModelType == "PFINFET");

                    // Dispatch BSIM-CMG Multi-Gate FET (or any .MODEL card with L/W/VTH0 overrides)
                    if (isMultiGateKeyword || modelIt != modelRegistry.end()) {
                        BsimCmgParameters bsimParams;
                        bsimParams.isNChannel = (resolvedModelType != "PMOS" &&
                                                 resolvedModelType != "PFINFET" &&
                                                 modelTokenUpper.front() != 'P');

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

                        // 1. Inherit .MODEL card parameters first (if present)
                        if (modelIt != modelRegistry.end()) {
                            applyModelCardToBsim(modelTokenUpper, bsimParams);
                        }

                        // 2. Apply inline instance overrides (e.g., TYPE=PMOS L={L_MIN} NFIN=4)
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
	    // =================================================================
            // BEHAVIORAL SOURCE PARSING:
            //   Bname N+ N- V=<expression>
            //   Bname N+ N- I=<expression>
            // =================================================================
            else if (compTypeIndicator == 'b' && tokens.size() >= 4) {
                std::vector<NetNode> nodes = parseNodes(tokens, 2, netlist);

                // Concatenate all tokens starting from index 3 to capture spaced formulas
                std::string rhsSpec;
                for (size_t k = 3; k < tokens.size(); ++k) {
                    if (!rhsSpec.empty()) rhsSpec += " ";
                    rhsSpec += tokens.at(k);
                }

                size_t eqPos = rhsSpec.find('=');
                if (eqPos == std::string::npos || eqPos == 0 || eqPos + 1 >= rhsSpec.size()) {
                    throw std::runtime_error("Malformed Behavioral B-source syntax (expected V=<expr> or I=<expr>).");
                }

                std::string modeStr = rhsSpec.substr(0, eqPos);
                modeStr.erase(std::remove_if(modeStr.begin(), modeStr.end(), ::isspace), modeStr.end());
                std::transform(modeStr.begin(), modeStr.end(), modeStr.begin(),
                               [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

                BehavioralSource::Mode bMode = BehavioralSource::Mode::CURRENT;
                if (modeStr == "V") {
                    bMode = BehavioralSource::Mode::VOLTAGE;
                } else if (modeStr == "I") {
                    bMode = BehavioralSource::Mode::CURRENT;
                } else {
                    throw std::runtime_error("Unsupported Behavioral B-source mode '" + modeStr + "' (use V= or I=).");
                }

                std::string exprBody = rhsSpec.substr(eqPos + 1);

                // Bind node name resolver so V(nodeName) inside expression maps to NetNode
                auto resolver = [this, &netlist](const std::string& rawNode) -> NetNode {
                    return this->getOrAssignNode(rawNode, netlist);
                };

                netlist.addComponent(std::make_shared<BehavioralSource>(
                    compName, nodes, bMode, exprBody, paramRegistry, resolver
                ));
            }

        }
        catch (const std::exception& e) {
            std::cerr << "Parser Warning: Skipping malformed line entry element ["
                      << compName << "]: " << e.what() << std::endl;
        }
    }

    // Propagate global .TEMP to BSIM-CMG OSDI devices once after Pass 3 finishes
    double globalTempC = netlist.getTemperatureCelsius();
    netlist.inspectComponents([globalTempC](Component& comp, MnaIndex) {
        if (auto* finfet = dynamic_cast<BsimCmgFet*>(&comp)) {
            finfet->setAmbientTemperature(globalTempC);
        }
    });
}
