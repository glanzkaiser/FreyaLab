#pragma once

#include "CircuitComponents.h"
#include <string>
#include <vector>
#include <memory>

class FileParser {
private:
    // Parses numeric string tokens with SPICE suffixes (e.g., 10k, 2.5meg, 100n)
    static double parseValue(const std::string& valueStr);

    // Strips out non-standard punctuation delimiters like commas, equals, and parentheses
    static std::string cleanAndTokenize(const std::string& line);

public:
    /**
     * @brief Parses a standard SPICE-like netlist file.
     * 
     * @param filename Path to the source netlist text file.
     * @param outComponents Vector populated with polymorphic smart pointers to parsed components.
     * @param numNodes System matrix size variable, outputs total nodes including ground.
     * @return true if file opened and parsed successfully, false otherwise.
     */
    static bool parseNetlist(const std::string& filename, 
                             std::vector<std::shared_ptr<Component>>& outComponents, 
                             int& numNodes);
};

