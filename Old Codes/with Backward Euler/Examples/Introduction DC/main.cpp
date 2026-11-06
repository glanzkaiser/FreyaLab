#include "freyalab/CircuitComponents.h"
#include "freyalab/math.h"
#include "freyalab/DCCircuit.h"
#include "freyalab/FileParser.h"
#include <iostream>
#include <vector>
#include <memory>
#include <string>

int main() {
    std::string netlistFilename = "netlist.txt";
    std::vector<std::shared_ptr<Component>> components;
    int numNodes = 0;

    std::cout << "[Freyalab Simulator Core v2.0]\n";
    std::cout << "Loading circuit topology matrix from: " << netlistFilename << "...\n";

    // 1. Invoke the sanitized FileParser pipeline to build polymorphic pointers
    if (!FileParser::parseNetlist(netlistFilename, components, numNodes)) {
        std::cerr << "[Fatal Error] Unable to open or parse standard netlist target file: " 
                  << netlistFilename << "\n";
        return 1;
    }

    std::cout << "Successfully parsed netlist schema structural elements.\n";
    std::cout << "Detected Node Array Size: " << numNodes << " (including Ground Reference Node 0).\n";
    std::cout << "Total Active Physical Components Loaded: " << components.size() << "\n\n";

    // 2. Validate circuit integrity topology bounds
    if (components.empty()) {
        std::cerr << "[Fatal Error] Empty hardware description layout. Execution halted.\n";
        return 1;
    }

    // 3. Construct the direct MNA operating matrix engine framework
    std::cout << "Allocating network topology equations and finalizing CSR matrix sparsity patterns...\n";
    DCCircuit circuit(components, numNodes);

    // 4. Verify system operation mode boundaries
    if (!circuit.isDcCircuit()) {
        std::cout << "[Warning] Non-static sources encountered. "
                  << "Engine will treat reactive elements (Capacitors/Inductors) under DC operating conditions.\n";
    }

    // 5. Run the high-stability Operating Point Analysis
    std::cout << "Executing Newton-Raphson Modified Nodal Analysis (MNA)...\n";
    circuit.runOperatingPointAnalysis();

    std::cout << "\nSimulation execution task completed successfully.\n";
    return 0;
}
