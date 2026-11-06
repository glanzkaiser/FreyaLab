
#include "freyalab/FileParser.h"
#include "freyalab/DCCircuit.h"
#include "freyalab/ACCircuit.h"
#include "freyalab/PulseCircuit.h"
#include <iostream>
#include <vector>
#include <memory>


int main() {
    std::string netlistFile = "netlist.txt";
    std::vector<std::shared_ptr<Component>> components;
    int numNodes = 0;

	 // Create an instance of the parser class
    FileParser parser; 
    std::cout << "[Solver Engine] Initializing CSR Sparse Netlist Parser..." << std::endl;
    if (!parser.parseNetlist(netlistFile, components, numNodes)) {
        std::cerr << "[Error] Failed to read or parse: " << netlistFile << std::endl;
        return 1;
    }

    DCCircuit dcAnalyst(components, numNodes);
    ACCircuit acAnalyst(components, numNodes);
    PulseCircuit pulseAnalyst(components, numNodes);

    if (pulseAnalyst.isPulseCircuit()) {
        std::cout << "[Status] PULSE Netlist configuration signature verified." << std::endl;
        pulseAnalyst.runPulseTransientAnalysis(25.0, 50.0);
    } 
    else if (acAnalyst.isAcCircuit()) {
        std::cout << "[Status] AC Netlist configuration signature verified." << std::endl;
        acAnalyst.runTransientWaveformAnalysis(2e-3, 5e-6);
        acAnalyst.runFrequencySweepAnalysis(10.0, 100000.0, 10);
    } 
    else if (dcAnalyst.isDcCircuit()) {
        std::cout << "[Status] DC Netlist configuration signature verified." << std::endl;
        dcAnalyst.runOperatingPointAnalysis();
        dcAnalyst.runDcSweepAnalysis(0.0, 12.0, 0.5);
    }

    return 0;
}
