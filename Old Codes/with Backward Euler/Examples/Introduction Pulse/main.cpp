#include "freyalab/CircuitComponents.h"
#include "freyalab/math.h"
#include "freyalab/PulseCircuit.h"
#include "freyalab/FileParser.h"
#include <iostream>
#include <vector>
#include <memory>
#include <string>

int main() {
    std::string netlistFilename = "netlist.txt";
    std::vector<std::shared_ptr<Component>> components;
    int numNodes = 0;

    std::cout << "[Freyalab Simulator Core v2.0 - Pulse Transient Engine]\n";
    std::cout << "Loading circuit topology matrix from: " << netlistFilename << "...\n";

    // 1. Invoke the clean FileParser line tokenization engine
    if (!FileParser::parseNetlist(netlistFilename, components, numNodes)) {
        std::cerr << "[Fatal Error] Unable to open or parse standard netlist target file: " 
                  << netlistFilename << "\n";
        return 1;
    }

    std::cout << "Successfully parsed netlist schema structural elements.\n";
    std::cout << "Detected Node Array Size: " << numNodes << " (including Ground Reference Node 0).\n";
    std::cout << "Total Active Physical Components Loaded: " << components.size() << "\n\n";

    // 2. Validate circuit integrity layout parameters
    if (components.empty()) {
        std::cerr << "[Fatal Error] Empty hardware description layout. Execution halted.\n";
        return 1;
    }

    // 3. Construct the PulseCircuit simulation engine framework
    std::cout << "Allocating system topology MNA equations and finalizing sparse matrix patterns...\n";
    PulseCircuit circuit(components, numNodes);

    // 4. Verify operating mode boundaries
    if (!circuit.isPulseCircuit()) {
        std::cout << "[Warning] No dynamic PULSE source component flags identified in netlist.\n"
                  << "Engine will run transient sweeps using default step methods.\n\n";
    }

    // 5. Run Pulse Transient Analysis
    // Sweeps a 20 microsecond trace window using precise 50 nanosecond time-steps
    double durationMicroseconds = 20.0;
    double timeStepNanoseconds = 50.0;

    std::cout << "Executing high-velocity Pulse Transient Analysis...\n";
    std::cout << "Parameters: Duration = " << durationMicroseconds << " us, dt = " << timeStepNanoseconds << " ns\n";
    
    circuit.runPulseTransientAnalysis(durationMicroseconds, timeStepNanoseconds);

    std::cout << "\nPulse transient simulation run completed successfully.\n";
    std::cout << "Time-series waveform data arrays written out cleanly to 'pulse_transient_data.csv'.\n";
    
    return 0;
}
