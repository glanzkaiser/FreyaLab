#include "freyalab/CircuitComponents.h"
#include "freyalab/math.h"
#include "freyalab/ACCircuit.h"
#include "freyalab/FileParser.h"
#include <iostream>
#include <vector>
#include <memory>
#include <string>

int main() {
    std::string netlistFilename = "netlist.txt";
    std::vector<std::shared_ptr<Component>> components;
    int numNodes = 0;

    std::cout << "[Freyalab Simulator Core v2.0 - AC/Transient Engine]\n";
    std::cout << "Loading circuit topology matrix from: " << netlistFilename << "...\n";

    // 1. Invoke the sanitized FileParser parsing pipeline
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

    // 3. Construct the ACCircuit simulation engine matrix
    std::cout << "Allocating complex system topology network equations...\n";
    ACCircuit circuit(components, numNodes);

    // 4. Verify system operation mode boundaries
    if (!circuit.isAcCircuit()) {
        std::cout << "[Warning] No independent dynamic alternating current (AC) sources detected.\n"
                  << "Engine will execute sweeps using default parameters.\n\n";
    }

    // 5. Run AC Transient Waveform Analysis (Time Domain)
    // Simulates an execution profile (e.g., from t = 0 to 2ms with 5us steps)
    double tEnd = 2e-3;
    double dt = 5e-6;
    std::cout << "Executing time-domain AC Transient Waveform Analysis...\n";
    std::cout << "Simulation Parameters: End Time = " << tEnd << "s, Step Size (dt) = " << dt << "s\n";
    circuit.runTransientWaveformAnalysis(tEnd, dt);

    // 6. Run AC Frequency Sweep Analysis (Bode Data Generation)
    // Sweeps from 10Hz to 100kHz logarithmically at 10 points per decade
    double startFreq = 10.0;
    double endFreq = 100000.0;
    int pointsPerDecade = 10;
    std::cout << "\nExecuting logarithmic AC Frequency Sweep Analysis (Bode Plotting Generation)...\n";
    std::cout << "Sweep Range: " << startFreq << " Hz to " << endFreq << " Hz (" << pointsPerDecade << " pts/dec)\n";
    circuit.runFrequencySweepAnalysis(startFreq, endFreq, pointsPerDecade);

    std::cout << "\nAll alternating current simulation operations completed successfully.\n";
    std::cout << "Bode plot visualization matrices dumped out to 'bode_plot_data.csv'.\n";
    return 0;
}
