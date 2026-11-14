#include "freyalab/Netlist.h"
#include "freyalab/FileParser.h"
#include "freyalab/DCCircuit.h"
#include "freyalab/ACCircuit.h"
#include "freyalab/PulseCircuit.h"
#include "freyalab/CircuitComponents.h"
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

int main() {
    std::cout << "====================================================================\n";
    std::cout << "             FREYALAB HIERARCHICAL EXTENSION BOOT ENGINE            \n";
    std::cout << "====================================================================\n\n";

    const std::string netlistFilename = "circuit.txt";

    Netlist circuitWorkspace;
    FileParser circuitParser;

    std::cout << "[Parser] Opening external target file: " << netlistFilename << "\n";
    std::cout << "[Parser] Ingesting line data streams and flattening hierarchical subcircuits...\n";

    try {
        circuitParser.parseNetlistFile(netlistFilename, circuitWorkspace);
        std::cout << "[Parser] Successfully flattened macro layers into core primitive objects.\n";
    } catch (const std::exception& e) {
        std::cerr << "\n[Fatal Error] Subcircuit Expansion Failure: " << e.what() << "\n";
        std::cerr << "Ensure that '" << netlistFilename << "' is located in the working directory.\n";
        return -1;
    }

    // Build sparse matrix topology and assign type-safe MnaIndex rows
    std::cout << "[Solver] Executing structural topology layout diagnostics...\n";
    circuitWorkspace.buildSparsityPattern();

    std::cout << "\nChecking structural component parameters allocation rules below:";
    circuitWorkspace.printComponentInputOutputNodes();

    std::cout << "--------------------------------------------------------------------\n";
    std::cout << "Total Active Solution Matrix Size: " << circuitWorkspace.getTotalSystemSize() << "\n";
    std::cout << "Total Resolved Isolated Nodes Count: " << circuitWorkspace.getNumNodes() << "\n";
    std::cout << "====================================================================\n\n";

    // =================================================================
    // 1. RUN DC OPERATING POINT & SWEEP ANALYSIS
    // =================================================================
    std::cout << "[DC ANALYSIS] Initializing steady-state operating point engine...\n";
    DCCircuit dcEngine(circuitWorkspace);

    dcEngine.runOperatingPointAnalysis();

    // Print bounds-checked NetNode potentials and MnaIndex extra branch values
    const std::vector<double>& dcState = circuitWorkspace.getCurrentState();
    std::cout << "[DC TELEMETRY] Resolved NetNode Voltages:\n";
    for (int n = 1; n <= circuitWorkspace.getNumNodes(); ++n) {
        NetNode node{n};
        std::cout << "  V(" << circuitWorkspace.getNodeName(node) << " / Node_" << node.id << ") = "
                  << std::fixed << std::setprecision(6)
                  << getSafeVoltage(dcState, node) << " V\n";
    }

    circuitWorkspace.inspectComponents([&dcState](const Component& comp, MnaIndex extraIdx) {
        if (extraIdx.isValid()) {
            std::cout << "  Extra MNA Row [" << extraIdx.idx << "] (" << comp.getName() << ") = "
                      << std::scientific << std::setprecision(4)
                      << getSafeExtraValue(dcState, extraIdx) << "\n";
        }
    });

    try {
        dcEngine.runDcSweepAnalysis("V_SOURCE", 0.0, 10.0, 0.5);
    } catch (const std::exception& e) {
        std::cerr << "[DC Sweep Warning] " << e.what() << "\n";
    }

    // =================================================================
    // 2. RUN AC SMALL-SIGNAL FREQUENCY DOMAIN ANALYSIS
    // =================================================================
    std::cout << "\n[AC ANALYSIS] Initializing small-signal phasor sweep engine...\n";
    ACCircuit acEngine(circuitWorkspace);

    if (acEngine.isAcCircuit()) {
        std::cout << "[AC Engine] AC source detected. Initiating logarithmic frequency sweep...\n";
        acEngine.runFrequencySweepAnalysis(10.0, 10.0e6, 100);
    } else {
        std::cout << "[AC Engine] Skip: No dynamic AC magnitude flags found on circuit sources.\n";
    }

    // =================================================================
    // 3. RUN PULSE HIGH-SPEED TRANSIENT ANALYSIS (TR-BDF2)
    // =================================================================
    std::cout << "\n[PULSE TRANSIENT] Initializing two-stage composite integration engine...\n";
    PulseCircuit pulseEngine(circuitWorkspace);

    if (pulseEngine.isPulseCircuit()) {
        PulseSimulationOptions pulseOptions;
        pulseOptions.startTime = 0.0;
        pulseOptions.stopTime  = 100e-6;   // 100 microseconds total simulation horizon
        pulseOptions.initialDt = 10e-9;    // 10 nanoseconds starting step size
        pulseOptions.minDt     = 1e-12;    // 1 picosecond adaptive timestep floor boundary
        pulseOptions.maxDt     = 1e-6;     // 1 microsecond timestep ceiling capping

        // Automated scanning gathers Cause/Effect NetNodes and binds measurements relative to Ground (NetNode{0})
        std::vector<std::pair<NetNode, NetNode>> monitorPairs;
        pulseEngine.autoGenerateMonitorPairs(monitorPairs);

        pulseEngine.runPulseTransientAnalysis(pulseOptions, monitorPairs);
    } else {
        std::cout << "[Pulse Engine] Skip: No active non-linear switching devices discovered.\n";
    }

    std::cout << "====================================================================\n";
    std::cout << "SIMULATION SUCCESSFUL: All waveform output CSV streams locked onto disk.\n";
    std::cout << "====================================================================\n";

    return 0;
}
