#include "freyalab/Netlist.h"
#include "freyalab/FileParser.h"
#include "freyalab/DCCircuit.h"
#include "freyalab/ACCircuit.h"
#include "freyalab/PulseCircuit.h"
#include <iostream>
#include <memory>
#include <string>

int main() {
    std::cout << "====================================================================\n";
    std::cout << "             FREYALAB HIERARCHICAL EXTENSION BOOT ENGINE            \n";
    std::cout << "====================================================================\n\n";

    // Target the external independent structural netlist file
    const std::string netlistFilename = "circuit.txt";
    
    // Initialize clean solver containers and engine parser instances
    Netlist circuitWorkspace;
    FileParser circuitParser;

    std::cout << "[Parser] Opening external target vector: " << netlistFilename << "\n";
    std::cout << "[Parser] Ingesting line data streams and processing macro matrices...\n";
    
    try {
        // Run the multi-pass hierarchical macro compilation pipeline directly on the file
        circuitParser.parseNetlistFile(netlistFilename, circuitWorkspace);
        std::cout << "[Parser] Successfully flattened macro layers into core primitive objects.\n";
    } catch (const std::exception& e) {
        std::cerr << "\n[Fatal Error] Subcircuit Expansion Failure: " << e.what() << "\n";
        std::cerr << "Ensure that '" << netlistFilename << "' is located in the working directory.\n";
        return -1;
    }

    // Dynamically construct sparse patterns to sync unshifted 0-based layout offsets
    std::cout << "[Solver] Executing structural topology layout diagnostics...\n";
    circuitWorkspace.buildSparsityPattern();

    // Execute the input-output terminal printing tracking metrics 
    std::cout << "\nChecking structural component parameters allocation rules below:";
    circuitWorkspace.printComponentInputOutputNodes();

    std::cout << "--------------------------------------------------------------------\n";
    std::cout << "Total Active Solution Matrix Size: " << circuitWorkspace.getTotalSystemSize() << "\n";
    std::cout << "Total Resolved Isolated Nodes Count: " << circuitWorkspace.getNumNodes() << "\n";
    std::cout << "====================================================================\n";

    std::cout << "[Parser] Loading and flattening hierarchical subcircuits...\n";
    try {
        circuitParser.parseNetlistFile(netlistFilename, circuitWorkspace);
        std::cout << "[Parser] Hierarchy successfully flattened into base primitives.\n\n";
    } catch (const std::exception& e) {
        std::cerr << "[Fatal Error] Subcircuit expansion failure: " << e.what() << "\n";
        return -1;
    }

    // =================================================================
    // 1. RUN DC OPERATING POINT & SWEEP ANALYSIS
    // =================================================================
    std::cout << "[DC ANALYSIS] Initializing steady-state operating point engine...\n";
    DCCircuit dcEngine(circuitWorkspace);
    
    // Solves the background Newton-Raphson system and outputs a formatted terminal matrix table
    dcEngine.runOperatingPointAnalysis();

    // Executes a parametric source sweep on V_SOURCE from 0V to 10V in 0.5V steps
    // Automatically writes data to: dc_sweep_nodes.csv, dc_sweep_branches.csv, dc_sweep_power.csv
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
        // Sweeps frequencies from 10 Hz to 10 MHz across 100 points
        // Automatically captures linearized reactive parameters and outputs 'ac_frequency_sweep.csv'
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
        pulseOptions.stopTime = 100e-6;    // 100 microseconds total simulation horizon
        pulseOptions.initialDt = 10e-9;   // 10 nanoseconds starting step size
        pulseOptions.minDt = 1e-12;       // 1 picosecond adaptive timestep floor boundary
        pulseOptions.maxDt = 1e-6;        // 1 microsecond timestep ceiling capping

        std::vector<std::pair<int, int>> monitorPairs;
        // Automated scanning gathers Cause/Effect branches and binds measurements relative to Ground (0)
        pulseEngine.autoGenerateMonitorPairs(monitorPairs);

        // Runs the non-dissipative TR-BDF2 time integration sequence
        // Automatically manages local truncation error steps and outputs 'pulse_transient.csv'
        pulseEngine.runPulseTransientAnalysis(pulseOptions, monitorPairs);
    } else {
        std::cout << "[Pulse Engine] Skip: No active non-linear switching devices discovered.\n";
    }

    std::cout << "====================================================================\n";
    std::cout << "SIMULATION SUCCESSFUL: All waveform output CSV streams locked onto disk.\n";
    std::cout << "====================================================================\n";

    return 0;
}
