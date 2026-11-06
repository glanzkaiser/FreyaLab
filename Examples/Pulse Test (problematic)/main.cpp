#include "freyalab/FileParser.h"
#include "freyalab/Netlist.h"
#include "freyalab/PulseCircuit.h"
#include <iostream>
#include <memory>
#include <stdexcept>

int main() {
    std::cout << "========================================================================\n"
              << "              FREYALAB PULSE TRANSIENT SIMULATION TESTBENCH             \n"
              << "========================================================================\n" << std::endl;

    // 1. Instantiate the core polymorphic containers
    Netlist circuitNetlist;
    FileParser netlistParser;
    std::string netlistFilename = "netlist.txt";

    try {
        std::cout << "Step 1: Ingesting circuit architecture file [ " << netlistFilename << " ]..." << std::endl;
        netlistParser.parseNetlistFile(netlistFilename, circuitNetlist);
        
        std::cout << "  -> Netlist elements compiled successfully." << std::endl;
        std::cout << "  -> Total raw component nodes found: " << circuitNetlist.getTotalSystemSize()  << std::endl;

        // 2. Map topological layout using Option B (Permanently bind matrix counters)
        std::cout << "\nStep 2: Securing internal sparse matrix topology footprint..." << std::endl;
        circuitNetlist.buildSparsityPattern();

        int systemSize = circuitNetlist.getTotalSystemSize();
        if (systemSize == 0) {
            throw std::runtime_error("Empty Structural Exception: Sizing index calculation failed.");
        }

        std::cout << "  -> Linear System Equation Dimensions: " << systemSize << "x" << systemSize << std::endl;

        // 3. Structural Pulse Validation Checks
        PulseCircuit pulseEngine(circuitNetlist);
        std::cout << "\nStep 3: Checking target circuit validation state..." << std::endl;
        
        if (!pulseEngine.isPulseCircuit()) {
            std::cerr << "  [ABORT] System matrix state failed basic size validation." << std::endl;
            return EXIT_FAILURE;
        }
        std::cout << "  -> Validation Succeeded: Active structural network verified." << std::endl;

        // 4. Configure Adaptive TR-BDF2 Integration Options
        std::cout << "\nStep 4: Designing simulation execution boundaries..." << std::endl;
        
        PulseSimulationOptions pulseOptions;
        pulseOptions.startTime     = 0.0;       // t = 0 us
        pulseOptions.stopTime      = 50e-6;     // Total simulation timeline sweep (50 us)
        pulseOptions.initialDt     = 1e-9;      // Baseline initial delta-time window (1 ns)
        pulseOptions.minDt         = 1e-15;     // Absolute timestep truncation floor (1 fg)
        pulseOptions.maxDt         = 5e-7;      // Maximum timestep upper cap boundaries (500 ns)
        pulseOptions.maxIterations = 80;        // Newton-Raphson iteration ceiling cap per substep
        pulseOptions.relTol        = 1e-3;      // Relative integration error parameter
        pulseOptions.absTol        = 1e-6;      // Absolute integration error parameter

        // 5. Execute 2nd-Order L-Stable TR-BDF2 Transient Pulse Analysis
        // Trace Input node 1 (Base drive pulse line) vs Output node 3 (BJT collector node)
        int inputDriverNode = 1;
        int outputProbedNode = 3;

        std::cout << "\nStep 5: Executing transient solver pipeline..." << std::endl;
        pulseEngine.runPulseTransientAnalysis(pulseOptions, inputDriverNode, outputProbedNode);

        std::cout << "\n========================================================================\n"
                  << "   SUCCESS: PULSE TRANSIENT WAVEFORM ACQUISITION IN LOOP TERMINATED     \n"
                  << "========================================================================" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "\n[CRITICAL ENGINE EXCEPTION] Transient execution aborted: " 
                  << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}