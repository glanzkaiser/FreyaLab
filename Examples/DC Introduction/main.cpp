#include "freyalab/Netlist.h"
#include "freyalab/DCCircuit.h"
#include "freyalab/ACCircuit.h"
#include "freyalab/PulseCircuit.h"
#include "freyalab/FileParser.h"
#include <iostream>
#include <string>
#include <memory>
#include <algorithm>
#include <iomanip>
void printFullLinearSystem(Netlist& netlist) {
    int size = netlist.getTotalSystemSize();
    const CsrMatrix& G = netlist.getGMatrix();
    const std::vector<double>& B = netlist.getBVector();

    // Reconstruct dense layout from CSR representation vectors for debugging visualization
    std::vector<std::vector<double>> denseG(size, std::vector<double>(size, 0.0));
    for (int r = 0; r < size; ++r) {
        int start = G.rowPointers[r];
        int end = G.rowPointers[r + 1];
        for (int p = start; p < end; ++p) {
            int c = G.colIndices[p];
            denseG[r][c] = G.values[p];
        }
    }

    std::cout << "\n========================================================================================\n";
    std::cout << "                 MNA SYSTEM MATRIX LAYER STATE (INITIAL EXPERIMENTAL STAMP)               \n";
    std::cout << "========================================================================================\n\n";
    
    // Header labels tracker
    std::cout << "      ";
    for (int c = 0; c < size; ++c) {
        std::cout << "  Col " << std::left << std::setw(2) << c << "    ";
    }
    std::cout << "  |   RHS Vector B\n";
    std::cout << std::string(12, '-') << std::string(size * 11, '-') << "---|------------------\n";

    std::cout << std::scientific << std::setprecision(3);
    for (int r = 0; r < size; ++r) {
        std::cout << "Row " << std::left << std::setw(2) << r << ":";
        for (int c = 0; c < size; ++c) {
            std::cout << " " << (denseG[r][c] >= 0 ? " " : "") << denseG[r][c];
        }
        std::cout << "  |  " << (B[r] >= 0 ? " " : "") << B[r] << "\n";
    }
    std::cout << "========================================================================================\n";
}

void printUsage() {
    std::cout << "========================================================================\n"
              << "              SPICE Engine (TR-BDF2 + VBIC Left-Looking LU)             \n"
              << "========================================================================\n"
              << "Usage: ./spice_solver <netlist.txt> <mode> [options]\n\n"
              << "Available Execution Modes:\n"
              << "  --dc-op        : Run static DC Operating Point analysis\n"
              << "  --dc-sweep     : Run DC Parametric Source Sweep analysis\n"
              << "                   Requires: <src_name> <start> <stop> <step>\n"
              << "  --ac-trans     : Run AC Time-Domain Transient Waveform simulation\n"
              << "                   Requires: <input_node> <output_node>\n"
              << "  --ac-sweep     : Run phasor small-signal Frequency AC Sweep (Bode Data)\n"
              << "                   Requires: <start_freq> <stop_freq> <points>\n"
              << "  --pulse-trans  : Run high-speed Pulse Transient switching analysis\n"
              << "                   Requires: <input_node> <output_node>\n"
              << "========================================================================\n";
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        printUsage();
        return 1;
    }

    std::string netlistFile = argv[1];
    std::string runMode = argv[2];
    std::transform(runMode.begin(), runMode.end(), runMode.begin(), ::tolower);

    // 1. Initialize network structures and parse the input target circuit
    Netlist netlist;
    FileParser parser;
    try {
        std::cout << "Parsing netlist schema target file: " << netlistFile << "..." << std::endl;
        parser.parseNetlistFile(netlistFile, netlist);
        
        // FIX: You must invoke Option B's structural allocation pass to secure row indices
        std::cout << "Locking circuit network graph matrix coordinates..." << std::endl;
        netlist.buildSparsityPattern();
        
        std::cout << "Matrix Equation Allocation complete. Total system size: " 
                  << netlist.getTotalSystemSize() << "x" << netlist.getTotalSystemSize() << std::endl;
                  
    } catch (const std::exception& e) {
        std::cerr << "Fatal Initialization Error: " << e.what() << std::endl;
        return 1;
    }

    // 2. Route execution flow paths dynamically to the appropriate analysis sub-block wrapper
    try {
        // ---------------------------------------------------------------------
        // MODE A: DC OPERATING POINT ANALYSIS
        // ---------------------------------------------------------------------
        if (runMode == "--dc-op") {
            // =================================================================
            // EXPERIMENTAL BASELINE STAMPING FOR INITIAL MATRIX DEBUG PEEKING
            // =================================================================
            // 1. Buat status tebakan awal 0V/0A untuk menghitung nilai konduktansi non-linier awal
            std::vector<double> dummyState(netlist.getTotalSystemSize(), 0.0);
            
            // 2. Bungkus matriks G dan vektor B ke dalam konteks MNA
            MnaMatrixContext peekCtx{
                netlist.getGMatrix(),
                netlist.getBVector(),
                &dummyState, 
                &netlist.getHistoricV(), 
                &netlist.getMidpointV(),
                0.0,   // time
                1e-6,  // dt
                0      // orderStep: DC operating point
            };
            
            // 3. Paksa komponen menuliskan nilai konduktansi/kapasitansi awal mereka ke matriks G
            netlist.inspectComponents([&peekCtx](const Component& comp, int extraIdx) {
                // Gunakan const_cast untuk memanggil non-const method .stamp() secara aman
                const_cast<Component&>(comp).stamp(peekCtx);
            });

            // 4. Cetak visualisasi tata letak matriks MNA padat yang sudah terisi data
            printFullLinearSystem(netlist);
            
            // 5. Jalankan solver sirkuit Newton-Raphson yang sesungguhnya hingga konvergen
            DCCircuit dcAnalyzer(netlist);
            dcAnalyzer.runOperatingPointAnalysis();
        }
        // ---------------------------------------------------------------------
        // MODE B: DC PARAMETRIC SOURCE SWEEP ANALYSIS
        // ---------------------------------------------------------------------
        else if (runMode == "--dc-sweep") {
            if (argc < 7) {
                std::cerr << "Error: --dc-sweep requires: <source_name> <start> <stop> <step>\n";
                return 1;
            }
            std::string srcName = argv[3];
            double startVal = std::stod(argv[4]);
            double stopVal  = std::stod(argv[5]);
            double stepVal  = std::stod(argv[6]);

            DCCircuit dcAnalyzer(netlist);
            dcAnalyzer.runDcSweepAnalysis(srcName, startVal, stopVal, stepVal);
        }
        // ---------------------------------------------------------------------
        // MODE C: AC TIME-DOMAIN TRANSIENT SIMULATION
        // ---------------------------------------------------------------------
        else if (runMode == "--ac-trans") {
            if (argc < 5) {
                std::cerr << "Error: --ac-trans requires: <input_node_idx> <output_node_idx>\n";
                return 1;
            }
            int inNode  = std::stoi(argv[3]);
            int outNode = std::stoi(argv[4]);

            ACCircuit acAnalyzer(netlist);
            if (!acAnalyzer.isAcCircuit()) {
                std::cout << "Warning: Engine detected no dynamic states or AC sources.\n";
            }
            
            // FIX: Explicitly set safe parameter options to prevent uninitialized data reads
            AcSimulationOptions transOptions;
            transOptions.startTime = 0.0;
            transOptions.stopTime  = 5e-3; // 5 ms
            transOptions.initialDt = 1e-6; // 1 us
            transOptions.minDt     = 1e-12;
            transOptions.maxDt     = 1e-4;
            transOptions.relTol    = 1e-3;
            transOptions.absTol    = 1e-6;

            acAnalyzer.runTransientWaveformAnalysis(transOptions, inNode, outNode);
        }
        // ---------------------------------------------------------------------
        // MODE D: PHASOR SMALL-SIGNAL AC FREQUENCY SWEEP
        // ---------------------------------------------------------------------
        else if (runMode == "--ac-sweep") {
            if (argc < 6) {
                std::cerr << "Error: --ac-sweep requires: <start_freq> <stop_freq> <points_count>\n";
                return 1;
            }
            double startFreq = std::stod(argv[3]);
            double stopFreq  = std::stod(argv[4]);
            int pointsCount  = std::stoi(argv[5]);

            ACCircuit acAnalyzer(netlist);
            acAnalyzer.runFrequencySweepAnalysis(startFreq, stopFreq, pointsCount);
        }
        // ---------------------------------------------------------------------
        // MODE E: PULSE TRANSIENT SIMULATION
        // ---------------------------------------------------------------------
        else if (runMode == "--pulse-trans") {
            if (argc < 5) {
                std::cerr << "Error: --pulse-trans requires: <input_node_idx> <output_node_idx>\n";
                return 1;
            }
            int inNode  = std::stoi(argv[3]);
            int outNode = std::stoi(argv[4]);

            PulseCircuit pulseAnalyzer(netlist);
            if (!pulseAnalyzer.isPulseCircuit()) {
                std::cout << "Warning: Netlist structural layout contains no active switching metrics.\n";
            }
            
            // FIX: Explicitly set default switching boundary conditions
            PulseSimulationOptions pulseOptions;
            pulseOptions.startTime     = 0.0;
            pulseOptions.stopTime      = 200e-6; // 200 us
            pulseOptions.initialDt     = 1e-9;   // 1 ns
            pulseOptions.minDt         = 1e-15;
            pulseOptions.maxDt         = 1e-6;
            pulseOptions.maxIterations = 100;
            pulseOptions.relTol        = 1e-3;
            pulseOptions.absTol        = 1e-6;

            pulseAnalyzer.runPulseTransientAnalysis(pulseOptions, inNode, outNode);
        }
        else {
            std::cerr << "Error: Unrecognized execution parameter mode flag encountered: " << runMode << "\n";
            printUsage();
            return 1;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "\nSimulation Execution Failure: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "\nExecution task complete. Processing termination cleanups..." << std::endl;
    return 0;
}
