#include "freyalab/Netlist.h"
#include "freyalab/DCCircuit.h"
#include "freyalab/ACCircuit.h"
#include "freyalab/PulseCircuit.h"
#include "freyalab/FileParser.h"
#include "freyalab/CircuitComponents.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

void printFullLinearSystem(Netlist& netlist) {
    int size = netlist.getTotalSystemSize();
    const CsrMatrix& G = netlist.getGMatrix();
    const std::vector<double>& B = netlist.getBVector();

    std::vector<std::vector<double>> denseG(size, std::vector<double>(size, 0.0));
    for (int r = 0; r < size; ++r) {
        int start = G.rowPointers[r];
        int end = G.rowPointers[r + 1];
        for (int p = start; p < end; ++p) {
            int c = G.colIndices[p];
            if (c >= 0 && c < size) {
                denseG[r][c] = G.values[p];
            }
        }
    }

    std::cout << "\n========================================================================================\n";
    std::cout << "                 MNA SYSTEM MATRIX LAYER STATE (INITIAL EXPERIMENTAL STAMP)               \n";
    std::cout << "========================================================================================\n\n";

    std::cout << "      ";
    for (int c = 0; c < size; ++c) {
        std::cout << "  Col " << std::left << std::setw(2) << c << "    ";
    }
    std::cout << "  |   RHS Vector B\n";
    std::cout << std::string(12, '-') << std::string(size * 11, '-') << "---|------------------\n";

    std::cout << std::scientific << std::setprecision(3);
    for (int r = 0; r < size; ++r) {
        double rhsVal = getSafeExtraValue(B, MnaIndex{r});
        std::cout << "Row " << std::left << std::setw(2) << r << ":";
        for (int c = 0; c < size; ++c) {
            std::cout << " " << (denseG[r][c] >= 0 ? " " : "") << denseG[r][c];
        }
        std::cout << "  |  " << (rhsVal >= 0 ? " " : "") << rhsVal << "\n";
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
              << "                   Optional: <input_node> <output_node>\n"
              << "========================================================================\n";
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        printUsage();
        return 1;
    }

    // KUNCI PERBAIKAN: Menggunakan indeks eksplisit [1] dan [2] untuk argumen dasar
    std::string netlistFile = argv[1];
    std::string runMode = argv[2];
    std::transform(runMode.begin(), runMode.end(), runMode.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    Netlist netlist;
    FileParser parser;
    try {
        std::cout << "Parsing netlist schema target file: " << netlistFile << "..." << std::endl;
        parser.parseNetlistFile(netlistFile, netlist);

        std::cout << "Locking circuit network graph matrix coordinates..." << std::endl;
        netlist.buildSparsityPattern();

        std::cout << "Matrix Equation Allocation complete. Total system size: "
                  << netlist.getTotalSystemSize() << "x" << netlist.getTotalSystemSize() << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Fatal Initialization Error: " << e.what() << std::endl;
        return 1;
    }

    try {
        // ---------------------------------------------------------------------
        // MODE A: DC OPERATING POINT ANALYSIS
        // ---------------------------------------------------------------------
        if (runMode == "--dc-op") {
            std::vector<double> dummyState(netlist.getTotalSystemSize(), 0.0);

            MnaMatrixContext peekCtx{
                netlist.getGMatrix(),
                netlist.getBVector(),
                &dummyState,
                &netlist.getHistoricV(),
                &netlist.getMidpointV(),
                0.0,        // time
                1e-6,       // dt
                0,          // orderStep: DC operating point
                MnaIndex{}
            };

            netlist.inspectComponents([&peekCtx](Component& comp, MnaIndex extraIdx) {
                peekCtx.assignedExtraIdx = extraIdx;
                comp.stamp(peekCtx);
            });

            printFullLinearSystem(netlist);

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
            // KUNCI PERBAIKAN: Mengindeks array posisi 3, 4, 5, dan 6 secara berurutan
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
            // KUNCI PERBAIKAN: Mengindeks array posisi 3 dan 4 untuk nomor node
            NetNode inNode{ std::stoi(argv[3]) };
            NetNode outNode{ std::stoi(argv[4]) };

            ACCircuit acAnalyzer(netlist);
            if (!acAnalyzer.isAcCircuit()) {
                std::cout << "Warning: Engine detected no dynamic states or AC sources.\n";
            }

            AcSimulationOptions transOptions;
            transOptions.startTime = 0.0;
            transOptions.stopTime  = 5e-3;  // 5 ms
            transOptions.initialDt = 1e-6;  // 1 us
            transOptions.minDt     = 1e-12; // 1 ps
            transOptions.maxDt     = 1e-4;  // 100 us

            acAnalyzer.runTransientWaveformAnalysis(transOptions, inNode, outNode);
        }
        // ---------------------------------------------------------------------
        // MODE D: SMALL-SIGNAL AC FREQUENCY SWEEP (BODE PLOT)
        // ---------------------------------------------------------------------
        else if (runMode == "--ac-sweep") {
            if (argc < 6) {
                std::cerr << "Error: --ac-sweep requires: <start_freq> <stop_freq> <points>\n";
                return 1;
            }
            // KUNCI PERBAIKAN: Mengindeks array posisi 3, 4, dan 5 untuk variabel reaktif
            double startFreq = std::stod(argv[3]);
            double stopFreq  = std::stod(argv[4]);
            int pointsCount  = std::stoi(argv[5]);

            ACCircuit acAnalyzer(netlist);
            acAnalyzer.runFrequencySweepAnalysis(startFreq, stopFreq, pointsCount);
        }
        // ---------------------------------------------------------------------
        // MODE E: PULSE TRANSIENT SWITCHING ANALYSIS
        // ---------------------------------------------------------------------
        else if (runMode == "--pulse-trans") {
            PulseCircuit pulseAnalyzer(netlist);
            if (!pulseAnalyzer.isPulseCircuit()) {
                std::cout << "Warning: Engine detected no active nonlinear switching devices.\n";
            }

            PulseSimulationOptions pulseOpts;
            std::vector<std::pair<NetNode, NetNode>> monitorPairs;

            if (argc >= 5) {
                // KUNCI PERBAIKAN: Mengindeks array posisi 3 dan 4 untuk pemantau pulsa transien
                NetNode inNode{ std::stoi(argv[3]) };
                NetNode outNode{ std::stoi(argv[4]) };
                monitorPairs.push_back({inNode, NetNode{0}});
                monitorPairs.push_back({outNode, NetNode{0}});
            } else {
                pulseAnalyzer.autoGenerateMonitorPairs(monitorPairs);
            }

            pulseAnalyzer.runPulseTransientAnalysis(pulseOpts, monitorPairs);
        }
        else {
            std::cerr << "Error: Unknown execution mode '" << runMode << "'\n";
            printUsage();
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Simulation Execution Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
