#include "freyalab/ACCircuit.h"
#include "freyalab/math.h"
#include "freyalab/CircuitComponents.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <algorithm>
#include <stdexcept>

ACCircuit::ACCircuit(Netlist& netlistInstance) : netlist(netlistInstance) {}

bool ACCircuit::isAcCircuit() {
    return netlist.getTotalSystemSize() > 0;
}

bool ACCircuit::solveTransientStep(MnaMatrixContext& ctx, double targetTime, double currentDt, int stepPhase) {
    ctx.time = targetTime;
    ctx.dt = currentDt;
    ctx.orderStep = stepPhase;
    return netlist.runNewtonRaphson(ctx, 50, 1e-3, 1e-6);
}

void ACCircuit::runTransientWaveformAnalysis(const AcSimulationOptions& options, int inputNode, int outputNode) {
    std::cout << "\n====================================================================================\n"
              << "                        AC TRANSIENT WAVEFORM TABLE (2ND-ORDER TR-BDF2)            \n"
              << "====================================================================================\n"
              << std::left << std::setw(12) << "Time(s)" << std::setw(15) << "V(input)" << std::setw(15) << "V(output)" << "\n"
              << "------------------------------------------------------------------------------------\n";

    if (netlist.getTotalSystemSize() == 0) {
        netlist.buildSparsityPattern();
    }

    double t = options.startTime;
    double dt = options.initialDt;
    double stopT = options.stopTime;
    
    double gamma = 2.0 - std::sqrt(2.0); // Optimal L-stable TR-BDF2 configuration parameter

    int currentExtraCount = netlist.getTotalSystemSize() - netlist.getNumNodes();
    
    // Allocate context with direct vector references to prevent internal parsing de-synchronization
    MnaMatrixContext ctx{ 
        const_cast<CsrMatrix&>(netlist.getGMatrix()), 
        const_cast<std::vector<double>&>(netlist.getHistoricV()), 
        &netlist.getStateV(), 
        &netlist.getHistoricV(), 
        &netlist.getMidpointV(), 
        t, dt, 0, currentExtraCount 
    };

    try {
        // ====================================================================
        // STEP 1: INITIAL DC OPERATING POINT ANALYSIS PHASE
        // ====================================================================
        ctx.orderStep = 0;
        ctx.time = t;
        ctx.dt = 0.0;
        ctx.currentIterV = &netlist.getStateV();
        ctx.historicStepV = &netlist.getHistoricV();
        ctx.midpointStepV = &netlist.getMidpointV();

        if (!netlist.runNewtonRaphson(ctx, 100, options.relTol, options.absTol)) {
            throw std::runtime_error("AC Transient baseline Operating Point diverged.");
        }
        
        // Synchronize operational baseline states
        netlist.getHistoricV() = netlist.getStateV();
        
        // Critical: Force components to calculate baseline history (q_hist_n / i_hist_n) at t=0
        netlist.inspectComponents([&ctx](const Component& comp, int extraIdx) {
            (void)extraIdx;
            const_cast<Component&>(comp).updateHistory(ctx);
        });

        // Print true baseline t=0 starting coordinates
        const auto& initialState = netlist.getCurrentState();
        double vIn = (inputNode > 0 && static_cast<size_t>(inputNode - 1) < initialState.size()) ? initialState[inputNode - 1] : 0.0;
        double vOut = (outputNode > 0 && static_cast<size_t>(outputNode - 1) < initialState.size()) ? initialState[outputNode - 1] : 0.0;
        std::cout << std::left << std::setw(12) << t << std::setw(15) << vIn << std::setw(15) << vOut << "\n";

        // ====================================================================
        // STEP 2: TIME STEPPING SIMULATION LOOP
        // ====================================================================
        while (t < stopT) {
            if (t + dt > stopT) {
                dt = stopT - t;
            }

            std::vector<double> backupState = netlist.getStateV();

            // Freeze charge layers across dynamic primitives before running integration steps
            netlist.inspectComponents([](const Component& comp, int extraIdx) {
                (void)extraIdx;
                const_cast<Component&>(comp).saveCheckpoint();
            });

            // --------------------------------------------------------------------
            // PHASE 1: TRAPEZOIDAL COMPANION SUB-STEP (t -> t + gamma * dt)
            // --------------------------------------------------------------------
            ctx.time = t + gamma * dt;
            ctx.dt = dt;
            ctx.orderStep = 1;
            
            // Re-bind references to defend against pointer alterations
            ctx.currentIterV = &netlist.getStateV();
            ctx.historicStepV = &netlist.getHistoricV();
            ctx.midpointStepV = &netlist.getMidpointV();

            if (!netlist.runNewtonRaphson(ctx, 50, options.relTol, options.absTol)) {
                dt = std::max(dt * 0.5, 1e-15); 
                netlist.getStateV() = backupState; 
                netlist.inspectComponents([](const Component& comp, int extraIdx) {
                    (void)extraIdx;
                    const_cast<Component&>(comp).restoreCheckpoint();
                });
                continue;
            }
            
            // Capture midpoint state solutions
            netlist.getMidpointV() = netlist.getStateV();
            
            // Populate dynamic parameters for TR step (stores q_hist_nc / i_hist_nc)
            netlist.inspectComponents([&ctx](const Component& comp, int extraIdx) {
                (void)extraIdx;
                const_cast<Component&>(comp).updateHistory(ctx);
            });

            // --------------------------------------------------------------------
            // PHASE 2: BDF2 COMPANION BOUNDARY STEP (t + gamma * dt -> t + dt)
            // --------------------------------------------------------------------
            ctx.time = t + dt;
            ctx.orderStep = 2;
            
            // Re-bind operational vector grid positions
            ctx.currentIterV = &netlist.getStateV();
            ctx.historicStepV = &netlist.getHistoricV();
            ctx.midpointStepV = &netlist.getMidpointV();

            if (!netlist.runNewtonRaphson(ctx, 50, options.relTol, options.absTol)) {
                dt = std::max(dt * 0.5, 1e-15);
                netlist.getStateV() = backupState;
                netlist.inspectComponents([](const Component& comp, int extraIdx) {
                    (void)extraIdx;
                    const_cast<Component&>(comp).restoreCheckpoint();
                });
                continue;
            }

            // Local Truncation Error Validation Check
            double error = netlist.estimateLocalTruncationError(dt, gamma,options.relTol, options.absTol);
            if (error > options.relTol && dt > 1e-15) {
                dt = std::max(dt * 0.5, 1e-15);
                netlist.getStateV() = backupState;
                netlist.inspectComponents([](const Component& comp, int extraIdx) {
                    (void)extraIdx;
                    const_cast<Component&>(comp).restoreCheckpoint();
                });
                continue;
            }

            // Step accepted: Commit full point history profiles permanently
            netlist.inspectComponents([&ctx](const Component& comp, int extraIdx) {
                (void)extraIdx;
                const_cast<Component&>(comp).updateHistory(ctx);
            });

            netlist.getHistoricV() = netlist.getStateV();
            t += dt;

            // Stream step outputs progressively
            const auto& state = netlist.getCurrentState();
            vIn = (inputNode > 0 && static_cast<size_t>(inputNode - 1) < state.size()) ? state[inputNode - 1] : 0.0;
            vOut = (outputNode > 0 && static_cast<size_t>(outputNode - 1) < state.size()) ? state[outputNode - 1] : 0.0;

            std::cout << std::left << std::setw(12) << t 
                      << std::setw(15) << vIn 
                      << std::setw(15) << vOut << "\n";

            // Progressive adaptive timestep width adjustment bounds
            if (error < options.relTol * 0.5) {
                dt = std::min(dt * 1.5, (stopT - t));
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Transient Waveform Generation Aborted: " << e.what() << std::endl;
    }
    std::cout << "====================================================================================\n" << std::endl;
}

void ACCircuit::runFrequencySweepAnalysis(double startFreq, double stopFreq, int pointsCount) {
    std::cout << "\n==============================================================================================\n"
              << "                                  AC SIMULATION OUTPUT TABLE                                  \n"
              << "==============================================================================================\n"
              << std::left << std::setw(12) << "Freq(Hz)" 
              << std::setw(16) << "Node_2 Mag(V)" << std::setw(16) << "Node_2 Phase(°)" 
              << std::setw(16) << "Node_3 Mag(V)" << std::setw(16) << "Node_3 Phase(°)" << "\n"
              << "----------------------------------------------------------------------------------------------\n";

    std::ofstream bodeCsv("bode_plot_data.csv");
    if (!bodeCsv.is_open()) {
        throw std::runtime_error("Fatal Error: Unable to instantiate Bode data file handle.");
    }
    bodeCsv << "Frequency(Hz),V2_Magnitude(dB),V2_Phase(deg),V3_Magnitude(dB),V3_Phase(deg)\n";

    int systemSize = netlist.getTotalSystemSize();
    if (systemSize == 0) {
        netlist.buildSparsityPattern();
        systemSize = netlist.getTotalSystemSize();
    }
    // Ambil kondisi DC Operating Point sebagai fondasi linearisasi komponen non-linear (seperti Kapasitansi Dioda/MOSFET)
    const std::vector<double>& dcState = netlist.getCurrentState();

    // Initialize complex sparse matrices dynamically from the netlist structural pattern
    CsrMatrixComplex A_complex, L_complex, U_complex;
    std::vector<int> P_complex;
    A_complex.initializeFromPattern(netlist.getGMatrix());

    std::vector<std::complex<double>> B_complex(systemSize, std::complex<double>(0.0, 0.0));

    double logStart = std::log10(std::max(startFreq, 1e-3));
    double logStop = std::log10(std::max(stopFreq, 1e-3));
    double logStep = (pointsCount > 1) ? (logStop - logStart) / (pointsCount - 1) : 0.0;

    // Automated Multi-Frequency Sweep Integration Loop
    for (int i = 0; i < pointsCount; ++i) {
        double freq = std::pow(10.0, logStart + i * logStep);
        double omega = 2.0 * M_PI * freq;

        A_complex.zero();
        std::fill(B_complex.begin(), B_complex.end(), std::complex<double>(0.0, 0.0));

        // Fully Automated AC Phasor Stamping Pass across all netlist elements
        netlist.inspectComponents([&A_complex, &B_complex, omega, &dcState](const Component& comp, int extraIdx) {
            (void)extraIdx; // Bypass legacy unallocated matrix counters cleanly
            const_cast<Component&>(comp).stampAC(A_complex, B_complex, omega, dcState);
        });

        try {
            // Re-initialize factorization layers dynamically on each loop to prevent column dimension drift
            L_complex.initializeFromPattern(netlist.getGMatrix());
            U_complex.initializeFromPattern(netlist.getGMatrix());

            if (!CircuitMath::luFactorizeComplex(A_complex, L_complex, U_complex, P_complex)) {
                throw std::runtime_error("Complex LU Factorization failed at frequency: " + std::to_string(freq));
            }

            std::vector<std::complex<double>> X_complex = CircuitMath::luSolveComplex(L_complex, U_complex, P_complex, B_complex);

            // Extract target node scalar phasors cleanly (Node 2 -> Index 1, Node 3 -> Index 2)
            std::complex<double> v2_phasor = (systemSize >= 2) ? X_complex[1] : std::complex<double>(0.0, 0.0);
            std::complex<double> v3_phasor = (systemSize >= 3) ? X_complex[2] : std::complex<double>(0.0, 0.0);

            double v2_mag = std::abs(v2_phasor);
            double v2_deg = std::arg(v2_phasor) * 180.0 / M_PI;
            double v2_db  = 20.0 * std::log10(std::max(v2_mag, 1e-12));

            double v3_mag = std::abs(v3_phasor);
            double v3_deg = std::arg(v3_phasor) * 180.0 / M_PI;
            double v3_db  = 20.0 * std::log10(std::max(v3_mag, 1e-12));

            std::cout << std::left << std::setw(12) << freq
                      << std::setw(16) << v2_mag
                      << std::setw(16) << v2_deg
                      << std::setw(16) << v3_mag
                      << std::setw(16) << v3_deg << "\n";

            bodeCsv << freq << ","
                    << v2_db << "," << v2_deg << ","
                    << v3_db << "," << v3_deg << "\n";

        } catch (const std::exception& e) {
            std::cerr << "Warning: Frequency step failure at " << freq 
                      << " Hz due to solver abort: " << e.what() << std::endl;
        }
    }

    std::cout << "----------------------------------------------------------------------------------------------\n";
    std::cout << "Automated AC Sweep analysis finished. Bode datasets exported to bode_plot_data.csv" << std::endl;
}

