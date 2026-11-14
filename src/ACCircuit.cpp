#include "freyalab/ACCircuit.h"
#include "freyalab/math.h"
#include "freyalab/CircuitComponents.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
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

void ACCircuit::runTransientWaveformAnalysis(const AcSimulationOptions& options,
                                             NetNode inputNode,
                                             NetNode outputNode) {
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

    MnaMatrixContext ctx{
        netlist.getGMatrix(),
        netlist.getBVector(),
        &netlist.getStateV(),
        &netlist.getHistoricV(),
        &netlist.getMidpointV(),
        t, dt, 0, MnaIndex{}
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

        netlist.getHistoricV() = netlist.getStateV();

        netlist.inspectComponents([&ctx](Component& comp, MnaIndex extraIdx) {
            (void)extraIdx;
            comp.updateHistory(ctx);
        });

        const auto& initialState = netlist.getCurrentState();
        double vIn  = getSafeVoltage(initialState, inputNode);
        double vOut = getSafeVoltage(initialState, outputNode);
        std::cout << std::left << std::setw(12) << t << std::setw(15) << vIn << std::setw(15) << vOut << "\n";

        // ====================================================================
        // STEP 2: TIME STEPPING SIMULATION LOOP
        // ====================================================================
        while (t + options.minDt * 0.5 < stopT) {
            if (t + dt > stopT) {
                dt = stopT - t;
            }

            if (dt < options.minDt) {
                throw std::runtime_error("Timestep crashed below minimum safety limit (stagnant matrix grid configuration).");
            }

            std::vector<double> backupState = netlist.getStateV();

            netlist.inspectComponents([](Component& comp, MnaIndex extraIdx) {
                (void)extraIdx;
                comp.saveCheckpoint();
            });

            // --------------------------------------------------------------------
            // PHASE 1: TRAPEZOIDAL COMPANION SUB-STEP (t -> t + gamma * dt)
            // --------------------------------------------------------------------
            ctx.time = t + gamma * dt;
            ctx.dt = dt;
            ctx.orderStep = 1;

            ctx.currentIterV = &netlist.getStateV();
            ctx.historicStepV = &netlist.getHistoricV();
            ctx.midpointStepV = &netlist.getMidpointV();

            if (!netlist.runNewtonRaphson(ctx, options.maxIterations, options.relTol, options.absTol)) {
                dt *= 0.5;
                netlist.getStateV() = backupState;
                netlist.inspectComponents([](Component& comp, MnaIndex extraIdx) {
                    (void)extraIdx;
                    comp.restoreCheckpoint();
                });
                continue;
            }

            netlist.getMidpointV() = netlist.getStateV();

            netlist.inspectComponents([&ctx](Component& comp, MnaIndex extraIdx) {
                (void)extraIdx;
                comp.updateHistory(ctx);
            });

            // --------------------------------------------------------------------
            // PHASE 2: BDF2 COMPANION BOUNDARY STEP (t + gamma * dt -> t + dt)
            // --------------------------------------------------------------------
            ctx.time = t + dt;
            ctx.orderStep = 2;

            ctx.currentIterV = &netlist.getStateV();
            ctx.historicStepV = &netlist.getHistoricV();
            ctx.midpointStepV = &netlist.getMidpointV();

            if (!netlist.runNewtonRaphson(ctx, options.maxIterations, options.relTol, options.absTol)) {
                dt *= 0.5;
                netlist.getStateV() = backupState;
                netlist.inspectComponents([](Component& comp, MnaIndex extraIdx) {
                    (void)extraIdx;
                    comp.restoreCheckpoint();
                });
                continue;
            }

            // Local Truncation Error Validation Check
            double error = netlist.estimateLocalTruncationError(dt, gamma, options.relTol, options.absTol);
            if (error > options.relTol) {
                dt *= 0.5;
                netlist.getStateV() = backupState;
                netlist.inspectComponents([](Component& comp, MnaIndex extraIdx) {
                    (void)extraIdx;
                    comp.restoreCheckpoint();
                });
                continue;
            }

            // Step accepted: Commit step updates permanently
            netlist.inspectComponents([&ctx](Component& comp, MnaIndex extraIdx) {
                (void)extraIdx;
                comp.updateHistory(ctx);
            });

            netlist.getHistoricV() = netlist.getStateV();
            t += dt;

            // Stream step outputs progressively via bounds-checked NetNode lookups
            const auto& state = netlist.getCurrentState();
            vIn  = getSafeVoltage(state, inputNode);
            vOut = getSafeVoltage(state, outputNode);

            std::cout << std::left << std::setw(12) << t
                      << std::setw(15) << vIn
                      << std::setw(15) << vOut << "\n";

            if (error < options.relTol * 0.5) {
                dt = std::min(dt * 1.5, options.maxDt);
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

    // Solve the DC Operating Point first so nonlinear devices (MOSFET, BJT, Diode,
    // and BSIM-CMG FinFET/RibbonFET/GAA-FET) are linearized around their true DC bias
    MnaMatrixContext dcCtx{
        netlist.getGMatrix(),
        netlist.getBVector(),
        &netlist.getStateV(),
        &netlist.getHistoricV(),
        &netlist.getMidpointV(),
        0.0,   // time
        0.0,   // dt (0.0 for static DC)
        0,     // orderStep = 0 (DC Operating Point)
        MnaIndex{}
    };

    if (!netlist.runNewtonRaphson(dcCtx, 150, 1e-3, 1e-6)) {
        std::cerr << "Warning: DC Operating Point did not fully converge prior to AC Sweep; "
                  << "linearizing around best available state.\n";
    }

    const std::vector<double> dcState = netlist.getCurrentState();

    // Initialize complex sparse matrices dynamically from the netlist structural pattern
    CsrMatrixComplex A_complex, L_complex, U_complex;
    std::vector<int> P_complex;
    A_complex.initializeFromPattern(netlist.getGMatrix());
    L_complex.initializeFromPattern(netlist.getGMatrix());
    U_complex.initializeFromPattern(netlist.getGMatrix());
    std::vector<std::complex<double>> B_complex(systemSize, std::complex<double>(0.0, 0.0));

    double logStart = std::log10(std::max(startFreq, 1e-3));
    double logStop  = std::log10(std::max(stopFreq, 1e-3));
    double logStep  = (pointsCount > 1) ? (logStop - logStart) / (pointsCount - 1) : 0.0;

    // Pre-resolve type-safe MnaIndex coordinates for Node 2 and Node 3
    const MnaIndex mNode2 = MnaIndex::fromNetNode(NetNode{2});
    const MnaIndex mNode3 = MnaIndex::fromNetNode(NetNode{3});

    // Automated Multi-Frequency Sweep Integration Loop
    for (int i = 0; i < pointsCount; ++i) {
        double freq  = std::pow(10.0, logStart + i * logStep);
        double omega = 2.0 * M_PI * freq;

        A_complex.zero();
        L_complex.zero();
        U_complex.zero();
        std::fill(B_complex.begin(), B_complex.end(), std::complex<double>(0.0, 0.0));

        // Stamp small-signal AC phasor admittances across all netlist elements
        netlist.inspectComponents([&A_complex, &B_complex, omega, &dcState](Component& comp, MnaIndex extraIdx) {
            (void)extraIdx;
            comp.stampAC(A_complex, B_complex, omega, dcState);
        });

        try {
            L_complex.initializeFromPattern(netlist.getGMatrix());
            U_complex.initializeFromPattern(netlist.getGMatrix());

            if (!CircuitMath::luFactorizeComplex(A_complex, L_complex, U_complex, P_complex)) {
                throw std::runtime_error("Complex LU Factorization failed at frequency: " + std::to_string(freq));
            }

            std::vector<std::complex<double>> X_complex = CircuitMath::luSolveComplex(L_complex, U_complex, P_complex, B_complex);

            int numNodes = netlist.getNumNodes();
            std::complex<double> v2_phasor = (mNode2.isValid() && mNode2.idx < numNodes && static_cast<size_t>(mNode2.idx) < X_complex.size())
                                           ? X_complex[mNode2.idx] : std::complex<double>(0.0, 0.0);
            std::complex<double> v3_phasor = (mNode3.isValid() && mNode3.idx < numNodes && static_cast<size_t>(mNode3.idx) < X_complex.size())
                                           ? X_complex[mNode3.idx] : std::complex<double>(0.0, 0.0);

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


