#include "freyalab/ACCircuit.h"
#include "freyalab/math.h"
#include "freyalab/KLU.h"
#include "freyalab/CircuitComponents.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

ACCircuit::ACCircuit(Netlist& netlistInstance) : netlist(netlistInstance) {}

bool ACCircuit::isAcCircuit() {
    if (netlist.getTotalSystemSize() == 0) {
        netlist.buildSparsityPattern();
    }
    return netlist.getTotalSystemSize() > 0;
}

void ACCircuit::autoGenerateMonitorPairs(std::vector<std::pair<NetNode, NetNode>>& monitorPairs) {
    if (netlist.getTotalSystemSize() == 0) {
        netlist.buildSparsityPattern();
    }

    std::set<NetNode> causeNodes;
    std::set<NetNode> effectNodes;

    // 1. Pindai topologi untuk mendeteksi node Cause (sumber stimulus) dan Effect (beban/output)
    netlist.inspectComponents([&](const Component& comp, MnaIndex extraIdx) {
        (void)extraIdx;
        ComponentType type = comp.getType();
        const auto& nodes = comp.getNodes();
        if (nodes.empty()) return;

        // Cause nodes: sumber AC/DC tegangan & arus atau output gerbang logika
        if (type == ComponentType::VOLTAGE_SOURCE ||
            type == ComponentType::AC_SOURCE_DEVICE ||
            type == ComponentType::CURRENT_SOURCE) {
            if (!nodes[0].isGround()) causeNodes.insert(nodes[0]);
            if (nodes.size() >= 2 && !nodes[1].isGround()) causeNodes.insert(nodes[1]);
        } else if (type == ComponentType::DIGITAL_IC) {
            if (nodes.size() >= 2 && !nodes.back().isGround()) {
                causeNodes.insert(nodes.back());
            }
        }

        // Effect nodes: elemen pasif/reaktif (C, L, R) atau terminal output transistor (Drain/Collector)
        if (type == ComponentType::CAPACITOR ||
            type == ComponentType::INDUCTOR ||
            type == ComponentType::RESISTOR) {
            for (NetNode n : nodes) {
                if (!n.isGround()) effectNodes.insert(n);
            }
        } else if (type == ComponentType::MOSFET ||
                   type == ComponentType::FINFET ||
                   type == ComponentType::RIBBONFET ||
                   type == ComponentType::GAA_FET ||
                   type == ComponentType::BJT) {
            if (!nodes[0].isGround()) effectNodes.insert(nodes[0]);
        }
    });

    // 2. Prioritaskan Cause nodes agar tidak duplikat di daftar Effect
    for (NetNode c : causeNodes) {
        effectNodes.erase(c);
    }

    // 3. Susun pasangan pengukuran diferensial terhadap Ground (NetNode{0})
    monitorPairs.clear();
    for (NetNode c : causeNodes) {
        monitorPairs.push_back({c, NetNode{0}});
    }
    for (NetNode e : effectNodes) {
        monitorPairs.push_back({e, NetNode{0}});
    }

    // Fallback: jika tidak ada node spesifik terdeteksi, ukur seluruh node eksternal 1..numNodes
    if (monitorPairs.empty()) {
        int numNodes = netlist.getNumNodes();
        for (int i = 1; i <= numNodes; ++i) {
            monitorPairs.push_back({NetNode{i}, NetNode{0}});
        }
    }
}

bool ACCircuit::solveTransientStep(MnaMatrixContext& ctx, double targetTime, double currentDt, int stepPhase) {
    ctx.time = targetTime;
    ctx.dt = currentDt;
    ctx.orderStep = stepPhase;
    return netlist.runNewtonRaphson(ctx, 50, 1e-3, 1e-6);
}

void ACCircuit::runTransientWaveformAnalysis(const AcSimulationOptions& options) {
    std::vector<std::pair<NetNode, NetNode>> autoPairs;
    autoGenerateMonitorPairs(autoPairs);
    runTransientWaveformAnalysis(options, autoPairs);
}

void ACCircuit::runTransientWaveformAnalysis(const AcSimulationOptions& options,
                                             NetNode inputNode,
                                             NetNode outputNode) {
    std::vector<std::pair<NetNode, NetNode>> pairs = {
        {inputNode, NetNode{0}},
        {outputNode, NetNode{0}}
    };
    runTransientWaveformAnalysis(options, pairs);
}

void ACCircuit::runTransientWaveformAnalysis(const AcSimulationOptions& options,
                                             const std::vector<std::pair<NetNode, NetNode>>& monitorPairs) {
    if (netlist.getTotalSystemSize() == 0) {
        netlist.buildSparsityPattern();
    }

    auto formatPairLabel = [this](const std::pair<NetNode, NetNode>& pair) -> std::string {
        std::string posName = pair.first.isGround() ? "GND" : netlist.getNodeName(pair.first);
        if (pair.second.isGround()) {
            return "V(" + posName + ")";
        }
        return "V(" + posName + "," + netlist.getNodeName(pair.second) + ")";
    };

    int tableWidth = 14 + static_cast<int>(monitorPairs.size()) * 18;
    std::string borderLine(std::max(tableWidth, 60), '=');
    std::string dashLine(std::max(tableWidth, 60), '-');

    std::cout << "\n" << borderLine << "\n"
              << "                        AC TRANSIENT WAVEFORM TABLE (2ND-ORDER TR-BDF2)            \n"
              << borderLine << "\n"
              << std::left << std::setw(14) << "Time(s)";

    for (const auto& pair : monitorPairs) {
        std::cout << std::setw(18) << formatPairLabel(pair);
    }
    std::cout << "\n" << dashLine << "\n";

    std::ofstream acTransCsv("ac_transient_waveform_data.csv");
    if (acTransCsv.is_open()) {
        acTransCsv << "Time(s)";
        for (const auto& pair : monitorPairs) {
            acTransCsv << "," << formatPairLabel(pair);
        }
        acTransCsv << "\n";
    }

    double t = options.startTime;
    double dt = options.initialDt;
    double stopT = options.stopTime;
    const double gamma = 2.0 - std::sqrt(2.0);

    MnaMatrixContext ctx{
        netlist.getGMatrix(),
        netlist.getBVector(),
        &netlist.getStateV(),
        &netlist.getHistoricV(),
        &netlist.getMidpointV(),
        t, dt, 0, MnaIndex{},
	netlist.getTemperatureCelsius(),
        netlist.getThermalVoltage()
    };

    try {
        // STEP 1: INITIAL DC OPERATING POINT
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
            ctx.assignedExtraIdx = extraIdx;
            comp.updateHistory(ctx);
        });

        const auto& initialState = netlist.getCurrentState();
        std::cout << std::left << std::setw(14) << t;
        if (acTransCsv.is_open()) acTransCsv << t;
        for (const auto& pair : monitorPairs) {
            double vDiff = getSafeVoltage(initialState, pair.first) - getSafeVoltage(initialState, pair.second);
            std::cout << std::setw(18) << vDiff;
            if (acTransCsv.is_open()) acTransCsv << "," << vDiff;
        }
        std::cout << "\n";
        if (acTransCsv.is_open()) acTransCsv << "\n";

        // STEP 2: TR-BDF2 TIME STEPPING LOOP
        while (t + options.minDt * 0.5 < stopT) {
            if (t + dt > stopT) {
                dt = stopT - t;
            }
            if (dt < options.minDt) {
                throw std::runtime_error("Timestep crashed below minimum safety limit.");
            }

            std::vector<double> backupState = netlist.getStateV();
            netlist.inspectComponents([](Component& comp, MnaIndex extraIdx) {
                (void)extraIdx;
                comp.saveCheckpoint();
            });

            // PHASE 1: TRAPEZOIDAL MIDPOINT (t -> t + gamma * dt)
            ctx.time = t + gamma * dt;
            ctx.dt = dt;
            ctx.orderStep = 1;
            ctx.currentIterV = &netlist.getStateV();
            ctx.historicStepV = &netlist.getHistoricV();
            ctx.midpointStepV = nullptr;

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
                ctx.assignedExtraIdx = extraIdx;
                comp.updateHistory(ctx);
            });

            // PHASE 2: BDF2 COMPOSITE STEP (t + gamma * dt -> t + dt)
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

            // Normalized LTE check (<= 1.0 accepted)
            double error = netlist.estimateLocalTruncationError(dt, gamma, options.relTol, options.absTol);
            if (error > 1.0) {
                dt *= 0.5;
                netlist.getStateV() = backupState;
                netlist.inspectComponents([](Component& comp, MnaIndex extraIdx) {
                    (void)extraIdx;
                    comp.restoreCheckpoint();
                });
                continue;
            }

            netlist.inspectComponents([&ctx](Component& comp, MnaIndex extraIdx) {
                ctx.assignedExtraIdx = extraIdx;
                comp.updateHistory(ctx);
            });

            netlist.getHistoricV() = netlist.getStateV();
            t += dt;

            const auto& state = netlist.getCurrentState();
            std::cout << std::left << std::setw(14) << t;
            if (acTransCsv.is_open()) acTransCsv << t;
            for (const auto& pair : monitorPairs) {
                double vDiff = getSafeVoltage(state, pair.first) - getSafeVoltage(state, pair.second);
                std::cout << std::setw(18) << vDiff;
                if (acTransCsv.is_open()) acTransCsv << "," << vDiff;
            }
            std::cout << "\n";
            if (acTransCsv.is_open()) acTransCsv << "\n";

            if (error < 0.5) {
                dt = std::min(dt * 1.5, options.maxDt);
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Transient Waveform Generation Aborted: " << e.what() << std::endl;
    }
    std::cout << dashLine << "\n" << std::endl;
}

void ACCircuit::runFrequencySweepAnalysis(double startFreq, double stopFreq, int pointsCount) {
    std::vector<std::pair<NetNode, NetNode>> autoPairs;
    autoGenerateMonitorPairs(autoPairs);
    runFrequencySweepAnalysis(startFreq, stopFreq, pointsCount, autoPairs);
}

void ACCircuit::runFrequencySweepAnalysis(double startFreq,
                                          double stopFreq,
                                          int pointsCount,
                                          const std::vector<std::pair<NetNode, NetNode>>& monitorPairs) {
    int systemSize = netlist.getTotalSystemSize();
    if (systemSize == 0) {
        netlist.buildSparsityPattern();
        systemSize = netlist.getTotalSystemSize();
    }

    // Helper untuk membentuk label kolom: "Node_2" (single-ended) atau "V(2,3)" (diferensial)
    auto formatPairLabel = [this](const std::pair<NetNode, NetNode>& pair) -> std::string {
        std::string posName = pair.first.isGround() ? "GND" : netlist.getNodeName(pair.first);
        if (pair.second.isGround()) {
            return posName;
        }
        std::string negName = netlist.getNodeName(pair.second);
        return "V(" + posName + "," + negName + ")";
    };

    int tableWidth = 14 + static_cast<int>(monitorPairs.size()) * 32;
    std::string borderLine(std::max(tableWidth, 60), '=');
    std::string dashLine(std::max(tableWidth, 60), '-');

    std::cout << "\n" << borderLine << "\n"
              << "                                  AC SIMULATION OUTPUT TABLE                                  \n"
              << borderLine << "\n"
              << std::left << std::setw(14) << "Freq(Hz)";

    for (const auto& pair : monitorPairs) {
        std::string lbl = formatPairLabel(pair);
        std::cout << std::setw(16) << (lbl + " Mag(V)")
                  << std::setw(16) << (lbl + " Phase(deg)");
    }
    std::cout << "\n" << dashLine << "\n";

    std::ofstream bodeCsv("bode_plot_data.csv");
    if (!bodeCsv.is_open()) {
        throw std::runtime_error("Fatal Error: Unable to instantiate Bode data file handle.");
    }

    bodeCsv << "Frequency(Hz)";
    for (const auto& pair : monitorPairs) {
        std::string lbl = formatPairLabel(pair);
        bodeCsv << "," << lbl << "_Magnitude(dB)," << lbl << "_Phase(deg)";
    }
    bodeCsv << "\n";

    // 1. Selesaikan DC Operating Point terlebih dahulu untuk linearisasi komponen nonlinier
    MnaMatrixContext dcCtx{
        netlist.getGMatrix(),
        netlist.getBVector(),
        &netlist.getStateV(),
        &netlist.getHistoricV(),
        &netlist.getMidpointV(),
        0.0,   // time
        0.0,   // dt (0.0 untuk DC statis)
        0,     // orderStep = 0 (DC Operating Point)
        MnaIndex{},
        netlist.getTemperatureCelsius(),
        netlist.getThermalVoltage()
    };

    if (!netlist.runNewtonRaphson(dcCtx, 150, 1e-3, 1e-6)) {
        std::cerr << "Warning: DC Operating Point did not fully converge prior to AC Sweep; "
                  << "linearizing around best available state.\n";
    }

    const std::vector<double> dcState = netlist.getCurrentState();

    // 2. Inisialisasi matriks kompleks dari pola spasial G
    CsrMatrixComplex A_complex, L_complex, U_complex;
    std::vector<int> P_complex;
    A_complex.initializeFromPattern(netlist.getGMatrix());

    bool useKlu = netlist.isUsingKluBackend();
    if (!useKlu) {
        L_complex.initializeFromPattern(netlist.getGMatrix());
        U_complex.initializeFromPattern(netlist.getGMatrix());
        P_complex.resize(systemSize, 0);
    }

    std::vector<std::complex<double>> B_complex(systemSize, std::complex<double>(0.0, 0.0));

    double logStart = std::log10(std::max(startFreq, 1e-3));
    double logStop  = std::log10(std::max(stopFreq, 1e-3));
    double logStep  = (pointsCount > 1) ? (logStop - logStart) / (pointsCount - 1) : 0.0;

    int numNodes = netlist.getNumNodes();
    auto getSafePhasor = [numNodes](const std::vector<std::complex<double>>& X, NetNode node) -> std::complex<double> {
        if (node.isGround()) return std::complex<double>(0.0, 0.0);
        MnaIndex idx = MnaIndex::fromNetNode(node);
        if (idx.isValid() && idx.idx < numNodes && static_cast<size_t>(idx.idx) < X.size()) {
            return X[idx.idx];
        }
        return std::complex<double>(0.0, 0.0);
    };

    // 3. Loop sapuan frekuensi logaritmik
    for (int i = 0; i < pointsCount; ++i) {
        double freq  = std::pow(10.0, logStart + i * logStep);
        double omega = 2.0 * M_PI * freq;

        A_complex.zero();
        std::fill(B_complex.begin(), B_complex.end(), std::complex<double>(0.0, 0.0));

        netlist.inspectComponents([&A_complex, &B_complex, omega, &dcState](Component& comp, MnaIndex extraIdx) {
            (void)extraIdx;
            comp.stampAC(A_complex, B_complex, omega, dcState);
        });

        try {
            std::vector<std::complex<double>> X_complex;

            if (useKlu) {
                if (!KluSolver::factorizeNumericComplex(A_complex, netlist.getKluCache())) {
                    throw std::runtime_error("Complex KLU Factorization failed at frequency: " + std::to_string(freq));
                }
                X_complex = KluSolver::solveNumericComplex(netlist.getKluCache(), B_complex);
            } else {
                P_complex.assign(systemSize, 0);
                if (!CircuitMath::luFactorizeComplex(A_complex, L_complex, U_complex, P_complex)) {
                    throw std::runtime_error("Complex LU Factorization failed at frequency: " + std::to_string(freq));
                }
                X_complex = CircuitMath::luSolveComplex(L_complex, U_complex, P_complex, B_complex);
            }

            std::cout << std::left << std::setw(14) << freq;
            bodeCsv << freq;

            for (const auto& pair : monitorPairs) {
                std::complex<double> vDiff = getSafePhasor(X_complex, pair.first) - getSafePhasor(X_complex, pair.second);
                double mag = std::abs(vDiff);
                double deg = std::arg(vDiff) * 180.0 / M_PI;
                double db  = 20.0 * std::log10(std::max(mag, 1e-12));

                std::cout << std::setw(16) << mag << std::setw(16) << deg;
                bodeCsv << "," << db << "," << deg;
            }

            std::cout << "\n";
            bodeCsv << "\n";

        } catch (const std::exception& e) {
            std::cerr << "Warning: Frequency step failure at " << freq
                      << " Hz due to solver abort: " << e.what() << std::endl;
        }
    }

    std::cout << dashLine << "\n";
    std::cout << "Automated AC Sweep analysis finished. Bode datasets exported to bode_plot_data.csv" << std::endl;
}


