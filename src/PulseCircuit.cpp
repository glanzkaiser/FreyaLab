#include "freyalab/PulseCircuit.h"
#include "freyalab/CircuitComponents.h"
#include "freyalab/math.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

PulseCircuit::PulseCircuit(Netlist& netlistInstance) : netlist(netlistInstance) {}

bool PulseCircuit::isPulseCircuit() {
    if (netlist.getTotalSystemSize() == 0) {
        netlist.buildSparsityPattern();
    }

    bool hasPulseFeatures = false;
    netlist.inspectComponents([&](const Component& comp, MnaIndex extraIdx) {
        (void)extraIdx;
        ComponentType cType = comp.getType();
        if (cType == ComponentType::DIGITAL_IC ||
            cType == ComponentType::MOSFET ||
            cType == ComponentType::FINFET ||
            cType == ComponentType::RIBBONFET ||
            cType == ComponentType::GAA_FET ||
            cType == ComponentType::BJT) {
            hasPulseFeatures = true;
        }
    });

    return (netlist.getTotalSystemSize() > 0) && hasPulseFeatures;
}

void PulseCircuit::detectAndPrintCauseEffectNodes(std::vector<NetNode>& causevector,
                                                  std::vector<NetNode>& effectvector) {
    std::set<NetNode> causeNodes;
    std::set<NetNode> effectNodes;

    // 1. Inspect the topology and populate unique Cause and Effect NetNode sets
    netlist.inspectComponents([&](const Component& comp, MnaIndex extraIdx) {
        (void)extraIdx;
        ComponentType type = comp.getType();
        const auto& nodes = comp.getNodes();

        if (nodes.empty()) return;

        // Cause node discovery rule (active stimulus drivers & gate outputs)
        if (type == ComponentType::VOLTAGE_SOURCE || type == ComponentType::CURRENT_SOURCE) {
            if (!nodes[0].isGround()) causeNodes.insert(nodes[0]);
        }
        else if (type == ComponentType::DIGITAL_IC ||
                 (!comp.getName().empty() && (comp.getName()[0] == 'U' || comp.getName()[0] == 'X'))) {
            if (nodes.size() >= 2 && !nodes.back().isGround()) {
                causeNodes.insert(nodes.back());
            }
        }

        // Effect node discovery rule (reactive energy tanks + active transistor output drains)
        if (type == ComponentType::CAPACITOR || type == ComponentType::INDUCTOR) {
            for (NetNode node : nodes) {
                if (!node.isGround()) effectNodes.insert(node);
            }
        }
        else if (type == ComponentType::MOSFET ||
                 type == ComponentType::FINFET ||
                 type == ComponentType::RIBBONFET ||
                 type == ComponentType::GAA_FET ||
                 type == ComponentType::BJT) {
            // Terminal 0 is Drain (FET) or Collector (BJT)
            if (!nodes[0].isGround()) effectNodes.insert(nodes[0]);
        }
    });

    // 2. Clean up conflicts: a cause node overrides an effect node
    for (NetNode cause : causeNodes) {
        effectNodes.erase(cause);
    }

    // 3. Transfer discovered sets into output vectors
    causevector.clear();
    effectvector.clear();

    for (NetNode node : causeNodes) causevector.push_back(node);
    for (NetNode node : effectNodes) effectvector.push_back(node);

    // 4. Diagnostic display output
    std::cout << "\n========================================================================================\n";
    std::cout << "                 FREYALAB TOPOLOGICAL CAUSE & EFFECT CACHE DETECTOR                     \n";
    std::cout << "========================================================================================\n";

    std::cout << " [CAUSE NODES] (Active Forcing Elements / Instant Switch Drivers):\n";
    if (causevector.empty()) {
        std::cout << "   -> None Detected\n";
    } else {
        for (NetNode node : causevector) {
            std::cout << "   -> Node " << std::left << std::setw(3) << node.id
                      << " (" << netlist.getNodeName(node) << ") [Driven by active command stimulus]\n";
        }
    }

    std::cout << "\n [EFFECT NODES] (Dynamic Storage Networks / Active Switching Outputs):\n";
    if (effectvector.empty()) {
        std::cout << "   -> None Detected\n";
    } else {
        for (NetNode node : effectvector) {
            std::cout << "   -> Node " << std::left << std::setw(3) << node.id
                      << " (" << netlist.getNodeName(node) << ") [Dynamic energy tank / switching output]\n";
        }
    }
    std::cout << "========================================================================================\n\n";
}

void PulseCircuit::autoGenerateMonitorPairs(std::vector<std::pair<NetNode, NetNode>>& monitorPairs) {
    std::vector<NetNode> causes;
    std::vector<NetNode> effects;
    detectAndPrintCauseEffectNodes(causes, effects);

    monitorPairs.clear();

    std::cout << "========================================================================================\n"
              << "                 FREYALAB AUTOMATED TOPOLOGICAL MONITOR PAIR COMPILER                    \n"
              << "========================================================================================\n";

    std::cout << " [GENERATED CAUSE PAIRS]:\n";
    for (NetNode node : causes) {
        monitorPairs.push_back({node, NetNode{0}});
        std::cout << "   -> Paired: V(Node_" << node.id << ", GND)\n";
    }

    std::cout << "\n [GENERATED EFFECT PAIRS]:\n";
    for (NetNode node : effects) {
        monitorPairs.push_back({node, NetNode{0}});
        std::cout << "   -> Paired: V(Node_" << node.id << ", GND)\n";
    }
    std::cout << "========================================================================================\n\n";
}

bool PulseCircuit::solveTransientStep(MnaMatrixContext& ctx, double targetTime, double currentDt, int stepPhase) {
    ctx.time = targetTime;
    ctx.dt = currentDt;
    ctx.orderStep = stepPhase; // 1 = Trapezoidal Midpoint Substep, 2 = BDF2 Final Step Boundary
    constexpr int maxIterations = 50;
    constexpr double relTol = 1e-3;
    constexpr double absTol = 1e-6;
    return netlist.runNewtonRaphson(ctx, maxIterations, relTol, absTol);
}

void PulseCircuit::runPulseTransientAnalysis(const PulseSimulationOptions& options,
                                             const std::vector<std::pair<NetNode, NetNode>>& monitorPairs) {
    if (netlist.getTotalSystemSize() == 0) {
        netlist.buildSparsityPattern();
    }

    int columnsWidth = static_cast<int>(monitorPairs.size()) * 28 + 80;

    std::cout << "\n=" << std::string(columnsWidth, '=') << "\n"
              << "                        PULSE TRANSIENT SIMULATION METRIC TABLE (TR-BDF2)                        \n"
              << "=" << std::string(columnsWidth, '=') << "\n";

    std::cout << std::left << std::setw(12) << "Time(us)";
    for (const auto& pair : monitorPairs) {
        std::string label = "V(" + std::to_string(pair.first.id) + ")[V]";
        std::cout << std::setw(28) << label;
    }
    std::cout << std::setw(16) << "I_Base(A)"
              << std::setw(16) << "I_Collector(A)"
              << std::setw(16) << "I_Emitter(A)"
              << std::setw(20) << "Fill Vel.(V/us)" << "\n"
              << std::string(columnsWidth + 1, '-') << "\n";

    std::ofstream pulseCsv("pulse_transient_comprehensive_data.csv");
    if (!pulseCsv.is_open()) {
        throw std::runtime_error("Fatal Error: Unable to instantiate comprehensive pulse data output file handle.");
    }

    pulseCsv << "Time(s)";
    for (const auto& pair : monitorPairs) {
        pulseCsv << ",V_" << std::to_string(pair.first.id) << "(V)";
    }
    pulseCsv << ",I_Base(A),I_Collector(A),I_Emitter(A),Fill_Velocity(V/us)\n";

    double t = options.startTime;
    double dt = options.initialDt;
    double stopT = options.stopTime;
    const double gamma = 2.0 - std::sqrt(2.0);

    auto& G_mat = netlist.getGMatrix();
    auto& B_vec = netlist.getBVector();
    auto& state_V = netlist.getStateV();
    auto& hist_V = netlist.getHistoricV();
    auto& mid_V = netlist.getMidpointV();

    MnaMatrixContext ctx{ G_mat, B_vec, &state_V, &hist_V, &mid_V, t, dt, 0, MnaIndex{} };

    auto calculateVbicCurrents = [&](const std::vector<double>& targetedState, double& ib, double& ic, double& ie) {
        ib = 0.0;
        ic = 0.0;
        ie = 0.0;
        netlist.inspectComponents([&](const Component& comp, MnaIndex extraIdx) {
            (void)extraIdx;
            if (comp.getType() == ComponentType::BJT) {
                const auto& bjtNodes = comp.getNodes();
                if (bjtNodes.size() >= 3) {
                    NetNode nC = bjtNodes[0];
                    NetNode nB = bjtNodes[1];
                    NetNode nE = bjtNodes[2];
                    double vB = getSafeVoltage(targetedState, nB);
                    double vC = getSafeVoltage(targetedState, nC);
                    double vE = getSafeVoltage(targetedState, nE);
                    double ibe = 1e-16 * (std::exp(std::min((vB - vE) / 0.02585, 40.0)) - 1.0);
                    double ibbc = 1e-16 * (std::exp(std::min((vB - vC) / 0.02585, 40.0)) - 1.0);
                    ib = ibe + ibbc;
                    ic = ibe - ibbc;
                    ie = -(ib + ic);
                }
            }
        });
    };

    try {
        // --- PHASE 0: INITIAL OPERATING POINT STATE ---
        ctx.orderStep = 0;
        ctx.time = t;
        ctx.dt = 0.0;
        if (!netlist.runNewtonRaphson(ctx, 100, options.relTol, options.absTol)) {
            throw std::runtime_error("DC operating point convergence failed.");
        }
        hist_V = state_V;
        mid_V = state_V;
        netlist.inspectComponents([&](Component& comp, MnaIndex extraIdx) {
            (void)extraIdx;
            comp.updateHistory(ctx);
        });

        double i_b = 0.0, i_c = 0.0, i_e = 0.0;
        calculateVbicCurrents(state_V, i_b, i_c, i_e);
        double initialFillVelocity = 0.0;

        std::cout << std::left << std::setw(12) << (t * 1e6);
        pulseCsv << t;
        for (const auto& pair : monitorPairs) {
            double vDiff = getSafeVoltage(state_V, pair.first) - getSafeVoltage(state_V, pair.second);
            std::cout << std::setw(28) << vDiff;
            pulseCsv << "," << vDiff;
        }
        std::cout << std::setw(16) << i_b
                  << std::setw(16) << i_c
                  << std::setw(16) << i_e
                  << std::setw(20) << initialFillVelocity << "\n";
        pulseCsv << "," << i_b << "," << i_c << "," << i_e << "," << initialFillVelocity << "\n";

        // --- TIME INTEGRATION TRANSITION TIMELINE LOOP ---
        constexpr double TIME_SLACK = 1e-15;
        while (t < (stopT - TIME_SLACK)) {
            if (t + dt > stopT) {
                dt = stopT - t;
            }

            std::vector<double> backupState = state_V;
            std::vector<double> backupMid   = mid_V;

            double vOutOld = 0.0;
            if (!monitorPairs.empty()) {
                vOutOld = getSafeVoltage(state_V, monitorPairs.back().first) -
                          getSafeVoltage(state_V, monitorPairs.back().second);
            }

            netlist.inspectComponents([](Component& comp, MnaIndex extraIdx) {
                (void)extraIdx;
                comp.saveCheckpoint();
            });

            // PHASE 1: TRAPEZOIDAL SUB-STEP (t -> t + gamma * dt)
            ctx.time = t + gamma * dt;
            ctx.dt = dt;
            ctx.orderStep = 1;
            if (!netlist.runNewtonRaphson(ctx, options.maxIterations, options.relTol, options.absTol)) {
                dt = std::max(dt * 0.5, options.minDt);
                state_V = backupState;
                mid_V = backupMid;
                netlist.inspectComponents([](Component& comp, MnaIndex extraIdx) {
                    (void)extraIdx;
                    comp.restoreCheckpoint();
                });
                continue;
            }
            mid_V = state_V;
            netlist.inspectComponents([&ctx](Component& comp, MnaIndex extraIdx) {
                (void)extraIdx;
                comp.updateHistory(ctx);
            });

            // PHASE 2: BDF2 STEP BOUNDARY (t + gamma * dt -> t + dt)
            ctx.time = t + dt;
            ctx.dt = dt;
            ctx.orderStep = 2;
            if (!netlist.runNewtonRaphson(ctx, options.maxIterations, options.relTol, options.absTol)) {
                dt = std::max(dt * 0.5, options.minDt);
                state_V = backupState;
                mid_V = backupMid;
                netlist.inspectComponents([](Component& comp, MnaIndex extraIdx) {
                    (void)extraIdx;
                    comp.restoreCheckpoint();
                });
                continue;
            }

            // LTE Validation Pass
            double error = netlist.estimateLocalTruncationError(dt, gamma, options.relTol, options.absTol);
            if (error > 1.0) {
                dt = std::max(dt * 0.5, options.minDt);
                state_V = backupState;
                mid_V = backupMid;
                netlist.inspectComponents([](Component& comp, MnaIndex extraIdx) {
                    (void)extraIdx;
                    comp.restoreCheckpoint();
                });
                continue;
            }

            // Step Accepted: Commit history updates
            netlist.inspectComponents([&ctx](Component& comp, MnaIndex extraIdx) {
                (void)extraIdx;
                comp.updateHistory(ctx);
            });
            hist_V = state_V;

            double vOutNew = 0.0;
            if (!monitorPairs.empty()) {
                vOutNew = getSafeVoltage(state_V, monitorPairs.back().first) -
                          getSafeVoltage(state_V, monitorPairs.back().second);
            }
            double fillVelocity = (dt > 0.0) ? (((vOutNew - vOutOld) / dt) * 1e-6) : 0.0;

            t += dt;
            ctx.time = t;
            calculateVbicCurrents(state_V, i_b, i_c, i_e);

            std::cout << std::left << std::setw(12) << (t * 1e6);
            pulseCsv << t;
            for (const auto& pair : monitorPairs) {
                double vDiff = getSafeVoltage(state_V, pair.first) - getSafeVoltage(state_V, pair.second);
                std::cout << std::setw(28) << vDiff;
                pulseCsv << "," << vDiff;
            }
            std::cout << std::setw(16) << i_b
                      << std::setw(16) << i_c
                      << std::setw(16) << i_e
                      << std::setw(20) << fillVelocity << "\n";
            pulseCsv << "," << i_b << "," << i_c << "," << i_e << "," << fillVelocity << "\n";

            if (error > 1e-12) {
                double hFactor = 0.9 * std::pow(1.0 / error, 1.0 / 3.0);
                dt = std::min(dt * std::max(0.5, std::min(1.5, hFactor)), options.maxDt);
            } else {
                dt = std::min(dt * 1.5, options.maxDt);
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "\n[FATAL] Pulse Transient Matrix Sweep Loop Aborted: " << e.what() << std::endl;
    }
    std::cout << "=" << std::string(columnsWidth, '=') << "\n" << std::endl;
}

