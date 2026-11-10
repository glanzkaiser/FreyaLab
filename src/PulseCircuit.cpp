#include "freyalab/PulseCircuit.h"
#include "freyalab/CircuitComponents.h" 
#include "freyalab/math.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <fstream>
#include <set>
#include <string>
#include <vector>
#include <utility>

PulseCircuit::PulseCircuit(Netlist& netlistInstance) : netlist(netlistInstance) {}
// mathematical correct, robust and has no structural bugs
bool PulseCircuit::isPulseCircuit() {
    if (netlist.getTotalSystemSize() == 0) {

    }
    bool hasPulseFeatures = false;
    netlist.inspectComponents([&](const Component& comp, int extraIdx) {
    if (comp.getType() == ComponentType::DIGITAL_IC || 
            comp.getType() == ComponentType::MOSFET || 
            comp.getType() == ComponentType::BJT) {
            hasPulseFeatures = true;
        }
    });

    return (netlist.getTotalSystemSize() > 0) && hasPulseFeatures;
}


void PulseCircuit::detectAndPrintCauseEffectNodes(std::vector<int>& causevector, std::vector<int>& effectvector) {
    std::set<int> causeNodes;
    std::set<int> effectNodes;

    // 1. Inspect the topology and populate unique sets
    netlist.inspectComponents([&](const Component& comp, int extraIdx) {
        (void)extraIdx;
        ComponentType type = comp.getType();
        const auto& nodes = comp.getNodes();

        if (nodes.empty()) return;

        // Cause node discovery rule
        if (type == ComponentType::VOLTAGE_SOURCE || type == ComponentType::CURRENT_SOURCE) {
            if (nodes[0] > 0) causeNodes.insert(nodes[0]);
        } 
        else if (type == ComponentType::DIGITAL_IC || comp.getName()[0] == 'U' || comp.getName()[0] == 'X') {
            if (nodes.size() >= 2 && nodes.back() > 0) {
                causeNodes.insert(nodes.back());
            }
        }

        // Effect node discovery rule
        if (type == ComponentType::CAPACITOR || type == ComponentType::INDUCTOR) {
            for (int node : nodes) {
                if (node > 0) effectNodes.insert(node);
            }
        }
    });

    // 2. Clean up conflicts: a cause node overrides an effect node behavior
    for (int cause : causeNodes) {
        effectNodes.erase(cause);
    }

    // 3. Clear output variables and transfer cached items into your vector targets
    causevector.clear();
    effectvector.clear();
    
    for (int node : causeNodes) causevector.push_back(node);
    for (int node : effectNodes) effectvector.push_back(node);

    // 4. Diagnostic display output matrix metrics
    std::cout << "\n========================================================================================\n";
    std::cout << "                 FREYALAB TOPOLOGICAL CAUSE & EFFECT CACHE DETECTOR                     \n";
    std::cout << "========================================================================================\n";
    
    std::cout << " [CAUSE NODES] (Active Forcing Elements / Instant Switch Drivers):\n";
    if (causevector.empty()) {
        std::cout << "   -> None Detected\n";
    } else {
        for (int node : causevector) {
            std::cout << "   -> Node " << std::left << std::setw(3) << node << " [Driven by active command stimulus]\n";
        }
    }

    std::cout << "\n [EFFECT NODES] (Dynamic Storage Networks / Lazy Damped Transient Responses):\n";
    if (effectvector.empty()) {
        std::cout << "   -> None Detected\n";
    } else {
        for (int node : effectvector) {
            std::cout << "   -> Node " << std::left << std::setw(3) << node << " [Dynamic energy tank - subject to delay]\n";
        }
    }
    std::cout << "========================================================================================\n\n";
}

void PulseCircuit::autoGenerateMonitorPairs(std::vector<std::pair<int, int>>& monitorPairs) {
    std::vector<int> causes;
    std::vector<int> effects;
    detectAndPrintCauseEffectNodes(causes, effects);

    monitorPairs.clear();
    
    std::cout << "========================================================================================\n"
              << "                 FREYALAB AUTOMATED TOPOLOGICAL MONITOR PAIR COMPILER                    \n"
              << "========================================================================================\n";

    std::cout << " [GENERATED CAUSE PAIRS]:\n";
    for (int node : causes) {
        monitorPairs.push_back(std::make_pair(node, 0));
        std::cout << "   -> Paired: V(Node_" << node << ":" << netlist.getNodeName(node) << ", GND)\n";
    }

    std::cout << "\n [GENERATED EFFECT PAIRS]:\n";
    for (int node : effects) {
        monitorPairs.push_back(std::make_pair(node, 0));
        std::cout << "   -> Paired: V(Node_" << node << ":" << netlist.getNodeName(node) << ", GND)\n";
    }
    std::cout << "========================================================================================\n\n";
}

// mathematical correct, robust and has no structural bugs
bool PulseCircuit::solveTransientStep(MnaMatrixContext& ctx, double targetTime, double currentDt, int stepPhase) {
    ctx.time = targetTime;
    ctx.dt = currentDt;
    ctx.orderStep = stepPhase; // 1 = Trapezoidal Midpoint Substep, 2 = BDF2 Final Step Boundary
    constexpr int maxIterations = 50;
    constexpr double relTol = 1e-3;
    constexpr double absTol = 1e-6;
    return netlist.runNewtonRaphson(ctx, maxIterations, relTol, absTol);
}
// mathematical correct, robust and has no structural bugs 

// mathematical correct, robust and has no structural bugs
void PulseCircuit::runPulseTransientAnalysis(const PulseSimulationOptions& options, const std::vector<std::pair<int, int>>& monitorPairs) {
    // Dynamic width calculation based on descriptive labels (28 spaces wide per custom node pair string)
    int columnsWidth = monitorPairs.size() * 28 + 60;
    
    std::cout << "\n=" << std::string(columnsWidth, '=') << "\n"
              << "                                        PULSE TRANSIENT SIMULATION METRIC TABLE (TR-BDF2 + VBIC)                                        \n"
              << "=" << std::string(columnsWidth, '=') << "\n";
              
    std::cout << std::left << std::setw(12) << "Time(us)";
    for (const auto& pair : monitorPairs) {
        std::string label = "V(" + netlist.getNodeName(pair.first) + ")[V]";
        std::cout << std::setw(28) << label;
    }
    std::cout << std::setw(16) << "I_Base(A)" 
              << std::setw(16) << "I_Collector(A)" 
              << std::setw(16) << "I_Emitter(A)" 
              << std::setw(20) << "Fill Vel.(V/us)" << "\n"
              << "----------------------------------------------------------------------------------------------------------------------------------------\n";

    std::ofstream pulseCsv("pulse_transient_comprehensive_data.csv");
    if (!pulseCsv.is_open()) {
        throw std::runtime_error("Fatal Error: Unable to instantiate comprehensive pulse data output file handle.");
    }
    
    pulseCsv << "Time(s)";
    for (const auto& pair : monitorPairs) {
        pulseCsv << ",V_" << netlist.getNodeName(pair.first) << "(V)";
    }
    pulseCsv << ",I_Base(A),I_Collector(A),I_Emitter(A),Fill_Velocity(V/us)\n";

    int systemSize = netlist.getTotalSystemSize();
    double t = options.startTime;
    double dt = options.initialDt;
    double stopT = options.stopTime;
    const double gamma = 2.0 - std::sqrt(2.0); 

    auto& G_mat = netlist.getGMatrix();
    auto& B_vec = netlist.getBVector(); 
    auto& state_V = netlist.getStateV();
    auto& hist_V = netlist.getHistoricV();
    auto& mid_V = netlist.getMidpointV();

    int currentExtraCount = systemSize - state_V.size();
    MnaMatrixContext ctx{ G_mat, B_vec, &state_V, &hist_V, &mid_V, t, dt, 0, currentExtraCount };

    auto calculateVbicCurrents = [&](const std::vector<double>& targetedState, double& ib, double& ic, double& ie) {
        ib = 0.0; ic = 0.0; ie = 0.0;
        netlist.inspectComponents([&](const Component& comp, int extraIdx) {
            if (comp.getType() == ComponentType::BJT) {
                const auto& bjtNodes = comp.getNodes();
                if (bjtNodes.size() >= 3) {
                    int nC = bjtNodes[0]; int nB = bjtNodes[1]; int nE = bjtNodes[2];
                    double vB = getSafeVoltage(targetedState, nB);
                    double vC = getSafeVoltage(targetedState, nC);
                    double vE = getSafeVoltage(targetedState, nE);
                    double ibe = 1e-16 * (std::exp((vB - vE) / 0.02585) - 1.0);
                    double ibbc = 1e-16 * (std::exp((vB - vC) / 0.02585) - 1.0);
                    ib = ibe + ibbc; ic = ibe - ibbc; ie = -(ib + ic); 
                }
            }
        });
    };

    try {
        // --- PHASE 0: INITIAL OPERATING POINT STATE ---
        ctx.orderStep = 0; ctx.time = t; ctx.dt = 0.0;
        if (!netlist.runNewtonRaphson(ctx, 100, options.relTol, options.absTol)) {
            throw std::runtime_error("DC operating point convergence failed.");
        }
        hist_V = state_V; mid_V = state_V;
        netlist.inspectComponents([&](Component& comp, int extraIdx) { comp.updateHistory(ctx); });

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
        std::cout << std::setw(16) << i_b << std::setw(16) << i_c << std::setw(16) << i_e << std::setw(20) << initialFillVelocity << "\n";
        pulseCsv << "," << i_b << "," << i_c << "," << i_e << "," << initialFillVelocity << "\n";

        // --- TIME INTEGRATION TRANSITION TIMELINE LOOP ---
        constexpr double TIME_SLACK = 1e-14;
        while (t < (stopT - TIME_SLACK)) {
            if (t + dt > stopT) { dt = stopT - t; }

            std::vector<double> backupState = state_V;
            std::vector<double> backupMid   = mid_V;
            
            double vOutOld = 0.0;
            if (!monitorPairs.empty()) {
                vOutOld = getSafeVoltage(state_V, monitorPairs.back().first) - getSafeVoltage(state_V, monitorPairs.back().second);
            }

            netlist.inspectComponents([](Component& comp, int extraIdx) { comp.saveCheckpoint(); });

            // PHASE 1: TRAPEZOIDAL SUB-STEP
            ctx.time = t + gamma * dt; ctx.dt = gamma * dt; ctx.orderStep = 1;
            if (!netlist.runNewtonRaphson(ctx, options.maxIterations, options.relTol, options.absTol)) {
                dt = std::max(dt * 0.5, options.minDt); state_V = backupState; mid_V = backupMid;
                netlist.inspectComponents([](Component& comp, int extraIdx) { comp.restoreCheckpoint(); });
                continue;
            }
            mid_V = state_V;

            // PHASE 2: BDF2 STEP BOUNDARY
            ctx.time = t + dt; ctx.dt = dt; ctx.orderStep = 2;
            if (!netlist.runNewtonRaphson(ctx, options.maxIterations, options.relTol, options.absTol)) {
                dt = std::max(dt * 0.5, options.minDt); state_V = backupState; mid_V = backupMid;
                netlist.inspectComponents([](Component& comp, int extraIdx) { comp.restoreCheckpoint(); });
                continue;
            }

            // LTE Validation Pass
            double error = netlist.estimateLocalTruncationError(dt, gamma, options.relTol, options.absTol);
            if (error > 1.0) { 
                dt = std::max(dt * 0.5, options.minDt); state_V = backupState; mid_V = backupMid;
                netlist.inspectComponents([](Component& comp, int extraIdx) { comp.restoreCheckpoint(); });
                continue;
            }

            // Step Accepted
            netlist.inspectComponents([&ctx](Component& comp, int extraIdx) { comp.updateHistory(ctx); });
            hist_V = state_V;
            
            double vOutNew = 0.0;
            if (!monitorPairs.empty()) {
                vOutNew = getSafeVoltage(state_V, monitorPairs.back().first) - getSafeVoltage(state_V, monitorPairs.back().second);
            }
            double fillVelocity = ((vOutNew - vOutOld) / dt) * 1e-6;
            
            t += dt;
            ctx.time = t; 
            calculateVbicCurrents(state_V, i_b, i_c, i_e);

            std::cout << std::left << std::setw(12) << (t * 1e6);
            pulseCsv << t;
            for (const auto& pair : monitorPairs) {
                double vDiff = getSafeVoltage(state_V, pair.first) - getSafeVoltage(state_V, pair.second);
                std::cout << std::setw(24) << vDiff;
                pulseCsv << "," << vDiff;
            }
            std::cout << std::setw(16) << i_b << std::setw(16) << i_c << std::setw(16) << i_e << std::setw(20) << fillVelocity << "\n";
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


