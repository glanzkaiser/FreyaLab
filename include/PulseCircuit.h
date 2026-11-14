// PulseCircuit.h
#pragma once
#include "Netlist.h"
#include "CircuitComponents.h"
#include <string>
#include <utility>
#include <vector>

struct PulseSimulationOptions {
    double startTime = 0.0;
    double stopTime = 200e-6;   // 200 us default time window
    double initialDt = 10e-9;   // 10 ns initial time step size
    double minDt = 1e-12;       // 1 ps timestep floor limit
    double maxDt = 1e-6;        // 1 us timestep ceiling capping
    double relTol = 1e-3;
    double absTol = 1e-6;
    int maxIterations = 20;
};

class PulseCircuit {
private:
    Netlist& netlist;

public:
    explicit PulseCircuit(Netlist& netlistInstance);

    void detectAndPrintCauseEffectNodes(std::vector<NetNode>& causevector,
                                        std::vector<NetNode>& effectvector);

    // Automatically scans the circuit topology, extracts Cause and Effect NetNodes,
    // and maps them as differential measurement pairs relative to Ground (NetNode{0}).
    void autoGenerateMonitorPairs(std::vector<std::pair<NetNode, NetNode>>& monitorPairs);

    // Dynamic inspection to verify if the model space contains valid active elements
    bool isPulseCircuit();

    // Resolves a single standalone nonlinear transient step handling dual-phase TR-BDF2 updates
    bool solveTransientStep(MnaMatrixContext& ctx, double targetTime, double currentDt, int stepPhase);

    // Generates a comprehensive transient metric table and logs data into a CSV file.
    // Evaluates custom differential NetNode pairs dynamically.
    void runPulseTransientAnalysis(const PulseSimulationOptions& options,
                                   const std::vector<std::pair<NetNode, NetNode>>& monitorPairs);
};

