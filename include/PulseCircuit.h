// PulseCircuit.h
#pragma once
#include "Netlist.h"
#include <string>
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

    // Dynamic inspection to verify if the model space contains valid active elements
    bool isPulseCircuit();

    // Resolves a single standalone non-linear transient step handling dual-phase TR-BDF2 updates
    bool solveTransientStep(MnaMatrixContext& ctx, double targetTime, double currentDt, int stepPhase);

    // Generates a fully generic pulse switching transient table tracking mapped nodes dynamically
    void runPulseTransientAnalysis(const PulseSimulationOptions& options, int inputNode, int outputNode);
};

