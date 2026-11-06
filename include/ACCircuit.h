// ACCircuit.h
#pragma once
#include "Netlist.h"
#include <string>
#include <vector>
#include <complex>

struct AcSimulationOptions {
    double startTime = 0.0;
    double stopTime = 5e-3;
    double initialDt = 1e-6;
    double minDt = 1e-12;       // ADDED: Floor simulation step limit
    double maxDt = 1e-4;        // ADDED: Ceiling simulation step limit
    double relTol = 1e-3;
    double absTol = 1e-6;
    int maxIterations = 50;     // ADDED: Non-linear iteration limits
};

class ACCircuit {
private:
    Netlist& netlist;

public:
    explicit ACCircuit(Netlist& netlistInstance);

    // Dynamic component inspection to determine AC capability
    bool isAcCircuit();

    // Integrates a single dual-phase non-linear time slice via Netlist orchestrator overrides
    bool solveTransientStep(MnaMatrixContext& ctx, double targetTime, double currentDt, int stepPhase);

    // Performs complete time integration, logging explicit physical values over custom time bounds
    void runTransientWaveformAnalysis(const AcSimulationOptions& options, int inputNode, int outputNode);

    // Generalized Complex Phasor-Domain Small-Signal Linear Multi-Frequency Solvers
    void runFrequencySweepAnalysis(double startFreq, double stopFreq, int pointsCount);
};

