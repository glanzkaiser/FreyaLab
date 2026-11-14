// ACCircuit.h
#pragma once
#include "Netlist.h"
#include "CircuitComponents.h"
#include <complex>
#include <string>
#include <vector>

struct AcSimulationOptions {
    double startTime = 0.0;
    double stopTime = 5e-3;
    double initialDt = 1e-6;
    double minDt = 1e-12;       // Floor simulation step limit
    double maxDt = 1e-4;        // Ceiling simulation step limit
    double relTol = 1e-3;
    double absTol = 1e-6;
    int maxIterations = 50;     // Nonlinear iteration limits
};

class ACCircuit {
private:
    Netlist& netlist;

public:
    explicit ACCircuit(Netlist& netlistInstance);

    // Dynamic component inspection to determine AC capability
    bool isAcCircuit();

    // Integrates a single dual-phase nonlinear time slice via Netlist orchestrator
    bool solveTransientStep(MnaMatrixContext& ctx, double targetTime, double currentDt, int stepPhase);

    // Performs complete time integration, logging explicit physical values over type-safe NetNode bounds
    void runTransientWaveformAnalysis(const AcSimulationOptions& options,
                                      NetNode inputNode,
                                      NetNode outputNode);

    // Generalized Complex Phasor-Domain Small-Signal Linear Multi-Frequency Solver
    void runFrequencySweepAnalysis(double startFreq, double stopFreq, int pointsCount);
};

