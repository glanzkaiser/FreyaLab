#pragma once

#include "CircuitComponents.h"
#include "math.h"
#include <vector>
#include <memory>
#include <map>
#include <utility> // Guarantee std::pair template instantiation across strict layout boundaries

class ACCircuit {
private:
    std::vector<std::shared_ptr<Component>> components;
    std::map<std::shared_ptr<Component>, int> componentExtraIndices;
    int numNodes = 0;
    int totalSize = 0;
    std::vector<std::pair<int, int>> sparsityPatterns;

    // High-precision iterative state solver mapping Option B context definitions
    bool solveTransientStep(double time, double dt, const std::vector<double>& oldV, std::vector<double>& outV);

public:
    // Dynamic constructor initializing indices and pre-compiling CSR patterns
    ACCircuit(const std::vector<std::shared_ptr<Component>>& parsedComponents, int nodeCount);

    // Structural loop check mapping whether network contains AC components
    bool isAcCircuit() const;

    // Time-domain transient step sweep utilizing Backward Euler companion representations
    void runTransientWaveformAnalysis(double tEnd = 2e-3, double dt = 5e-6);

    // Alternating current quadrature frequency sweep mapping vector outputs to disk
    void runFrequencySweepAnalysis(double startFreq = 10.0, double endFreq = 100000.0, int pointsPerDecade = 10);
};

