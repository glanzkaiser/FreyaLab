#pragma once

#include "CircuitComponents.h"
#include "math.h"
#include <vector>
#include <memory>
#include <map>
#include <utility> // Guarantee std::pair template instantiation across strict layout boundaries

class PulseCircuit {
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
    PulseCircuit(const std::vector<std::shared_ptr<Component>>& parsedComponents, int nodeCount);

    // Structural loop check mapping whether network contains a Pulse source
    bool isPulseCircuit() const;

    // High-speed time-domain transient analysis step sweep utilizing companion representations
    void runPulseTransientAnalysis(double durationMicroseconds = 20.0, double timeStepNanoseconds = 50.0);
};

