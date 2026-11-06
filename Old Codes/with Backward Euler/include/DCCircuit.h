#pragma once

#include "CircuitComponents.h"
#include "math.h"
#include <vector>
#include <memory>
#include <map>
#include <utility> 

class DCCircuit {
private:
    std::vector<std::shared_ptr<Component>> components;
    std::map<std::shared_ptr<Component>, int> componentExtraIndices;
    int numNodes = 0;
    int totalSize = 0;
    std::shared_ptr<VoltageSource> sweepSource = nullptr;
    std::vector<std::pair<int, int>> sparsityPatterns;

    bool solveMna(const std::vector<double>& oldV, std::vector<double>& outV);

public:
    DCCircuit(const std::vector<std::shared_ptr<Component>>& parsedComponents, int nodeCount);

    bool isDcCircuit() const;

    void runOperatingPointAnalysis();

    void runDcSweepAnalysis(double startV = 0.0, double endV = 12.0, double stepV = 0.5);
};
