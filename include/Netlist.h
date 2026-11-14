
// Netlist.h
#pragma once
#include "CircuitComponents.h"
#include "math.h"
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

struct NetlistSimulationOptions {
    double startTime = 0.0;
    double stopTime = 1e-3;
    double initialDt = 1e-6;
    double minDt = 1e-9;
    double maxDt = 1e-5;
    double relTol = 1e-3;
    double absTol = 1e-6;
    int maxIterations = 20;
};

class Netlist {
private:
    std::vector<std::shared_ptr<Component>> components;
    int numNodes = 0;
    int numExtraRows = 0;
    int totalSize = 0;

    CsrMatrix G;
    std::vector<double> B;
    std::vector<double> stateV;     // Solved state vector v_{n+1}
    std::vector<double> historicV;  // Converged past step state vector v_n
    std::vector<double> midpointV;  // Converged TR midpoint step state vector v_{n+gamma}

    void analyzeTopology();

    std::map<NetNode, std::string> nodeNameMap;

public:
    Netlist() = default;

    // Type-safe NetNode label registration and lookup
    void setNodeName(NetNode node, const std::string& name) {
        if (!node.isGround()) {
            nodeNameMap[node] = name;
        }
    }

    std::string getNodeName(NetNode node) const {
        if (node.isGround()) return "GND";
        auto it = nodeNameMap.find(node);
        if (it != nodeNameMap.end()) {
            return it->second;
        }
        return "Node_" + std::to_string(node.id);
    }

    void printComponentNodes() const;
    void printComponentInputOutputNodes() const;

    void addComponent(std::shared_ptr<Component> component);
    void buildSparsityPattern();
    bool runNewtonRaphson(MnaMatrixContext& ctx, int maxIter, double relTol, double absTol);
    double estimateLocalTruncationError(double dt, double gamma, double relTol, double absTol);
    void simulate(const NetlistSimulationOptions& options);

    // Public Accessor Infrastructure for Solver Bridges
    const std::vector<double>& getCurrentState() const { return stateV; }
    int getTotalSystemSize() const { return totalSize; }
    int getNumNodes() const { return numNodes; }

    const CsrMatrix& getGMatrix() const { return G; }
    CsrMatrix& getGMatrix() { return G; }

    // References to vectors for pointer-based context mapping
    std::vector<double>& getStateV() { return stateV; }
    std::vector<double>& getHistoricV() { return historicV; }
    std::vector<double>& getMidpointV() { return midpointV; }
    std::vector<double>& getBVector() { return B; }

    // Type-Safe Topology Inspection Interfaces (passing MnaIndex for extra rows)
    void inspectComponents(std::function<void(const Component&, MnaIndex extraIdx)> callback) const;

    void inspectComponents(std::function<void(Component&, MnaIndex extraIdx)> callback) {
        for (const auto& comp : components) {
            if (comp) {
                callback(*comp, comp->getAssignedExtraIdx());
            }
        }
    }
};

