// Netlist.h
#pragma once
#include "CircuitComponents.h"
#include "math.h"
#include "KLU.h"
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

struct AnalysisDirective {
    enum class Kind {
        DC_OP,
        DC_SWEEP,
        AC_SWEEP,
        TRANSIENT,
        NOISE
    };

    Kind kind = Kind::DC_OP;

    // .DC sourceName start stop step
    std::string dcSourceName;
    double dcStart = 0.0;
    double dcStop  = 0.0;
    double dcStep  = 0.1;

    // .AC DEC|OCT|LIN points fstart fstop
    std::string acScaleType = "DEC";
    int acPoints = 100;
    double acStartFreq = 10.0;
    double acStopFreq  = 100e6;

    // .TRAN tstep tstop [tstart] [tmax]
    double tranStep  = 1e-6;
    double tranStop  = 10e-3;
    double tranStart = 0.0;
    double tranMaxDt = 0.0;

    // .NOISE V(outPos [, outNeg]) inputSource DEC|OCT|LIN points fstart fstop
    NetNode noiseOutPos{0};
    NetNode noiseOutNeg{0};
    std::string noiseInputSource;
};

class Netlist {
private:
    std::vector<std::shared_ptr<Component>> components;
    std::vector<AnalysisDirective> analysisDirectives;
    int numNodes = 0;
    int numExtraRows = 0;
    int totalSize = 0;

    CsrMatrix G;
    std::vector<double> B;
    std::vector<double> stateV;     // Solved state vector v_{n+1}
    std::vector<double> historicV;  // Converged past step state vector v_n
    std::vector<double> midpointV;  // Converged TR midpoint step state vector v_{n+gamma}

    // Adaptive KLU sparse solver cache (activated when totalSize >= 300)
    KluSymbolicContext kluCache;
    bool useKluBackend = false;

    void analyzeTopology();

    std::map<NetNode, std::string> nodeNameMap;

    double ambientTempCelsius = 27.0; // Standard SPICE nominal temperature (300.15 K)

public:
    Netlist() = default;

    void setTemperatureCelsius(double tempC) {
        ambientTempCelsius = tempC;
    }

    double getTemperatureCelsius() const {
        return ambientTempCelsius;
    }

    double getTemperatureKelvin() const {
        return ambientTempCelsius + 273.15;
    }

    double getThermalVoltage() const {
        constexpr double K_BOLTZMANN = 1.380649e-23;
        constexpr double Q_ELECTRON  = 1.602176634e-19;
        return (K_BOLTZMANN * getTemperatureKelvin()) / Q_ELECTRON;
    }

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

    void addAnalysisDirective(const AnalysisDirective& cmd) {
        analysisDirectives.push_back(cmd);
    }
    const std::vector<AnalysisDirective>& getAnalysisDirectives() const {
        return analysisDirectives;
    }

    // Public Accessor Infrastructure for Solver Bridges
    const std::vector<double>& getCurrentState() const { return stateV; }
    int getTotalSystemSize() const { return totalSize; }
    int getNumNodes() const { return numNodes; }

    const CsrMatrix& getGMatrix() const { return G; }
    CsrMatrix& getGMatrix() { return G; }

    // KLU backend accessors for ACCircuit & FreyaLabVisualizer
    bool isUsingKluBackend() const { return useKluBackend; }
    KluSymbolicContext& getKluCache() { return kluCache; }

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

