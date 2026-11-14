
#include "freyalab/Netlist.h"
#include "freyalab/CircuitComponents.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

void Netlist::addComponent(std::shared_ptr<Component> component) {
    // 1. Fail-fast safety validation on structural pointers
    if (!component) {
        throw std::invalid_argument("Netlist::addComponent - Cannot add a null component pointer.");
    }

    // 2. Scan NetNode entities for topological and numerical sanity bounds checks
    for (NetNode node : component->getNodes()) {
        if (node.id < 0) {
            throw std::runtime_error("Netlist::addComponent - Negative NetNode ID '" +
                                     std::to_string(node.id) + "' detected in component: " +
                                     component->getName());
        }
    }

    // 3. Ingest the pointer cleanly; sizing calculations are deferred to topology analysis
    components.push_back(std::move(component));
}

void Netlist::printComponentInputOutputNodes() const {
    std::cout << "\n====================================================================================\n";
    std::cout << "                     COMPONENTS NODE INPUT / OUTPUT TERMINAL MAP                    \n";
    std::cout << "====================================================================================\n";
    std::cout << std::left << std::setw(15) << "Device"
              << std::setw(18) << "Type"
              << std::setw(22) << "Input Terminal"
              << std::setw(22) << "Output Terminal"
              << "Extra Control Pins\n";
    std::cout << "------------------------------------------------------------------------------------\n";

    this->inspectComponents([](const Component& comp, MnaIndex extraIdx) {
        (void)extraIdx;
        const auto& nodes = comp.getNodes();

        std::cout << std::left << std::setw(15) << comp.getName();

        std::string typeStr = "";
        std::string inStr = "N/A";
        std::string outStr = "N/A";
        std::string extraStr = "";

        switch (comp.getType()) {
            case ComponentType::RESISTOR:
            case ComponentType::CAPACITOR:
            case ComponentType::INDUCTOR:
            case ComponentType::DIODE:
            case ComponentType::VOLTAGE_SOURCE:
            case ComponentType::CURRENT_SOURCE: {
                if (comp.getType() == ComponentType::RESISTOR)       typeStr = "Resistor";
                else if (comp.getType() == ComponentType::CAPACITOR) typeStr = "Capacitor";
                else if (comp.getType() == ComponentType::INDUCTOR)  typeStr = "Inductor";
                else if (comp.getType() == ComponentType::DIODE)     typeStr = "Diode (A->K)";
                else                                                 typeStr = "Source";

                if (nodes.size() >= 2) {
                    inStr  = std::to_string(nodes[0].id) + " [Positive/Anode]";
                    outStr = std::to_string(nodes[1].id) + " [Negative/Cathode]";
                }
                break;
            }
            case ComponentType::MOSFET:
            case ComponentType::FINFET:
            case ComponentType::RIBBONFET:
            case ComponentType::GAA_FET: {
                if (comp.getType() == ComponentType::FINFET)         typeStr = "FinFET (D-G-S)";
                else if (comp.getType() == ComponentType::RIBBONFET) typeStr = "RibbonFET (D-G-S)";
                else if (comp.getType() == ComponentType::GAA_FET)   typeStr = "GAA-FET (D-G-S)";
                else                                                 typeStr = "MOSFET (D-G-S)";

                if (nodes.size() >= 3) {
                    inStr    = std::to_string(nodes[1].id) + " [Gate]";
                    outStr   = std::to_string(nodes[0].id) + " [Drain]";
                    extraStr = "Source: " + std::to_string(nodes[2].id);
                }
                break;
            }
            case ComponentType::BJT: {
                typeStr = "BJT (C-B-E)";
                if (nodes.size() >= 3) {
                    inStr    = std::to_string(nodes[1].id) + " [Base]";
                    outStr   = std::to_string(nodes[0].id) + " [Collector]";
                    extraStr = "Emitter: " + std::to_string(nodes[2].id);
                }
                break;
            }
            case ComponentType::DIGITAL_IC: {
                typeStr = "Logic Gate";
                if (nodes.size() >= 2) {
                    // Matches BehaviorLogicGateComponent layout: [InA, (InB...), Out]
                    outStr = std::to_string(nodes.back().id) + " [Out Pin]";

                    std::string inputs = "";
                    for (size_t i = 0; i + 1 < nodes.size(); ++i) {
                        inputs += std::to_string(nodes[i].id);
                        if (i + 2 < nodes.size()) inputs += ", ";
                    }
                    inStr = inputs.empty() ? "None" : inputs + " [In Pin(s)]";
                }
                break;
            }
            default: {
                typeStr = "Primitive Component";
                if (!nodes.empty()) inStr = std::to_string(nodes[0].id);
                if (nodes.size() >= 2) outStr = std::to_string(nodes[1].id);
                break;
            }
        }

        std::cout << std::left << std::setw(18) << typeStr
                  << std::setw(22) << inStr
                  << std::setw(22) << outStr
                  << extraStr << "\n";
    });

    std::cout << "====================================================================================\n\n";
}

void Netlist::printComponentNodes() const {
    std::cout << "\n==================================================\n";
    std::cout << "          NETLIST COMPONENT NODE MAPPING          \n";
    std::cout << "==================================================\n";
    std::cout << std::left << std::setw(20) << "Component Name" << "Connected Nodes (ID: Label)\n";
    std::cout << "--------------------------------------------------\n";

    // Traverse all instantiated primitives safely using type-safe MnaIndex & NetNode
    this->inspectComponents([this](const Component& comp, MnaIndex extraIdx) {
        (void)extraIdx;

        std::cout << std::left << std::setw(20) << comp.getName();

        const auto& associatedNodes = comp.getNodes();
        if (associatedNodes.empty()) {
            std::cout << "[No nodes connected]";
        } else {
            for (size_t i = 0; i < associatedNodes.size(); ++i) {
                NetNode node = associatedNodes[i];
                std::string nodeLabel = this->getNodeName(node);
                std::cout << node.id << ":" << nodeLabel;
                if (i + 1 < associatedNodes.size()) {
                    std::cout << ", ";
                }
            }
        }
        std::cout << "\n";
    });

    std::cout << "==================================================\n\n";
}

void Netlist::analyzeTopology() {
    // 1. Scan active network components to discover the true maximum 1-based NetNode ID
    int maxNodeId = 0;
    for (const auto& comp : components) {
        if (!comp) continue;
        for (NetNode node : comp->getNodes()) {
            if (!node.isGround() && node.id > maxNodeId) {
                maxNodeId = node.id;
            }
        }
    }

    // External voltage row count equals maxNodeId (1-based NetNode 1..maxNodeId -> 0-based MnaIndex 0..maxNodeId-1)
    numNodes = maxNodeId;

    // 2. Initialize the dynamic extra row counter directly at the next available 0-based MnaIndex boundary
    int globalExtraRowCounter = numNodes;

    numExtraRows = 0;
    for (auto& comp : components) {
        if (!comp) continue;

        // Reset component internal tracking records with an invalid MnaIndex (-1)
        comp->setAssignedExtraIdx(MnaIndex{});

        int req = comp->getExtraRowRequirement();
        if (req > 0) {
            comp->allocateRows(globalExtraRowCounter);
            numExtraRows += req;
        }
    }

    // 3. Total MNA matrix system dimension equals external nodes + allocated extra branch/internal rows
    totalSize = globalExtraRowCounter;
}

void Netlist::buildSparsityPattern() {
    // 1. Establish precise topological dimensions and allocate internal component indices
    analyzeTopology();

    std::vector<std::pair<MnaIndex, MnaIndex>> typedPattern;
    typedPattern.reserve(totalSize * 6);

    // 2. Map standard main structural diagonals firmly to guarantee matrix regular stability
    for (int i = 0; i < totalSize; ++i) {
        MnaIndex diag{ i };
        typedPattern.push_back({diag, diag});
    }

    // 3. Accumulate structural component coordinate footprints using type-safe MnaIndex
    for (auto& comp : components) {
        if (!comp) continue;

        MnaIndex assignedIdx = comp->getAssignedExtraIdx();
        comp->getSparsityPattern(typedPattern, assignedIdx);
    }

    // 4. Unwrap validated MnaIndex pairs into raw 0-based integers at the CsrMatrix boundary
    std::vector<std::pair<int, int>> rawPattern;
    rawPattern.reserve(typedPattern.size());
    for (const auto& entry : typedPattern) {
        if (entry.first.isValid() && entry.second.isValid()) {
            rawPattern.push_back({entry.first.idx, entry.second.idx});
        }
    }

    // 5. Seal footprint matrix models and allocate memory buffers to identical structural horizons
    G.initialize(totalSize, rawPattern);
    B.assign(totalSize, 0.0);
    stateV.assign(totalSize, 0.0);
    historicV.assign(totalSize, 0.0);
    midpointV.assign(totalSize, 0.0);
}

void Netlist::inspectComponents(std::function<void(const Component&, MnaIndex extraIdx)> callback) const {
    for (const auto& comp : components) {
        if (comp) {
            callback(*comp, comp->getAssignedExtraIdx());
        }
    }
}

bool Netlist::runNewtonRaphson(MnaMatrixContext& ctx, int maxIter, double relTol, double absTol) {
    std::vector<double> iterV = stateV;       // Seed guess from current system state vector
    std::vector<double> prevIterV = stateV;   // Prior iteration guess tracker for pnjlim references

    // Allocate matrix workspaces outside the loop to maximize execution pipeline speed
    CsrMatrix L, U;
    std::vector<int> P;

    L.initializeFromPattern(G);
    U.initializeFromPattern(G);

    for (int iter = 0; iter < maxIter; ++iter) {
        G.zero();
        std::fill(B.begin(), B.end(), 0.0);

        P.clear();
        P.resize(totalSize, 0);

        ctx.currentIterV = &iterV;
        ctx.historicStepV = &prevIterV;

        // 1. Device Physics Stamping Pass (passing strongly-typed MnaIndex)
        for (auto& comp : components) {
            if (comp) {
                ctx.assignedExtraIdx = comp->getAssignedExtraIdx();
                comp->stamp(ctx);
            }
        }

        // 2. Apply Gmin normalization across external voltage nodes
        for (int i = 0; i < numNodes; ++i) {
            G.addValue(i, i, 1e-12);
        }

        // 3. Factorize the system matrix using the pre-allocated L and U maps
        if (!CircuitMath::luFactorize(G, L, U, P)) {
            return false; // Matrix structural singularity hit
        }

        // 4. Resolve absolute new system voltage vector via LU back-substitution
        std::vector<double> solvedV = CircuitMath::luSolve(L, U, P, B);
        if (solvedV.size() != static_cast<size_t>(totalSize)) {
            return false;
        }

        // 5. Apply SPICE voltage stepping limits and check convergence
        bool converged = true;
        prevIterV = iterV;

        for (int i = 0; i < totalSize; ++i) {
            double rawDelta = solvedV[i] - prevIterV[i];
            double limitedDelta = rawDelta;

            // Damping logic: restrict node steps from swinging by more than 0.5V per turn
            if (i < numNodes) {
                if (limitedDelta > 0.5)  limitedDelta = 0.5;
                if (limitedDelta < -0.5) limitedDelta = -0.5;
            }

            iterV[i] = prevIterV[i] + limitedDelta;

            double currentAbsTol = (i < numNodes) ? absTol : 1e-3;
            double allowedTolerance = currentAbsTol + relTol * std::max(std::abs(solvedV[i]), std::abs(prevIterV[i]));
            if (std::abs(rawDelta) > allowedTolerance) {
                converged = false;
            }
        }

        if (converged) {
            stateV = iterV; // Commit successful convergence globally
            return true;
        }
    }

    return false; // Loop exceeded max iterations without converging
}

// mathematically correct and robust 
double Netlist::estimateLocalTruncationError(double dt, double gamma, double relTol, double absTol) {
    (void)dt;

    // Explicit TR-BDF2 third-order error factor formulation components
    const double factor   = 1.0 / (3.0 * (2.0 - gamma));
    const double termMid  = 1.0 / gamma;
    const double termHist = (1.0 - gamma) / gamma;

    // Derive a dedicated absolute current floor for auxiliary branch/internal rows
    const double absTolCurrent = absTol * 1e-3;

    double maxNormalizedError = 0.0;

    for (int i = 0; i < totalSize; ++i) {
        MnaIndex row{ i };
        if (!row.isValid()) continue;

        double vState = getSafeExtraValue(stateV, row);
        double vMid   = getSafeExtraValue(midpointV, row);
        double vHist  = getSafeExtraValue(historicV, row);

        // TR-BDF2 Local Truncation Error polynomial formulation
        double rawError = factor * std::abs(vState - (vMid * termMid) + (vHist * termHist));
        if (rawError <= 0.0) continue;

        double maxVal = std::max({std::abs(vState), std::abs(vMid), std::abs(vHist)});

        // Enforce distinct dimensional physical floors based on variable types (Voltages vs Currents/Internal)
        double currentAbsFloor = (row.idx < numNodes) ? absTol : absTolCurrent;

        double tol = currentAbsFloor + relTol * maxVal;
        if (tol < 1e-12) {
            tol = 1e-12;
        }

        double normalizedError = rawError / tol;
        if (normalizedError > maxNormalizedError) {
            maxNormalizedError = normalizedError;
        }
    }

    // Returns the maximum local distortion scale; if <= 1.0, the time step is accepted
    return maxNormalizedError;
}


// good and robust and no structural bugs
void Netlist::simulate(const NetlistSimulationOptions& options) {
    buildSparsityPattern();

    double t = options.startTime;
    double dt = options.initialDt;
    const double gamma = 2.0 - std::sqrt(2.0);

    // Pre-allocate backup state buffer to avoid heap thrashing during step rejections
    std::vector<double> backupState(totalSize, 0.0);

    // Build integration context with type-safe MnaIndex{} (default invalid = -1)
    MnaMatrixContext ctx{ G, B, &stateV, &historicV, &midpointV, t, dt, 0, MnaIndex{} };

    // ====================================================================
    // STEP 1: INITIAL STABLE DC OPERATING POINT
    // ====================================================================
    std::cout << "[SOLVER] Executing Initial Static Operating Point (DC Analysis)..." << std::endl;
    ctx.orderStep = 0;
    ctx.time = t;
    ctx.dt = 0.0;
    ctx.historicStepV = &stateV;
    ctx.midpointStepV = &stateV;
    ctx.currentIterV = &stateV;

    if (!runNewtonRaphson(ctx, 100, options.relTol, options.absTol)) {
        throw std::runtime_error("Critical Solver Error: DC Operating Point calculation diverged.");
    }

    historicV = stateV;
    midpointV = stateV;

    // Initialize dynamic companion states (q_hist_n / i_hist_n / OSDI states) at t=0
    for (auto& comp : components) {
        if (!comp) continue;
        ctx.assignedExtraIdx = comp->getAssignedExtraIdx();
        comp->updateHistory(ctx);
    }

    std::cout << "[SOLVER] DC Operating Point converged. Launching transient time loop..." << std::endl;

    // ====================================================================
    // STEP 2: TIME DOMAIN SIMULATION STEPPING LOOP
    // ====================================================================
    constexpr double TIME_SLACK = 1e-15;
    while (t < (options.stopTime - TIME_SLACK)) {
        if (t + dt > options.stopTime) {
            dt = options.stopTime - t;
        }

        backupState = stateV;

        // Freeze charge/flux/OSDI state layers across dynamic primitives before stepping
        for (auto& comp : components) {
            if (comp) comp->saveCheckpoint();
        }

        // --------------------------------------------------------------------
        // PHASE 1: TRAPEZOIDAL COMPANION SUBSTEP (t -> t + gamma*dt)
        // --------------------------------------------------------------------
        ctx.time = t + gamma * dt;
        ctx.dt = dt;
        ctx.orderStep = 1;
        ctx.historicStepV = &historicV;
        ctx.midpointStepV = nullptr;
        ctx.currentIterV = &stateV;

        if (!runNewtonRaphson(ctx, options.maxIterations, options.relTol, options.absTol)) {
            if (dt <= options.minDt) {
                throw std::runtime_error("Fatal Error: Newton-Raphson failed at the minimum timestep limit.");
            }
            dt = std::max(dt * 0.5, options.minDt);
            stateV = backupState;
            for (auto& comp : components) {
                if (comp) comp->restoreCheckpoint();
            }
            continue;
        }
        midpointV = stateV;

        // Capture midpoint history updates (q_hist_nc / i_hist_nc / OSDI midpoint state)
        for (auto& comp : components) {
            if (!comp) continue;
            ctx.assignedExtraIdx = comp->getAssignedExtraIdx();
            comp->updateHistory(ctx);
        }

        // --------------------------------------------------------------------
        // PHASE 2: BDF2 COMPANION BOUNDARY STEP (t + gamma*dt -> t + dt)
        // --------------------------------------------------------------------
        ctx.time = t + dt;
        ctx.dt = dt;
        ctx.orderStep = 2;
        ctx.midpointStepV = &midpointV;
        ctx.historicStepV = &historicV;

        if (!runNewtonRaphson(ctx, options.maxIterations, options.relTol, options.absTol)) {
            if (dt <= options.minDt) {
                throw std::runtime_error("Fatal Error: Newton-Raphson failed at the minimum timestep limit.");
            }
            dt = std::max(dt * 0.5, options.minDt);
            stateV = backupState;
            for (auto& comp : components) {
                if (comp) comp->restoreCheckpoint();
            }
            continue;
        }

        // --------------------------------------------------------------------
        // LOCAL TRUNCATION ERROR CONTROL CHECK
        // --------------------------------------------------------------------
        double error = estimateLocalTruncationError(dt, gamma, options.relTol, options.absTol);
        double targetError = 1.0;

        if (error > targetError) {
            if (dt <= options.minDt) {
                throw std::runtime_error("Fatal Error: Timestep truncation violation hit hard floor limit: " + std::to_string(dt));
            }
            dt = std::max(dt * 0.5, options.minDt);
            stateV = backupState;
            for (auto& comp : components) {
                if (comp) comp->restoreCheckpoint();
            }
            continue;
        }

        // Commit permanent history state variables
        for (auto& comp : components) {
            if (!comp) continue;
            ctx.assignedExtraIdx = comp->getAssignedExtraIdx();
            comp->updateHistory(ctx);
        }

        historicV = stateV;
        t += dt;

        std::cout << "[STEP] Time: " << t << " s | Step Width (dt): " << dt << " s | Normalized LTE: " << error << std::endl;

        // Smoothly adjust next timestep using an LTE-driven scaling window
        if (error > 0.0) {
            double safeFactor = 0.9 * std::pow(targetError / error, 1.0 / 3.0);
            safeFactor = std::max(0.3, std::min(2.0, safeFactor));
            dt = std::min(options.maxDt, std::max(options.minDt, dt * safeFactor));
        } else {
            dt = std::min(options.maxDt, dt * 2.0);
        }
    }
    std::cout << "[SOLVER] Transient simulation successfully completed." << std::endl;
}



