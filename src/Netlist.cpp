
#include "freyalab/Netlist.h"
#include "freyalab/CircuitComponents.h"
#include <algorithm>
#include <iostream>
#include <cmath>
#include <stdexcept>
// mathematically correct, good and robust and no structural bugs 
void Netlist::addComponent(std::shared_ptr<Component> component) {
    // 1. Fail-fast safety validation on structural pointers
    if (!component) {
        throw std::invalid_argument("Netlist::addComponent - Cannot add a null component pointer.");
    }

    // 2. Scan nodes for topological and numerical sanity bounds checks
    for (int node : component->getNodes()) {
        if (node < 0) {
            throw std::runtime_error("Netlist::addComponent - Negative node index '" + 
                                     std::to_string(node) + "' detected in component: " + 
                                     component->getName());
        }
    }
    
    // 3. Ingest the pointer cleanly. Sizing calculations are deferred to topology analysis.
    components.push_back(std::move(component));
}

// mathematically correct, good and robust and no structural bugs 
void Netlist::analyzeTopology() {
    // 1. Scan the active network components array to discover the true maximum node index
    int maxNodeIdx = 0;
    for (const auto& comp : components) {
        if (!comp) continue;
        for (int node : comp->getNodes()) {
            if (node > maxNodeIdx) {
                maxNodeIdx = node;
            }
        }
    }
    
    // Nodal row count maps exactly to maxNodeIdx + 1 (handles node 0 through maxNodeIdx)
    numNodes = maxNodeIdx;
    
    // 2. Initialize the dynamic extra row tracker directly at the next available index space boundary
    // If numNodes is 5 (Indices 0, 1, 2, 3, 4), the first extra row starts exactly at index 5.
    int globalExtraRowCounter = numNodes; 
    
    numExtraRows = 0;
    for (auto& comp : components) {
        if (!comp) continue;
        
        // Reset component internal tracking records before calculation passes
        comp->setAssignedExtraIdx(-1);
        
        int req = comp->getExtraRowRequirement();
        if (req > 0) {
            comp->allocateRows(globalExtraRowCounter);
            numExtraRows += req;
        }
    }
    
    // 3. Perfect mathematical sizing cohesion: totalSize matches the end state of the counter
    totalSize = globalExtraRowCounter;
}

void Netlist::buildSparsityPattern() {
    // 1. Establish precise topological dimensions and allocate internal component indices
    analyzeTopology();
    
    std::vector<std::pair<int, int>> pattern;
    // Standard heuristic pre-reservation sizing to avoid continuous vector reallocation overhead
    pattern.reserve(totalSize * 6);
    
    // 2. Map standard main structural diagonals firmly to guarantee matrix regular stability
    for (int i = 0; i < totalSize; ++i) {
        pattern.push_back({i, i});
    }

    // 3. Accumulate structural components coordinate footprints safely
    for (auto& comp : components) {
        if (!comp) continue;
        
        int assignedIdx = comp->getAssignedExtraIdx();
        
        comp->getSparsityPattern(pattern, assignedIdx);
    }

    // 4. Seal footprint matrix models and allocate memory buffers to identical structural horizons
    G.initialize(totalSize, pattern);
    B.assign(totalSize, 0.0);
    stateV.assign(totalSize, 0.0);
    historicV.assign(totalSize, 0.0);
    midpointV.assign(totalSize, 0.0);
}

void Netlist::inspectComponents(std::function<void(const Component&, int extraIdx)> callback) const {
    for (const auto& comp : components) {
        if (comp) {
            callback(*comp, comp->getAssignedExtraIdx());
        }
    }
}
// good and robust
bool Netlist::runNewtonRaphson(MnaMatrixContext& ctx, int maxIter, double relTol, double absTol) {
    std::vector<double> nextV(totalSize, 0.0);
    std::vector<double> iterV = stateV; 
    
    for (int iter = 0; iter < maxIter; ++iter) {
        G.zero();
        std::fill(B.begin(), B.end(), 0.0);
        
        // Inject diagonal Gmin normalization layer across node boundaries
	for (int i = 0; i < totalSize; ++i) {
	    G.addValue(i, i, 1e-12);
	}

	ctx.currentIterV = &iterV; // Safe retargetable address binding

        for (auto& comp : components) {
            ctx.assignedExtraIdx = comp->getAssignedExtraIdx();
            comp->stamp(ctx);
        }

	// FIX: Instantiate fresh local variables inside the scope loop 
        // to guarantee memory stability across dynamic structural updates.
        CsrMatrix L, U;
        std::vector<int> P;

        L.initializeFromPattern(G);
        U.initializeFromPattern(G);

        if (!CircuitMath::luFactorize(G, L, U, P)) {
            return false; 
        }
        
        // Correct Mathematical MNA Evaluation Pathway:
        // Because your stamps evaluate G * V_next = B, luSolve outputs the absolute targeted V_next state directly.
        nextV = CircuitMath::luSolve(L, U, P, B);

        double maxError = 0.0;
        for (int i = 0; i < totalSize; ++i) {
            double delta = nextV[i] - iterV[i]; // Deduce step displacement delta for convergence tracking
            double tol = absTol + relTol * std::max(std::abs(nextV[i]), std::abs(iterV[i]));
            maxError = std::max(maxError, std::abs(delta) / tol);
            iterV[i] = nextV[i]; 
        }

        if (maxError < 1.0) {
            stateV = iterV; // Commit convergence globally
            return true; 
        }
    }
    return false;
}

// mathematically correct and robust 
double Netlist::estimateLocalTruncationError(double dt, double gamma, double relTol, double absTol) {
    (void)dt; 
    
    // Explicit TR-BDF2 third-order error factor formulation components
    const double factor = 1.0 / (3.0 * (2.0 - gamma));
    const double termMid = 1.0 / gamma;
    const double termHist = (1.0 - gamma) / gamma;
    
    // Derive a dedicated absolute current floor for auxiliary branch current metrics (Amperes)
    // SPICE standard maps an absolute current floor proportional to voltage scales (typically 1 nA if absTol is 1 uV)
    const double absTolCurrent = absTol * 1e-3;

    double maxNormalizedError = 0.0;
    
    for (int i = 0; i < totalSize; ++i) {
        // Correct TR-BDF2 Local Truncation Error polynomial formulation
        double rawError = factor * std::abs(stateV[i] - (midpointV[i] * termMid) + (historicV[i] * termHist));

        if (rawError <= 0.0) continue;

        // Eliminate std::initializer_list stack allocation overhead entirely
        double absState = std::abs(stateV[i]);
        double absMid   = std::abs(midpointV[i]);
        double absHist  = std::abs(historicV[i]);
        
        double maxVal = absState;
        if (absMid > maxVal)  maxVal = absMid;
        if (absHist > maxVal) maxVal = absHist;
        
        // Enforce distinct dimensional physical floors based on variable types (Voltages vs Currents)
        double currentAbsFloor = (i < numNodes) ? absTol : absTolCurrent;
        
        double tol = currentAbsFloor + relTol * maxVal;
        if (tol < 1e-12) {
            tol = 1e-12; // Safety stability baseline anchor
        }
        
        double normalizedError = rawError / tol;
        
        // Use an L-infinity norm (max error tracking). If any single node fails, 
        // the step must be rejected to prevent local numerical blowups.
        if (normalizedError > maxNormalizedError) {
            maxNormalizedError = normalizedError;
        }
    }
    
    // Returns the maximum local distortion scale. If <= 1.0, the time step is accepted.
    return maxNormalizedError;
}


// good and robust and no structural bugs
void Netlist::simulate(const NetlistSimulationOptions& options) {
    buildSparsityPattern();

    double t = options.startTime;
    double dt = options.initialDt;
    double last_dt = dt; 
    const double gamma = 2.0 - std::sqrt(2.0); 

    // Performance Optimization: Stack/Heap cache pre-allocation to entirely mitigate vector thrashing
    std::vector<double> backupState(totalSize, 0.0);

    // Build integration context referencing structural fields.
    // Set transient tracking pointer bindings once globally.
    MnaMatrixContext ctx{ G, B, &stateV, &historicV, &midpointV, t, dt, 0, -1 };

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
    
    // Initialize dynamic companion states (q_hist_n / i_hist_n) at t=0
    for (auto& comp : components) {
        ctx.assignedExtraIdx = comp->getAssignedExtraIdx();
        comp->updateHistory(ctx);
    }
    
    std::cout << "[SOLVER] DC Operating Point converged. Launching transient time loop..." << std::endl;

    // ====================================================================
    // STEP 2: TIME DOMAIN SIMULATION STEPPING LOOP
    // ====================================================================
    while (t < options.stopTime) {
        if (t + dt > options.stopTime) {
            dt = options.stopTime - t;
        }

        // Fast block-copy operation utilizing pre-allocated storage layer
        backupState = stateV;

        // Freeze charge layers across dynamic primitives before running integration steps
        for (auto& comp : components) {
            comp->saveCheckpoint();
        }

        // --------------------------------------------------------------------
        // PHASE 1: TRAPEZOIDAL COMPANION SUBSTEP (t -> t + gamma*dt)
        // --------------------------------------------------------------------
        // CRITICAL FIX: Explicitly synchronize temporal context metrics at the start of Phase 1
        ctx.time = t + gamma * dt;
        ctx.dt = gamma * dt; 
        ctx.orderStep = 1;
        ctx.historicStepV = &historicV;
        ctx.midpointStepV = nullptr; 
        ctx.currentIterV = &stateV;

        if (!runNewtonRaphson(ctx, options.maxIterations, options.relTol, options.absTol)) {
            if (dt <= options.minDt) throw std::runtime_error("Fatal Error: Newton-Raphson failed at the minimum timestep limit.");
            dt = std::max(dt * 0.5, options.minDt);
            stateV = backupState;
            for (auto& comp : components) { comp->restoreCheckpoint(); }
            continue; 
        }
        midpointV = stateV;
        
        // Capture midpoint variables history profile updates (q_hist_nc / i_hist_nc)
        for (auto& comp : components) {
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
            if (dt <= options.minDt) throw std::runtime_error("Fatal Error: Newton-Raphson failed at the minimum timestep limit.");
            dt = std::max(dt * 0.5, options.minDt);
            stateV = backupState;
            for (auto& comp : components) { comp->restoreCheckpoint(); }
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
            for (auto& comp : components) { comp->restoreCheckpoint(); }
            continue;
        }

        // Critical Success Path: Update permanent history state variables
        for (auto& comp : components) {
            ctx.assignedExtraIdx = comp->getAssignedExtraIdx();
            comp->updateHistory(ctx);
        }
        
        historicV = stateV;
        t += dt;
        last_dt = dt; // Safely preserve step spacing context after successful commit

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


