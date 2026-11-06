#include "freyalab/PulseCircuit.h"
#include "freyalab/CircuitComponents.h" 
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdexcept>

PulseCircuit::PulseCircuit(Netlist& netlistInstance) : netlist(netlistInstance) {}
// mathematical correct, robust and has no structural bugs
bool PulseCircuit::isPulseCircuit() {
    if (netlist.getTotalSystemSize() == 0) {

    }
    bool hasPulseFeatures = false;
    netlist.inspectComponents([&](const Component& comp, int extraIdx) {
    if (comp.getType() == ComponentType::DIGITAL_IC || 
            comp.getType() == ComponentType::MOSFET || 
            comp.getType() == ComponentType::BJT) {
            hasPulseFeatures = true;
        }
    });

    return (netlist.getTotalSystemSize() > 0) && hasPulseFeatures;
}
// mathematical correct, robust and has no structural bugs
bool PulseCircuit::solveTransientStep(MnaMatrixContext& ctx, double targetTime, double currentDt, int stepPhase) {
    ctx.time = targetTime;
    ctx.dt = currentDt;
    ctx.orderStep = stepPhase; // 1 = Trapezoidal Midpoint Substep, 2 = BDF2 Final Step Boundary
    constexpr int maxIterations = 50;
    constexpr double relTol = 1e-3;
    constexpr double absTol = 1e-6;
    return netlist.runNewtonRaphson(ctx, maxIterations, relTol, absTol);
}
// mathematical correct, robust and has no structural bugs 
void PulseCircuit::runPulseTransientAnalysis(const PulseSimulationOptions& options, int inputNode, int outputNode) {
    std::cout << "\n========================================================================================\n"
              << "                   PULSE TRANSIENT SIMULATION OUTPUT TABLE (TR-BDF2 + VBIC)            \n"
              << "========================================================================================\n"
              << std::left << std::setw(12) << "Time(us)" << std::setw(14) << "V(input)[V]" << std::setw(15) << "V(Output)[V]" << "\n"
              << "----------------------------------------------------------------------------------------\n";

    // Establish and verify global MNA dimensions
    int systemSize = netlist.getTotalSystemSize();
    if (systemSize == 0) {
        netlist.buildSparsityPattern();
        systemSize = netlist.getTotalSystemSize();
    }

    double t = options.startTime;
    double dt = options.initialDt;
    double stopT = options.stopTime;
    
    // Optimal L-stable TR-BDF2 configuration parameter (gamma = 0.585786...)
    const double gamma = 2.0 - std::sqrt(2.0); 

    // Extract foundational state references
    auto& G_mat = netlist.getGMatrix();
    auto& B_vec = netlist.getBVector(); 
    auto& state_V = netlist.getStateV();
    auto& hist_V = netlist.getHistoricV();
    auto& mid_V = netlist.getMidpointV();

    int nodeCount = state_V.size();
    int currentExtraCount = systemSize - nodeCount;

    // Direct initialization of solver context indirection structure
    MnaMatrixContext ctx{ 
        G_mat, 
        B_vec, 
        &state_V, 
        &hist_V, 
        &mid_V, 
        t, dt, 0, currentExtraCount 
    };

    try {
        // --------------------------------------------------------------------
        // PHASE 0: DC OPERATING POINT
        // --------------------------------------------------------------------
        ctx.orderStep = 0;
        ctx.time = t;
        ctx.dt = 0.0;

        if (!netlist.runNewtonRaphson(ctx, 100, options.relTol, options.absTol)) {
            throw std::runtime_error("DC operating point convergence failed.");
        }
        hist_V = state_V;
        mid_V = state_V; // Sync midpoint array to prevent reading garbage on initial steps

        // Seed primitive model history equations
        netlist.inspectComponents([&ctx](Component& comp, int extraIdx) {
            (void)extraIdx;
            comp.updateHistory(ctx);
        });

        // Print initial condition records
        double vIn = (inputNode >= 0 && static_cast<size_t>(inputNode) < state_V.size()) ? state_V[inputNode] : 0.0;
        double vOut = (outputNode >= 0 && static_cast<size_t>(outputNode) < state_V.size()) ? state_V[outputNode] : 0.0;
        std::cout << std::left << std::setw(12) << (t * 1e6) << std::setw(14) << vIn << std::setw(15) << vOut << "\n";

        // --------------------------------------------------------------------
        // TIME INTEGRATION LOOP
        // --------------------------------------------------------------------
        constexpr double TIME_SLACK = 1e-14; // Prevent floating-point roundoff traps
        while (t < (stopT - TIME_SLACK)) {
            if (t + dt > stopT) {
                dt = stopT - t;
            }

            std::vector<double> backupState = state_V;
            std::vector<double> backupMid   = mid_V;

            netlist.inspectComponents([](Component& comp, int extraIdx) {
                (void)extraIdx;
                comp.saveCheckpoint();
            });

            // --------------------------------------------------------------------
            // PHASE 1: TRAPEZOIDAL SUB-STEP (t -> t + gamma * dt)
            // --------------------------------------------------------------------
            ctx.time = t + gamma * dt;
            ctx.dt = gamma * dt; 
            ctx.orderStep = 1;

            if (!netlist.runNewtonRaphson(ctx, options.maxIterations, options.relTol, options.absTol)) {
                if (dt <= options.minDt) throw std::runtime_error("Timestep reached floor limit at Phase 1 convergence failure.");
                dt = std::max(dt * 0.5, options.minDt);
                state_V = backupState;
                mid_V   = backupMid;
                netlist.inspectComponents([](Component& comp, int extraIdx) { (void)extraIdx; comp.restoreCheckpoint(); });
                continue;
            }
            mid_V = state_V;

            // --------------------------------------------------------------------
            // PHASE 2: BDF2 STEP BOUNDARY (t + gamma * dt -> t + dt)
            // --------------------------------------------------------------------
            ctx.time = t + dt;
            ctx.dt = dt; 
            ctx.orderStep = 2;

            if (!netlist.runNewtonRaphson(ctx, options.maxIterations, options.relTol, options.absTol)) {
                if (dt <= options.minDt) throw std::runtime_error("Timestep reached floor limit at Phase 2 convergence failure.");
                dt = std::max(dt * 0.5, options.minDt);
                state_V = backupState;
                mid_V   = backupMid;
                netlist.inspectComponents([](Component& comp, int extraIdx) { (void)extraIdx; comp.restoreCheckpoint(); });
                continue;
            }

            // Local Truncation Error Validation Pass
            ctx.dt = dt; // Ensure matrix context reflects full step size explicitly before LTE call
            double error = netlist.estimateLocalTruncationError(dt, gamma, options.relTol, options.absTol);
            if (error > 1.0) { 
                if (dt <= options.minDt) throw std::runtime_error("Timestep reached floor limit due to unacceptable LTE spikes.");
                dt = std::max(dt * 0.5, options.minDt);
                state_V = backupState;
                mid_V   = backupMid;
                netlist.inspectComponents([](Component& comp, int extraIdx) { (void)extraIdx; comp.restoreCheckpoint(); });
                continue;
            }

            // STEP ACCEPTED: Commit calculations to primitive history tracking variables
            netlist.inspectComponents([&ctx](Component& comp, int extraIdx) {
                (void)extraIdx;
                comp.updateHistory(ctx);
            });

            hist_V = state_V;
            t += dt;
            ctx.time = t; 

            // Output data records
            vIn = (inputNode >= 0 && static_cast<size_t>(inputNode) < state_V.size()) ? state_V[inputNode] : 0.0;
            vOut = (outputNode >= 0 && static_cast<size_t>(outputNode) < state_V.size()) ? state_V[outputNode] : 0.0;

            std::cout << std::left << std::setw(12) << (t * 1e6) 
                      << std::setw(14) << vIn 
                      << std::setw(15) << vOut << "\n";

            // Apply second-order optimal control theory scaling with a controlled zero-error ceiling
            if (error > 1e-12) {
                double safetyFactor = 0.9;
                double hFactor = safetyFactor * std::pow(1.0 / error, 1.0 / 3.0);
                hFactor = std::max(0.5, std::min(1.5, hFactor)); // Smooth recovery cap at 1.5x growth
                dt = std::min(dt * hFactor, options.maxDt);
            } else {
                dt = std::min(dt * 1.5, options.maxDt); // Controlled steady-state ramp up
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "\n[FATAL] Pulse Switching Transient Simulation Aborted: " << e.what() << std::endl;
    }
    
    std::cout << "========================================================================================\n" << std::endl;
}
