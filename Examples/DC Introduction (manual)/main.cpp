#include "freyalab/CircuitComponents.h"
#include "freyalab/Netlist.h"
#include "freyalab/DCCircuit.h"
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

// Helper implementation to inspect the assembled CSR matrix during debug views
void printFullLinearSystem(Netlist& netlist) {
    int size = netlist.getTotalSystemSize();
    const CsrMatrix& G = netlist.getGMatrix();
    const std::vector<double>& B = netlist.getBVector();

    // Reconstruct dense layout from CSR representation vectors for debugging visualization
    std::vector<std::vector<double>> denseG(size, std::vector<double>(size, 0.0));
    for (int r = 0; r < size; ++r) {
        int start = G.rowPointers[r];
        int end = G.rowPointers[r + 1];
        for (int p = start; p < end; ++p) {
            int c = G.colIndices[p];
            if (c >= 0 && c < size) {
                denseG[r][c] = G.values[p];
            }
        }
    }

    std::cout << "\n========================================================================================\n";
    std::cout << "                 MNA SYSTEM MATRIX LAYER STATE (INITIAL EXPERIMENTAL STAMP)               \n";
    std::cout << "========================================================================================\n\n";

    std::cout << "      ";
    for (int c = 0; c < size; ++c) {
        std::cout << "  Col " << std::left << std::setw(2) << c << "    ";
    }
    std::cout << "  |   RHS Vector B\n";
    std::cout << std::string(12, '-') << std::string(size * 11, '-') << "---|------------------\n";

    std::cout << std::scientific << std::setprecision(3);
    for (int r = 0; r < size; ++r) {
        double rhsVal = getSafeExtraValue(B, MnaIndex{r});
        std::cout << "Row " << std::left << std::setw(2) << r << ":";
        for (int c = 0; c < size; ++c) {
            std::cout << " " << (denseG[r][c] >= 0 ? " " : "") << denseG[r][c];
        }
        std::cout << "  |  " << (rhsVal >= 0 ? " " : "") << rhsVal << "\n";
    }
    std::cout << "========================================================================================\n";
}

int main() {
    try {
        std::cout << "[INFO]: Initializing MNA Simulation Engine Framework Components Layer...\n";
        Netlist circuitNetlist;

        // ==========================================
        // 1. INDEPENDENT ACTIVE SOURCE EMITTERS
        // ==========================================
        // V1 1 0 DC 5.0
        circuitNetlist.addComponent(std::make_shared<VoltageSource>("V1", std::vector<NetNode>{{1}, {0}}, 5.0));
        // V2 2 0 DC 2.5
        circuitNetlist.addComponent(std::make_shared<VoltageSource>("V2", std::vector<NetNode>{{2}, {0}}, 2.5));

        // ==========================================
        // 2. SIGNAL ROUTING & CURRENT BIAS CORRIDORS
        // ==========================================
        // R1 1 3 10k
        circuitNetlist.addComponent(std::make_shared<Resistor>("R1", std::vector<NetNode>{{1}, {3}}, 10000.0));
        // R2 1 4 1.2k
        circuitNetlist.addComponent(std::make_shared<Resistor>("R2", std::vector<NetNode>{{1}, {4}}, 1200.0));
        // R3 2 3 1k
        circuitNetlist.addComponent(std::make_shared<Resistor>("R3", std::vector<NetNode>{{2}, {3}}, 1000.0));

        // ==========================================
        // 3. ADVANCED ACTIVE VERTICAL BIPOLAR TRANSISTOR
        // ==========================================
        // Q1 4 3 0 npn -> Collector=4, Base=3, Emitter=0
        circuitNetlist.addComponent(std::make_shared<VbicBjtAdvanced>("Q1", std::vector<NetNode>{{4}, {3}, {0}}, true));

        // ==========================================
        // 4. HIGH-FREQUENCY FILTER CIRCUITS
        // ==========================================
        // C1 4 5 4.7u -> (Open-circuit during pure static DC OP operating sweeps)
        circuitNetlist.addComponent(std::make_shared<Capacitor>("C1", std::vector<NetNode>{{4}, {5}}, 4.7e-6));
        // L1 5 0 10m -> (Short-circuit structural row request in DC solver mode)
        circuitNetlist.addComponent(std::make_shared<Inductor>("L1", std::vector<NetNode>{{5}, {0}}, 10e-3));
        // R4 5 0 50
        circuitNetlist.addComponent(std::make_shared<Resistor>("R4", std::vector<NetNode>{{5}, {0}}, 50.0));

        // ==========================================
        // 5. SOLID-STATE CLAMPING PROTECTION RING & MOSFET
        // ==========================================
        // D1 4 1 -> Anode=4, Cathode=1
        circuitNetlist.addComponent(std::make_shared<Diode>("D1", std::vector<NetNode>{{4}, {1}}));
        // M1 5 4 0 nmos -> Drain=5, Gate=4, Source=0
        circuitNetlist.addComponent(std::make_shared<Mosfet>("M1", std::vector<NetNode>{{5}, {4}, {0}}, true));

        // --- STEP 1: Topological Matrix Grid Alignment & Pattern Compilation ---
        circuitNetlist.buildSparsityPattern();

        // --- STEP 2: Execute Prototype Baseline Stamping Pass for Structural Peeking ---
        CsrMatrix& G = circuitNetlist.getGMatrix();
        std::vector<double>& B = circuitNetlist.getBVector();
        std::vector<double>& stateV = circuitNetlist.getStateV();

        G.zero();
        std::fill(B.begin(), B.end(), 0.0);
        std::fill(stateV.begin(), stateV.end(), 0.0);

        MnaMatrixContext initialPeekCtx{
            G,
            B,
            &stateV,
            &circuitNetlist.getHistoricV(),
            &circuitNetlist.getMidpointV(),
            0.0,
            1e-6,
            0,
            MnaIndex{}
        };

        circuitNetlist.inspectComponents([&](Component& comp, MnaIndex extraIdx) {
            initialPeekCtx.assignedExtraIdx = extraIdx;
            comp.stamp(initialPeekCtx);
        });

        // Output structural layout allocations snapshot
        printFullLinearSystem(circuitNetlist);

        // --- STEP 3: Instantiate Core DC Solver Domain Block & Run Convergence Analysis ---
        DCCircuit dcSolver(circuitNetlist);
        dcSolver.runOperatingPointAnalysis();

        // --- STEP 4: Safe Node Voltage Summary via Type-Safe NetNode Helper ---
        const std::vector<double>& finalState = circuitNetlist.getCurrentState();
        std::cout << "Converged Node Voltages:\n";
        for (int n = 1; n <= circuitNetlist.getNumNodes(); ++n) {
            NetNode node{n};
            std::cout << "  V(Node_" << node.id << ") = "
                      << std::fixed << std::setprecision(6)
                      << getSafeVoltage(finalState, node) << " V\n";
        }

    } catch (const std::exception& error) {
        std::cerr << "\n[FATAL SYSTEM RUNTIME CORRUPTION]: " << error.what() << "\n";
        return 1;
    }

    return 0;
}
