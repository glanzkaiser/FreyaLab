#include "freyalab/PulseCircuit.h"
#include "freyalab/math.h"
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <iostream>

PulseCircuit::PulseCircuit(const std::vector<std::shared_ptr<Component>>& parsedComponents, int nodeCount)
    : components(parsedComponents), numNodes(nodeCount), totalSize(nodeCount) {
    int matrixNodeRows = numNodes - 1; // Node 0 (GND) is discarded from the active equation set
    int extraRows = 0;
    
    for (auto& comp : components) {
        int req = comp->getExtraRowRequirement();
        if (req > 0) { 
            componentExtraIndices[comp] = matrixNodeRows + extraRows; 
            extraRows += req; 
        }
    }
    
    totalSize = matrixNodeRows + extraRows;
    for (auto& comp : components) {
        int ext = componentExtraIndices.count(comp) ? componentExtraIndices[comp] : -1;
        comp->getSparsityPattern(sparsityPatterns, ext);
    }
}

bool PulseCircuit::isPulseCircuit() const {
    for (const auto& comp : components) {
        if (comp->getType() == ComponentType::VOLTAGE_SOURCE) {
            if (std::static_pointer_cast<VoltageSource>(comp)->getSourceType() == SourceType::PULSE) return true;
        }
    }
    return false;
}

bool PulseCircuit::solveTransientStep(double time, double dt, const std::vector<double>& oldV, std::vector<double>& outV) {
    int nrIterations = 0; 
    const int maxNrIterations = 100; 
    const double nrTolerance = 1e-6; 
    bool converged = false;
    
    std::vector<double> nrGuessV = oldV;
    if (nrGuessV.size() != static_cast<size_t>(totalSize)) {
        nrGuessV.assign(totalSize, 0.0);
    }
    
    outV.assign(totalSize, 0.0);
    CsrMatrix G; 
    G.numRows = totalSize; 
    MathSolver::finalizeStructure(G, sparsityPatterns);
    
    while (!converged && nrIterations < maxNrIterations) {
        G.clearNumericalValues(); 
        std::vector<double> B(totalSize, 0.0);
        
        for (auto& comp : components) {
            int assignedExtraIdx = componentExtraIndices.count(comp) ? componentExtraIndices.at(comp) : -1;
            // FIX: Perfect Type and Order matching against Option B MnaMatrixContext Signature
            MnaMatrixContext ctx{G, B, nrGuessV, oldV, time, dt, assignedExtraIdx};
            comp->stamp(ctx);
        }
        
        CsrMatrix L, U; 
        std::vector<int> pivot(totalSize, 0);
        if (!MathSolver::sparseLuFactorize(G, L, U, pivot)) {
            return false;
        }
        
        std::vector<double> nrSolV(totalSize, 0.0); 
        MathSolver::sparseLuSolve(L, U, pivot, B, nrSolV);
        
        double maxDelta = 0.0; 
        for (int i = 0; i < totalSize; ++i) {
            maxDelta = std::max(maxDelta, std::abs(nrSolV[i] - nrGuessV[i]));
        }
        
        nrGuessV = nrSolV; 
        nrIterations++; 
        
        if (maxDelta < nrTolerance) { 
            converged = true; 
            outV = nrSolV; 
        }
    }
    
    if (converged) {
        for (auto& comp : components) {
            if (comp->getType() == ComponentType::DIODE) {
                double vDrop = comp->calculateVoltageDrop(outV);
                std::static_pointer_cast<Diode>(comp)->vd_prev_iter = vDrop;
            }
            else if (comp->getType() == ComponentType::BJT) {
                auto bPtr = std::static_pointer_cast<BJT>(comp);
                const auto& nodes = bPtr->getNodes();
                double vC = (nodes[0] > 0) ? outV[nodes[0] - 1] : 0.0;
                double vB = (nodes[1] > 0) ? outV[nodes[1] - 1] : 0.0;
                double vE = (nodes[2] > 0) ? outV[nodes[2] - 1] : 0.0;
                
                bPtr->vbe_prev = bPtr->isNpn ? (vB - vE) : (vE - vB);
                bPtr->vbc_prev = bPtr->isNpn ? (vB - vC) : (vC - vB);
            }
        }
    }
    return converged;
}

void PulseCircuit::runPulseTransientAnalysis(double durationMicroseconds, double timeStepNanoseconds) {
    double tEnd = durationMicroseconds * 1e-6, dt = timeStepNanoseconds * 1e-9, time = 0.0;
    std::shared_ptr<VoltageSource> vSrc = nullptr; 
    std::shared_ptr<Capacitor> c1 = nullptr; 
    std::shared_ptr<Inductor> l1 = nullptr;
    
    for (auto& comp : components) {
        if (comp->getName() == "V1") vSrc = std::static_pointer_cast<VoltageSource>(comp);
        if (comp->getName() == "C1") c1 = std::static_pointer_cast<Capacitor>(comp);
        if (comp->getName() == "L1") l1 = std::static_pointer_cast<Inductor>(comp);
    }
    
    std::ofstream csvFile("pulse_transient_data.csv");
    if (csvFile.is_open()) {
        csvFile << "Time(us),V_input(V),V_Capacitor(V),I_Capacitor(A),V_Inductor(V),I_Inductor(A)\n";
    }
    
    std::cout << "\n========================================================================================\n"
              << "                         PULSE TRANSIENT SIMULATION OUTPUT TABLE                        \n"
              << "========================================================================================\n"
              << std::left << std::setw(12) << "Time(us)" << std::setw(14) << "V(input)[V]" << std::setw(15) << "V(Cap)[V]" << std::setw(15) << "I(Cap)[A]" << std::setw(15) << "V(Ind)[V]" << std::setw(15) << "I(Ind)[A]" << "\n"
              << "----------------------------------------------------------------------------------------\n";
              
    std::vector<double> oldV(totalSize, 0.0), currentV(totalSize, 0.0); 
    
    while (time <= tEnd + 1e-11) {
        if (solveTransientStep(time, dt, oldV, currentV)) {
            double vInVal = vSrc ? vSrc->getCurrentValue(time) : 0.0;
            
            // FIX: Corrected 0-based indexing layout parameters mapping (Node 2 -> index 1, Node 3 -> index 2)
            double vCVal = (numNodes >= 3) ? currentV[1] : 0.0;
            double vLVal = (numNodes >= 4) ? (currentV[1] - currentV[2]) : 0.0;
            
            int cExtra = (c1 && componentExtraIndices.count(c1)) ? componentExtraIndices[c1] : -1;
            double iCVal = c1 ? c1->calculateCurrent(currentV, cExtra) : 0.0;
            
            int lExtra = (l1 && componentExtraIndices.count(l1)) ? componentExtraIndices[l1] : -1;
            double iLVal = l1 ? l1->calculateCurrent(currentV, lExtra) : 0.0;
            
            std::cout << std::left << std::fixed << std::setprecision(3) << std::setw(12) << (time * 1e6) << std::setprecision(4)
                      << std::setw(14) << vInVal << std::setw(15) << vCVal << std::setw(15) << iCVal << std::setw(15) << vLVal << std::scientific << std::setw(15) << iLVal << "\n";
            
            if (csvFile.is_open()) {
                csvFile << (time * 1e6) << "," << vInVal << "," << vCVal << "," << iCVal << "," << vLVal << "," << iLVal << "\n";
            }
            
            oldV = currentV; 
            time += dt;
        } else {
            std::cerr << "[Fatal Error] Pulse transient convergence failed at t = " << (time * 1e6) << " us. Terminating run.\n";
            break;
        }
    }
    if (csvFile.is_open()) csvFile.close();
    std::cout << "========================================================================================\n";
}
