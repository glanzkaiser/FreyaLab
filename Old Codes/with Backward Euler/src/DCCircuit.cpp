#include "freyalab/DCCircuit.h"
#include "freyalab/math.h"
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <iostream>

DCCircuit::DCCircuit(const std::vector<std::shared_ptr<Component>>& parsedComponents, int nodeCount)
    : components(parsedComponents), numNodes(nodeCount), totalSize(nodeCount), sweepSource(nullptr) {
    int matrixNodeRows = numNodes - 1; 
    int extraRows = 0;
    
    for (auto& comp : components) {
        int req = comp->getExtraRowRequirement();
        if (req > 0) {
            componentExtraIndices[comp] = (numNodes - 1) + extraRows; 
            extraRows += req;
        }
        if (comp->getType() == ComponentType::VOLTAGE_SOURCE && comp->getName() == "V1") {
            sweepSource = std::static_pointer_cast<VoltageSource>(comp);
        }
    }
    
    totalSize = matrixNodeRows + extraRows; 
    for (auto& comp : components) {
        int ext = componentExtraIndices.count(comp) ? componentExtraIndices[comp] : -1;
        comp->getSparsityPattern(sparsityPatterns, ext);
    }
}

bool DCCircuit::isDcCircuit() const {
    for (const auto& comp : components) {
        if (comp->getType() == ComponentType::VOLTAGE_SOURCE) {
            if (std::static_pointer_cast<VoltageSource>(comp)->getSourceType() != SourceType::DC) return false;
        }
    }
    return true;
}

bool DCCircuit::solveMna(const std::vector<double>& oldV, std::vector<double>& outV) {
    int nrIterations = 0;
    const int maxNrIterations = 200; 
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
            int assignedExtraIdx = componentExtraIndices.count(comp) ? componentExtraIndices[comp] : -1;
            // Option B Configuration Strategy: 
            // 1st vector parameter points to active NR iteration layout vector (nrGuessV)
            // 2nd vector parameter maps explicitly down to the historic converged state sequence (oldV)
            MnaMatrixContext ctx{G, B, nrGuessV, oldV, 0.0, 0.0, assignedExtraIdx};
            comp->stamp(ctx);
        }
        
        CsrMatrix L, U; 
        std::vector<int> pivot(totalSize, 0);
        if (!MathSolver::sparseLuFactorize(G, L, U, pivot)) return false;
        
        std::vector<double> nrSolV(totalSize, 0.0);
        MathSolver::sparseLuSolve(L, U, pivot, B, nrSolV);
        
        double maxDelta = 0.0;
        for (int i = 0; i < totalSize; ++i) {
            maxDelta = std::max(maxDelta, std::abs(nrSolV[i] - nrGuessV[i]));
        }
        
        for (auto& comp : components) {
            if (comp->getType() == ComponentType::DIODE) {
                double vDrop = comp->calculateVoltageDrop(nrSolV);
                std::static_pointer_cast<Diode>(comp)->vd_prev_iter = vDrop;
            }
            else if (comp->getType() == ComponentType::BJT) {
                auto bPtr = std::static_pointer_cast<BJT>(comp);
                const auto& nodes = bPtr->getNodes();
                double vC = (nodes[0] > 0) ? nrSolV[nodes[0] - 1] : 0.0;
                double vB = (nodes[1] > 0) ? nrSolV[nodes[1] - 1] : 0.0;
                double vE = (nodes[2] > 0) ? nrSolV[nodes[2] - 1] : 0.0;
                
                bPtr->vbe_prev = bPtr->isNpn ? (vB - vE) : (vE - vB);
                bPtr->vbc_prev = bPtr->isNpn ? (vB - vC) : (vC - vB);
            }
        }

        nrGuessV = nrSolV;
        nrIterations++;
        
        if (maxDelta < nrTolerance) { 
            converged = true; 
            outV = nrSolV; 
        }
    }
    return converged;
}

void DCCircuit::runOperatingPointAnalysis() {
    std::vector<double> oldV(totalSize, 0.0), currentV(totalSize, 0.0);
    if (!solveMna(oldV, currentV)) {
        std::cerr << "[Error] DC Operating Point failed to converge.\n";
        return;
    }
    
    std::cout << "\n========================================================================\n"
              << "                      DC OPERATING POINT TABLE                          \n"
              << "========================================================================\n"
              << std::left << std::setw(10) << "Name" << std::setw(15) << "Nodes" << std::setw(15) << "V-Drop(V)" << std::setw(15) << "Current(A)" << std::setw(15) << "Power(W)" << "\n"
              << "------------------------------------------------------------------------\n";
              
    for (const auto& comp : components) {
        int extraIdx = componentExtraIndices.count(comp) ? componentExtraIndices.at(comp) : -1;
        
        std::string nStr = "";
        const auto& compNodes = comp->getNodes();
        if (compNodes.size() >= 3) { 
            nStr = std::to_string(compNodes[0]) + "-" + std::to_string(compNodes[1]) + "-" + std::to_string(compNodes[2]);
        } else if (compNodes.size() == 2) {
            nStr = std::to_string(compNodes[0]) + "-" + std::to_string(compNodes[1]);
        } else if (!compNodes.empty()) {
            nStr = std::to_string(compNodes[0]);
        }
        
        std::cout << std::left << std::setw(10) << comp->getName() << std::setw(15) << nStr << std::fixed << std::setprecision(4)
                  << std::setw(15) << comp->calculateVoltageDrop(currentV) 
                  << std::setw(15) << comp->calculateCurrent(currentV, extraIdx) 
                  << std::setw(15) << comp->calculatePower(currentV, extraIdx) << "\n";
    }
    std::cout << "========================================================================\n";
}

void DCCircuit::runDcSweepAnalysis(double startV, double endV, double stepV) {
    if (!sweepSource) return;
    std::ofstream csvV("dc_sweep_voltage.csv"), csvI("dc_sweep_current.csv"), csvP("dc_sweep_power.csv");
    
    csvV << "V1(V)"; 
    for (int i = 1; i < numNodes; ++i) csvV << ",Node_" << i; 
    csvV << "\n";
    
    csvI << "V1(V)"; csvP << "V1(V)"; 
    for (const auto& comp : components) { 
        csvI << "," << comp->getName(); 
        csvP << "," << comp->getName(); 
    } 
    csvI << "\n"; csvP << "\n";
    
    std::vector<double> currentV(totalSize, 0.0);
    for (double vVal = startV; vVal <= endV + 1e-9; vVal += stepV) {
        sweepSource->setDCValue(vVal); 
        std::vector<double> sweepV(totalSize, 0.0);
        
        for (auto& comp : components) {
            if (comp->getType() == ComponentType::DIODE) std::static_pointer_cast<Diode>(comp)->resetState();
            if (comp->getType() == ComponentType::BJT) std::static_pointer_cast<BJT>(comp)->resetState();
        }

        if (solveMna(currentV, sweepV)) {
            csvV << vVal; 
            for (int i = 0; i < numNodes - 1; ++i) {
                csvV << "," << sweepV[i];
            } 
            csvV << "\n";
            
            csvI << vVal; csvP << vVal;
            for (const auto& comp : components) {
                int ext = componentExtraIndices.count(comp) ? componentExtraIndices.at(comp) : -1;
                csvI << "," << comp->calculateCurrent(sweepV, ext); 
                csvP << "," << comp->calculatePower(sweepV, ext);
            }
            csvI << "\n"; csvP << "\n"; 
            currentV = sweepV;
        }
    }
}
