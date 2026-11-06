#include "freyalab/ACCircuit.h"
#include "freyalab/math.h"
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

ACCircuit::ACCircuit(const std::vector<std::shared_ptr<Component>>& parsedComponents, int nodeCount)
    : components(parsedComponents), numNodes(nodeCount), totalSize(nodeCount) {
    int matrixNodeRows = numNodes - 1; 
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

bool ACCircuit::isAcCircuit() const {
    for (const auto& comp : components) {
        if (comp->getType() == ComponentType::VOLTAGE_SOURCE) {
            if (std::static_pointer_cast<VoltageSource>(comp)->getSourceType() == SourceType::AC) return true;
        }
    }
    return false;
}

bool ACCircuit::solveTransientStep(double time, double dt, const std::vector<double>& oldV, std::vector<double>& outV) {
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
            int assignedExtraIdx = componentExtraIndices.count(comp) ? componentExtraIndices[comp] : -1;
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
    
    // FIX 3: Commit state boundary vectors exclusively when convergence validation passes
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

void ACCircuit::runTransientWaveformAnalysis(double tEnd, double dt) {
    std::cout << "\n====================================================================================\n"
              << "                            AC TRANSIENT WAVEFORM TABLE                             \n"
              << "====================================================================================\n"
              << std::left << std::setw(12) << "Time(s)" << std::setw(15) << "V(input)" << std::setw(15) << "V(midpoint)" << std::setw(15) << "V(output)" << std::setw(18) << "I(L1) current" << "\n"
              << "------------------------------------------------------------------------------------\n";
              
    std::vector<double> oldV(totalSize, 0.0), currentV(totalSize, 0.0); 
    double time = 0.0;
    
    std::shared_ptr<Inductor> l1 = nullptr;
    for (auto& comp : components) {
        if (comp->getName() == "L1") l1 = std::static_pointer_cast<Inductor>(comp);
    }
    
    while (time <= tEnd + 1e-9) {
        if (solveTransientStep(time, dt, oldV, currentV)) {
            int lExtra = (l1 && componentExtraIndices.count(l1)) ? componentExtraIndices[l1] : -1;
            
            // FIX 1: Explicit safe index protection blocks for Ground configurations
            double vInReal   = (numNodes >= 2) ? currentV[0] : 0.0;
            double vMidReal  = (numNodes >= 3) ? currentV[1] : 0.0;
            double vOutReal  = (numNodes >= 4) ? currentV[2] : 0.0;
            double iL1Real   = l1 ? l1->calculateCurrent(currentV, lExtra) : 0.0;

            std::cout << std::scientific << std::setprecision(4) << std::left << std::setw(12) << time << std::fixed << std::setprecision(5)
                      << std::setw(15) << vInReal 
                      << std::setw(15) << vMidReal
                      << std::setw(15) << vOutReal 
                      << std::scientific << std::setw(18) << iL1Real << "\n";
                      
            oldV = currentV; 
            time += dt;
        } else {
            std::cerr << "[Fatal Error] Transient convergence failed.\n";
            break;
        }
    }
    std::cout << "====================================================================================\n";
}

void ACCircuit::runFrequencySweepAnalysis(double startFreq, double endFreq, int pointsPerDecade) {
    std::cout << "\n==============================================================================================\n"
              << "                                  AC SIMULATION OUTPUT TABLE                                  \n"
              << "==============================================================================================\n"
              << std::left << std::setw(12) << "Freq(Hz)" << std::setw(14) << "V(2) Mag(V)" << std::setw(14) << "V(2) Phase(°)" << std::setw(14) << "V(3) Mag(V)" << std::setw(14) << "V(3) Phase(°)" << std::setw(15) << "I(L1) Mag(A)" << std::setw(15) << "I(L1) Phase(°)" << "\n"
              << "----------------------------------------------------------------------------------------------\n";
              
    std::ofstream bodeCsv("bode_plot_data.csv");
    bodeCsv << "Frequency(Hz),V2_Magnitude(dB),V2_Phase(deg),V3_Magnitude(dB),V3_Phase(deg),IL1_Magnitude(A),IL1_Phase(deg)\n";
    
    std::shared_ptr<Inductor> l1 = nullptr; 
    for (auto& comp : components) { 
        if (comp->getName() == "L1") l1 = std::static_pointer_cast<Inductor>(comp); 
    }
    int lExtra = (l1 && componentExtraIndices.count(l1)) ? componentExtraIndices[l1] : -1;
    
    double currentFreq = startFreq;
    while (currentFreq <= endFreq * 1.0001) {
        double period = 1.0 / currentFreq;
        
        // FIX: Maintain uniform fine-grain stepping to avoid time-step variance errors
        double settleDt = period / 100.0;
        std::vector<double> stateV(totalSize, 0.0);
        double time = 0.0;
        bool ok = true;
        
        // Settle for 5 complete periods to fully flush out start-up transients
        for (int step = 0; step < 500; ++step) {
            std::vector<double> nextV(totalSize, 0.0);
            if (!solveTransientStep(time + settleDt, settleDt, stateV, nextV)) {
                ok = false;
                break;
            }
            stateV = nextV;
            time += settleDt;
        }

        if (!ok) {
            std::cerr << "[Warning] Solver diverged at frequency: " << currentFreq << " Hz (Skipping)\n";
            currentFreq *= std::pow(10.0, 1.0 / pointsPerDecade);
            continue;
        }
        
        // Measure real time peak values directly by tracking extreme amplitudes over a full 6th cycle
        double v2_max = 0.0, v3_max = 0.0, il1_max = 0.0;
        double t_v2_peak = 0.0, t_v3_peak = 0.0, t_il1_peak = 0.0;

        for (int step = 0; step < 100; ++step) {
            std::vector<double> nextV(totalSize, 0.0);
            if (solveTransientStep(time + settleDt, settleDt, stateV, nextV)) {
                double v2_val = (totalSize >= 2) ? nextV[1] : 0.0; // Node 2 voltage
                double v3_val = (totalSize >= 3) ? nextV[2] : 0.0; // Node 3 voltage
                double il1_val = l1 ? l1->calculateCurrent(nextV, lExtra) : 0.0;

                if (std::abs(v2_val) > v2_max) { v2_max = std::abs(v2_val); t_v2_peak = time; }
                if (std::abs(v3_val) > v3_max) { v3_max = std::abs(v3_val); t_v3_peak = time; }
                if (std::abs(il1_val) > il1_max) { il1_max = std::abs(il1_val); t_il1_peak = time; }
                
                stateV = nextV;
                time += settleDt;
            }
        }

        // Calculate phase angles relative to the input AC source period mapping boundaries
        double v2_phase = std::fmod(t_v2_peak, period) / period * 360.0;
        double v3_phase = std::fmod(t_v3_peak, period) / period * 360.0;
        double il1_phase = std::fmod(t_il1_peak, period) / period * 360.0;
        
        if (v2_phase > 180.0) v2_phase -= 360.0;
        if (v3_phase > 180.0) v3_phase -= 360.0;
        if (il1_phase > 180.0) il1_phase -= 360.0;

        double v2_db = 20.0 * std::log10(v2_max + 1e-15);
        double v3_db = 20.0 * std::log10(v3_max + 1e-15);

        std::cout << std::scientific << std::setprecision(3) << std::left << std::setw(12) << currentFreq << std::fixed << std::setprecision(4)
                  << std::setw(14) << v2_max << std::setw(14) << v2_phase
                  << std::setw(14) << v3_max << std::setw(14) << v3_phase
                  << std::scientific << std::setw(15) << il1_max << std::fixed << std::setw(15) << il1_phase << "\n";
                  
        bodeCsv << currentFreq << "," << v2_db << "," << v2_phase << "," << v3_db << "," << v3_phase << "," << il1_max << "," << il1_phase << "\n";
        currentFreq *= std::pow(10.0, 1.0 / pointsPerDecade);
    }
    std::cout << "==============================================================================================\n";
    bodeCsv.close();
}
