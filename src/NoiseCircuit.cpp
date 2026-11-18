// NoiseCircuit.cpp
#include "freyalab/NoiseCircuit.h"
#include "freyalab/math.h"
#include "freyalab/KLU.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double K_BOLTZMANN = 1.380649e-23;     // J/K
constexpr double Q_ELECTRON  = 1.602176634e-19;  // Coulombs

// Builds the real CSR sparsity pattern transpose A^T from real CsrMatrix A
CsrMatrix transposeRealPattern(const CsrMatrix& A) {
    int n = A.numRows;
    CsrMatrix AT;
    AT.numRows = n;
    AT.rowPointers.assign(n + 1, 0);

    if (n <= 0) return AT;

    for (int r = 0; r < n; ++r) {
        for (int idx = A.rowPointers[r]; idx < A.rowPointers[r + 1]; ++idx) {
            int c = A.colIndices[idx];
            if (c >= 0 && c < n) {
                AT.rowPointers[c + 1]++;
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        AT.rowPointers[i + 1] += AT.rowPointers[i];
    }

    int nnz = AT.rowPointers[n];
    AT.colIndices.assign(nnz, 0);
    AT.values.assign(nnz, 0.0);

    std::vector<int> writeOffset(AT.rowPointers.begin(), AT.rowPointers.begin() + n);
    for (int r = 0; r < n; ++r) {
        for (int idx = A.rowPointers[r]; idx < A.rowPointers[r + 1]; ++idx) {
            int c = A.colIndices[idx];
            if (c >= 0 && c < n) {
                int dest = writeOffset[c]++;
                AT.colIndices[dest] = r;
                AT.values[dest] = A.values[idx];
            }
        }
    }

    return AT;
}

// Builds the exact transpose A^T of a square complex CSR matrix for adjoint noise sensitivity solves
CsrMatrixComplex transposeComplexCsr(const CsrMatrixComplex& A) {
    int n = A.numRows;
    CsrMatrixComplex AT;
    AT.numRows = n;
    AT.rowPointers.assign(n + 1, 0);

    if (n <= 0) return AT;

    // 1. Count nonzeros per transposed row (original column)
    for (int r = 0; r < n; ++r) {
        for (int idx = A.rowPointers[r]; idx < A.rowPointers[r + 1]; ++idx) {
            int c = A.colIndices[idx];
            if (c >= 0 && c < n) {
                AT.rowPointers[c + 1]++;
            }
        }
    }

    // 2. Prefix sum to form rowPointers
    for (int i = 0; i < n; ++i) {
        AT.rowPointers[i + 1] += AT.rowPointers[i];
    }

    int nnz = AT.rowPointers[n];
    AT.colIndices.assign(nnz, 0);
    AT.values.assign(nnz, std::complex<double>(0.0, 0.0));

    // 3. Scatter transposed entries
    std::vector<int> writeOffset(AT.rowPointers.begin(), AT.rowPointers.begin() + n);
    for (int r = 0; r < n; ++r) {
        for (int idx = A.rowPointers[r]; idx < A.rowPointers[r + 1]; ++idx) {
            int c = A.colIndices[idx];
            if (c >= 0 && c < n) {
                int dest = writeOffset[c]++;
                AT.colIndices[dest] = r;
                AT.values[dest] = A.values[idx];
            }
        }
    }

    return AT;
}

} // namespace

NoiseCircuit::NoiseCircuit(Netlist& netlistInstance) : netlist(netlistInstance) {}

std::vector<NoiseResultPoint> NoiseCircuit::runNoiseAnalysis(NetNode outPos,
                                                             NetNode outNeg,
                                                             const std::string& inputSourceName,
                                                             double startFreq,
                                                             double stopFreq,
                                                             int pointsCount) {
    int systemSize = netlist.getTotalSystemSize();
    if (systemSize == 0) {
        netlist.buildSparsityPattern();
        systemSize = netlist.getTotalSystemSize();
    }

    std::vector<NoiseResultPoint> results;
    if (systemSize <= 0 || pointsCount <= 0) return results;

    std::string outPosLabel = outPos.isGround() ? "0" : netlist.getNodeName(outPos);
    std::string outLabel = outNeg.isGround()
                               ? ("V(" + outPosLabel + ")")
                               : ("V(" + outPosLabel + "," + netlist.getNodeName(outNeg) + ")");

    std::cout << "\n========================================================================================\n"
              << "                  SMALL-SIGNAL ADJOINT .NOISE SPECTRAL DENSITY ANALYSIS                 \n"
              << "========================================================================================\n"
              << "Output Probe : " << outLabel
              << " | Reference Input : " << inputSourceName
              << " | Temp : " << std::fixed << std::setprecision(2) << netlist.getTemperatureCelsius() << " C\n"
              << "----------------------------------------------------------------------------------------\n"
              << std::left << std::setw(15) << "Freq(Hz)"
              << std::setw(22) << "Onoise(V/rtHz)"
              << std::setw(22) << "Inoise(V/rtHz)"
              << std::setw(18) << "Gain(dB)" << "\n"
              << "----------------------------------------------------------------------------------------\n";

    std::ofstream noiseCsv("noise_spectrum_data.csv");
    if (!noiseCsv.is_open()) {
        throw std::runtime_error("Fatal Error: Unable to open noise_spectrum_data.csv for output.");
    }
    noiseCsv << "Frequency(Hz),Onoise(V/sqrtHz),Inoise(V/sqrtHz),Gain(dB)\n";

    // 1. Solve DC Operating Point at global .TEMP prior to small-signal linearization
    MnaMatrixContext dcCtx{
        netlist.getGMatrix(),
        netlist.getBVector(),
        &netlist.getStateV(),
        &netlist.getHistoricV(),
        &netlist.getMidpointV(),
        0.0,
        0.0,
        0,
        MnaIndex{},
        netlist.getTemperatureCelsius(),
        netlist.getThermalVoltage()
    };

    if (!netlist.runNewtonRaphson(dcCtx, 150, 1e-3, 1e-6)) {
        std::cerr << "Warning: DC Operating Point did not fully converge prior to .NOISE analysis.\n";
    }

    const std::vector<double> dcState = netlist.getCurrentState();
    SimulationContext dcSimCtx{ dcState, 0.0, 0 };

    const double tempKelvin = netlist.getTemperatureKelvin();
    const double fourKT     = 4.0 * K_BOLTZMANN * tempKelvin;
    const double vtThermal  = netlist.getThermalVoltage();

    // Locate the specified input source extra MNA row / AC excitation amplitude
    std::string targetInputUpper = inputSourceName;
    std::transform(targetInputUpper.begin(), targetInputUpper.end(), targetInputUpper.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

    MnaIndex inputSourceExtraIdx{};
    NetNode inputPosNode{0};
    NetNode inputNegNode{0};
    ComponentType inputSourceType = ComponentType::VOLTAGE_SOURCE;

    netlist.inspectComponents([&](const Component& comp, MnaIndex extraIdx) {
        std::string cName = comp.getName();
        std::transform(cName.begin(), cName.end(), cName.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
        if (cName == targetInputUpper) {
            inputSourceExtraIdx = extraIdx;
            inputSourceType = comp.getType();
            const auto& nList = comp.getNodes();
            if (!nList.empty()) inputPosNode = nList[0];
            if (nList.size() >= 2) inputNegNode = nList[1];
        }
    });

    // 2. Initialize complex MNA matrix buffers
    bool useKlu = netlist.isUsingKluBackend();
    CsrMatrixComplex A_complex, L_complex, U_complex;
    std::vector<int> P_complex;
    A_complex.initializeFromPattern(netlist.getGMatrix());
    if (!useKlu) {
        L_complex.initializeFromPattern(netlist.getGMatrix());
        U_complex.initializeFromPattern(netlist.getGMatrix());
        P_complex.resize(systemSize, 0);
    }

    CsrMatrix patternAT = transposeRealPattern(netlist.getGMatrix());
    KluSymbolicContext adjSym;
    CsrMatrixComplex LT_complex, UT_complex;
    std::vector<int> PT_complex;

    if (useKlu) {
        adjSym = KluSolver::analyzeSymbolic(patternAT);
    } else {
        LT_complex.initializeFromPattern(patternAT);
        UT_complex.initializeFromPattern(patternAT);
        PT_complex.resize(systemSize, 0);
    }

    std::vector<std::complex<double>> B_fwd(systemSize, std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> E_adj(systemSize, std::complex<double>(0.0, 0.0));

    int numNodes = netlist.getNumNodes();
    auto getSafePhasor = [numNodes](const std::vector<std::complex<double>>& X, NetNode node) -> std::complex<double> {
        if (node.isGround()) return std::complex<double>(0.0, 0.0);
        MnaIndex idx = MnaIndex::fromNetNode(node);
        if (idx.isValid() && idx.idx < numNodes && static_cast<size_t>(idx.idx) < X.size()) {
            return X[idx.idx];
        }
        return std::complex<double>(0.0, 0.0);
    };

    // Build unit adjoint excitation vector e_out (+1 at outPos, -1 at outNeg)
    MnaIndex outPosIdx = MnaIndex::fromNetNode(outPos);
    MnaIndex outNegIdx = MnaIndex::fromNetNode(outNeg);
    if (outPosIdx.isValid() && outPosIdx.idx < numNodes) {
        E_adj[outPosIdx.idx] = std::complex<double>(1.0, 0.0);
    }
    if (outNegIdx.isValid() && outNegIdx.idx < numNodes) {
        E_adj[outNegIdx.idx] = std::complex<double>(-1.0, 0.0);
    }

    double logStart = std::log10(std::max(startFreq, 1e-3));
    double logStop  = std::log10(std::max(stopFreq, 1e-3));
    double logStep  = (pointsCount > 1) ? (logStop - logStart) / (pointsCount - 1) : 0.0;

    double integratedOnoiseSq = 0.0;
    double integratedInoiseSq = 0.0;
    double prevFreq = 0.0;
    double prevOnoiseDensitySq = 0.0;
    double prevInoiseDensitySq = 0.0;

    // 3. Logarithmic Frequency Sweep
    for (int i = 0; i < pointsCount; ++i) {
        double freq  = std::pow(10.0, logStart + i * logStep);
        double omega = 2.0 * M_PI * freq;

        A_complex.zero();
        std::fill(B_fwd.begin(), B_fwd.end(), std::complex<double>(0.0, 0.0));

        // Stamp small-signal AC conductances/admittances
        netlist.inspectComponents([&A_complex, &B_fwd, omega, &dcState](Component& comp, MnaIndex extraIdx) {
            (void)extraIdx;
            comp.stampAC(A_complex, B_fwd, omega, dcState);
        });

        // Override forward excitation with a clean 1.0 V (or 1.0 A) stimulus on inputSourceName
        std::fill(B_fwd.begin(), B_fwd.end(), std::complex<double>(0.0, 0.0));
        if (inputSourceExtraIdx.isValid() && static_cast<size_t>(inputSourceExtraIdx.idx) < B_fwd.size()) {
            B_fwd[inputSourceExtraIdx.idx] = std::complex<double>(1.0, 0.0);
        } else if (inputSourceType == ComponentType::CURRENT_SOURCE) {
            MnaIndex mInP = MnaIndex::fromNetNode(inputPosNode);
            MnaIndex mInN = MnaIndex::fromNetNode(inputNegNode);
            if (mInP.isValid() && mInP.idx < numNodes) B_fwd[mInP.idx] -= std::complex<double>(1.0, 0.0);
            if (mInN.isValid() && mInN.idx < numNodes) B_fwd[mInN.idx] += std::complex<double>(1.0, 0.0);
        }

        // Solve forward transfer H(jw) = Vout(jw) / Vin(jw)
        std::vector<std::complex<double>> X_fwd;
        if (useKlu) {
            if (!KluSolver::factorizeNumericComplex(A_complex, netlist.getKluCache())) {
                std::cerr << "Warning: Forward KLU factorization failed during .NOISE at " << freq << " Hz.\n";
                continue;
            }
            X_fwd = KluSolver::solveNumericComplex(netlist.getKluCache(), B_fwd);
        } else {
            P_complex.assign(systemSize, 0);
            if (!CircuitMath::luFactorizeComplex(A_complex, L_complex, U_complex, P_complex)) {
                std::cerr << "Warning: Forward AC factorization failed during .NOISE at " << freq << " Hz.\n";
                continue;
            }
            X_fwd = CircuitMath::luSolveComplex(L_complex, U_complex, P_complex, B_fwd);
        }

        std::complex<double> vOutPhasor = getSafePhasor(X_fwd, outPos) - getSafePhasor(X_fwd, outNeg);
        double forwardGainMag = std::max(std::abs(vOutPhasor), 1e-18);
        double forwardGainDb  = 20.0 * std::log10(forwardGainMag);

        // Solve adjoint transposed system A_c^T(jw) * Z_adj = E_adj
        CsrMatrixComplex AT_complex = transposeComplexCsr(A_complex);
        std::vector<std::complex<double>> Z_adj;

        if (useKlu) {
            if (!KluSolver::factorizeNumericComplex(AT_complex, adjSym)) {
                std::cerr << "Warning: Adjoint KLU factorization failed during .NOISE at " << freq << " Hz.\n";
                continue;
            }
            Z_adj = KluSolver::solveNumericComplex(adjSym, E_adj);
        } else {
            PT_complex.assign(systemSize, 0);
            if (!CircuitMath::luFactorizeComplex(AT_complex, LT_complex, UT_complex, PT_complex)) {
                std::cerr << "Warning: Adjoint transposed factorization failed during .NOISE at " << freq << " Hz.\n";
                continue;
            }
            Z_adj = CircuitMath::luSolveComplex(LT_complex, UT_complex, PT_complex, E_adj);
        }


        // 4. Accumulate uncorrelated device noise PSDs weighted by |Z_adj(n+) - Z_adj(n-)|^2
        double totalOutputNoisePsd = 0.0; // Units: V^2 / Hz

        netlist.inspectComponents([&](const Component& comp, MnaIndex extraIdx) {
            ComponentType cType = comp.getType();
            const auto& cNodes = comp.getNodes();
            if (cNodes.size() < 2) return;

            if (cType == ComponentType::RESISTOR) {
                double R = std::max(comp.getValue(), 1e-12);
                double siThermal = fourKT / R; // Johnson-Nyquist thermal current noise (A^2/Hz)

                std::complex<double> zDiff = getSafePhasor(Z_adj, cNodes[0]) - getSafePhasor(Z_adj, cNodes[1]);
                totalOutputNoisePsd += siThermal * std::norm(zDiff);
            }
            else if (cType == ComponentType::DIODE) {
                double id = std::abs(comp.calculateCurrent(dcSimCtx, extraIdx));
                double siShot = 2.0 * Q_ELECTRON * id; // Shot noise (A^2/Hz)

                std::complex<double> zDiff = getSafePhasor(Z_adj, cNodes[0]) - getSafePhasor(Z_adj, cNodes[1]);
                totalOutputNoisePsd += siShot * std::norm(zDiff);
            }
            else if (cType == ComponentType::BJT) {
                if (cNodes.size() < 3) return;
                double ic = std::abs(comp.calculateCurrent(dcSimCtx, extraIdx));
                double ib = ic / 100.0;

                // Collector-Emitter shot noise + Base-Emitter shot noise
                double siCollector = 2.0 * Q_ELECTRON * ic;
                double siBase      = 2.0 * Q_ELECTRON * ib;

                std::complex<double> zCE = getSafePhasor(Z_adj, cNodes[0]) - getSafePhasor(Z_adj, cNodes[2]);
                std::complex<double> zBE = getSafePhasor(Z_adj, cNodes[1]) - getSafePhasor(Z_adj, cNodes[2]);

                totalOutputNoisePsd += siCollector * std::norm(zCE) + siBase * std::norm(zBE);
            }
            else if (cType == ComponentType::MOSFET ||
                     cType == ComponentType::FINFET ||
                     cType == ComponentType::RIBBONFET ||
                     cType == ComponentType::GAA_FET) {
                if (cNodes.size() < 3) return;
                double id = std::abs(comp.calculateCurrent(dcSimCtx, extraIdx));
                double vov = std::max(std::abs(getSafeVoltage(dcState, cNodes[1]) - getSafeVoltage(dcState, cNodes[2])) - 0.3, 0.05);
                double gmApprox = (id > 1e-15) ? std::max((2.0 * id) / vov, id / (1.5 * vtThermal)) : 0.0;

                // Channel thermal noise (8/3 * kT * gm) + 1/f flicker noise approximation
                constexpr double KF_DEFAULT = 1e-24;
                double siChannelThermal = (8.0 / 3.0) * K_BOLTZMANN * tempKelvin * gmApprox;
                double siFlicker        = (KF_DEFAULT * id) / std::max(freq, 1.0);

                std::complex<double> zDS = getSafePhasor(Z_adj, cNodes[0]) - getSafePhasor(Z_adj, cNodes[2]);
                totalOutputNoisePsd += (siChannelThermal + siFlicker) * std::norm(zDS);
            }
        });

        double onoiseDensity = std::sqrt(std::max(totalOutputNoisePsd, 0.0));
        double inoiseDensity = onoiseDensity / forwardGainMag;

                // Trapezoidal integration across frequency bandwidth for total RMS noise
        if (i > 0) {
            double df = freq - prevFreq;
            integratedOnoiseSq += 0.5 * (prevOnoiseDensitySq + totalOutputNoisePsd) * df;
            double inoisePsd = inoiseDensity * inoiseDensity;
            integratedInoiseSq += 0.5 * (prevInoiseDensitySq + inoisePsd) * df;
        }
        prevFreq = freq;
        prevOnoiseDensitySq = totalOutputNoisePsd;
        prevInoiseDensitySq = inoiseDensity * inoiseDensity;

        NoiseResultPoint pt;
        pt.frequencyHz = freq;
        pt.onoiseDensity_V_sqrtHz = onoiseDensity;
        pt.inoiseDensity_V_sqrtHz = inoiseDensity;
        pt.transferGain_dB = forwardGainDb;
        results.push_back(pt);

        std::cout << std::left << std::scientific << std::setprecision(4)
                  << std::setw(15) << freq
                  << std::setw(22) << onoiseDensity
                  << std::setw(22) << inoiseDensity
                  << std::fixed << std::setprecision(2)
                  << std::setw(18) << forwardGainDb << "\n";

        noiseCsv << std::scientific << std::setprecision(6)
                 << freq << "," << onoiseDensity << "," << inoiseDensity << "," << forwardGainDb << "\n";
    }

    noiseCsv.close();

    double totalOnoiseRms = std::sqrt(std::max(integratedOnoiseSq, 0.0));
    double totalInoiseRms = std::sqrt(std::max(integratedInoiseSq, 0.0));

    std::cout << "----------------------------------------------------------------------------------------\n"
              << std::scientific << std::setprecision(4)
              << "Total Integrated Output Noise (V_rms) : " << totalOnoiseRms << " V\n"
              << "Total Integrated Input Noise  (V_rms) : " << totalInoiseRms << " V\n"
              << "========================================================================================\n"
              << "Adjoint .NOISE spectral analysis completed. Dataset exported to noise_spectrum_data.csv\n\n";

    return results;
}


