// freyalab/NoiseCircuit.h
#pragma once

#include "freyalab/Netlist.h"
#include "freyalab/CircuitComponents.h"
#include <string>
#include <utility>
#include <vector>

struct NoiseResultPoint {
    double frequencyHz = 0.0;
    double onoiseDensity_V_sqrtHz = 0.0; // Output noise voltage spectral density (V / sqrt(Hz))
    double inoiseDensity_V_sqrtHz = 0.0; // Equivalent input-referred noise density (V / sqrt(Hz))
    double transferGain_dB = 0.0;        // Small-signal forward transfer gain |Vout / Vin| (dB)
};

class NoiseCircuit {
private:
    Netlist& netlist;

public:
    explicit NoiseCircuit(Netlist& netlistInstance);

    /**
     * Executes small-signal .NOISE spectral density analysis across logarithmic frequency points
     * using the adjoint transposed MNA system A_c^T(omega) * z_adj = e_out.
     *
     * Computes:
     *   1. Output noise spectral density S_vout(f) in V / sqrt(Hz)
     *   2. Equivalent input-referred noise spectral density S_vin(f) in V / sqrt(Hz)
     *   3. Trapezoidally integrated total RMS output & input noise (V_rms)
     *   4. Exports spectrum dataset to "noise_spectrum_data.csv"
     */
    std::vector<NoiseResultPoint> runNoiseAnalysis(NetNode outPos,
                                                   NetNode outNeg,
                                                   const std::string& inputSourceName,
                                                   double startFreq,
                                                   double stopFreq,
                                                   int pointsCount);
};

