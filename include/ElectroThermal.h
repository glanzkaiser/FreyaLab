#pragma once
#include "Netlist.h"
#include "KLU.h"
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>


struct ThermalGridConfig {
    int nx = 48;
    int ny = 32;
    int nz = 2;                  // Layer 0: Silicon Active Die, Layer 1: Package / Heat Spreader
    double dieWidthM = 2.0e-3;   // 2.0 mm
    double dieHeightM = 1.5e-3;  // 1.5 mm
    double layerThickM = 100e-6; // 100 um per layer
    double kSilicon = 148.0;     // Thermal conductivity of Si [W / (m*K)]
    double kPackage = 380.0;     // Thermal conductivity of Cu spreader [W / (m*K)]
    double rhoCpSilicon = 1.63e6;// Volumetric heat capacity of Si [J / (m^3*K)]
    double rhoCpPackage = 3.45e6;// Volumetric heat capacity of Cu [J / (m^3*K)]
    double hConvection = 2500.0; // Bottom/heatsink convection coefficient [W / (m^2*K)]
    double ambientTempK = 300.15;// Reference ambient temperature [K] (27 C)
};

struct ThermalHeatSource {
    std::string componentName;
    double normX;      // Normalized [0.0, 1.0] X position on floorplan
    double normY;      // Normalized [0.0, 1.0] Y position on floorplan
    double powerWatts; // Dissipated Joule / switching power [W]
    int layerIdx = 0;  // Default: Layer 0 (Silicon Active Device Layer)
};

class ElectroThermalSolver {
public:
    explicit ElectroThermalSolver(const ThermalGridConfig& cfg = ThermalGridConfig{});

    // Steady-state 3D FVM heat conduction: G_th * dT = P
    std::vector<double> solveSteadyStateThermal(const std::vector<ThermalHeatSource>& sources);

    // Transient 3D FVM thermal step (Backward-Euler): (G_th + C_th/dt) * dT_new = P + (C_th/dt) * dT_prev
    std::vector<double> stepTransientThermal(const std::vector<ThermalHeatSource>& sources, double dt);

    // Resets temperature field back to ambientTempK
    void resetToAmbient();

    const ThermalGridConfig& getConfig() const { return config; }
    const std::vector<double>& getTemperatureFieldK() const { return tempFieldK; }
    double getMaxTempK() const { return maxTempK; }
    double getMinTempK() const { return minTempK; }

    // Helper: Extract 2D slice (nx * ny) for a given Z layer (0 = Si die, 1 = Package)
    std::vector<double> getLayerSliceK(int layerZ) const;

    // Inside ElectroThermalSolver in ElectroThermal.h:
    double sampleTemperatureAtK(double normX, double normY, int layerZ = 0) const {
    int ix = std::max(0, std::min(static_cast<int>(std::round(normX * (config.nx - 1))), config.nx - 1));
    int iy = std::max(0, std::min(static_cast<int>(std::round(normY * (config.ny - 1))), config.ny - 1));
    int iz = std::max(0, std::min(layerZ, config.nz - 1));
    int idx = (iz * config.ny + iy) * config.nx + ix;
    return (idx >= 0 && idx < static_cast<int>(tempFieldK.size())) ? tempFieldK[idx] : config.ambientTempK;
}

private:
    ThermalGridConfig config;
    std::vector<double> tempFieldK; // Current absolute temperature [K] per voxel (size nx*ny*nz)
    double minTempK = 300.15;
    double maxTempK = 300.15;

    // Shared 7-point FVM stencil builder & KLU/CSR-LU linear solver
    std::vector<double> assembleAndSolveFvm(const std::vector<ThermalHeatSource>& sources, double dtTransient);
};

