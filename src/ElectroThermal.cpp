#include "freyalab/ElectroThermal.h"
#include <algorithm>
#include <cmath>
#include <set>

ElectroThermalSolver::ElectroThermalSolver(const ThermalGridConfig& cfg)
    : config(cfg), minTempK(cfg.ambientTempK), maxTempK(cfg.ambientTempK)
{
    int totalCells = std::max(1, config.nx * config.ny * config.nz);
    tempFieldK.assign(totalCells, config.ambientTempK);
}

void ElectroThermalSolver::resetToAmbient() {
    std::fill(tempFieldK.begin(), tempFieldK.end(), config.ambientTempK);
    minTempK = config.ambientTempK;
    maxTempK = config.ambientTempK;
}

std::vector<double> ElectroThermalSolver::getLayerSliceK(int layerZ) const {
    int clampedZ = std::max(0, std::min(layerZ, config.nz - 1));
    int sliceSize = config.nx * config.ny;
    std::vector<double> slice(sliceSize, config.ambientTempK);
    int offset = clampedZ * sliceSize;
    for (int i = 0; i < sliceSize && (offset + i) < static_cast<int>(tempFieldK.size()); ++i) {
        slice[i] = tempFieldK[offset + i];
    }
    return slice;
}

std::vector<double> ElectroThermalSolver::solveSteadyStateThermal(const std::vector<ThermalHeatSource>& sources) {
    return assembleAndSolveFvm(sources, 0.0);
}

std::vector<double> ElectroThermalSolver::stepTransientThermal(const std::vector<ThermalHeatSource>& sources, double dt) {
    double safeDt = std::max(dt, 1e-15);
    return assembleAndSolveFvm(sources, safeDt);
}

std::vector<double> ElectroThermalSolver::assembleAndSolveFvm(const std::vector<ThermalHeatSource>& sources, double dtTransient) {
    const int nx = std::max(2, config.nx);
    const int ny = std::max(2, config.ny);
    const int nz = std::max(1, config.nz);
    const int totalCells = nx * ny * nz;

    if (static_cast<int>(tempFieldK.size()) != totalCells) {
        tempFieldK.assign(totalCells, config.ambientTempK);
    }

    // Voxel physical dimensions [m]
    const double dx = config.dieWidthM / static_cast<double>(nx);
    const double dy = config.dieHeightM / static_cast<double>(ny);
    const double dz = config.layerThickM;

    const double areaX = dy * dz; // Face area normal to X
    const double areaY = dx * dz; // Face area normal to Y
    const double areaZ = dx * dy; // Face area normal to Z
    const double voxelVol = dx * dy * dz;

    auto cellIndex = [nx, ny](int ix, int iy, int iz) -> int {
        return (iz * ny + iy) * nx + ix;
    };

    auto layerConductivity = [this](int iz) -> double {
        return (iz == 0) ? config.kSilicon : config.kPackage;
    };

    auto layerVolumetricCp = [this](int iz) -> double {
        return (iz == 0) ? config.rhoCpSilicon : config.rhoCpPackage;
    };

    // 1. Build 7-point FVM CSR sparsity pattern
    CsrMatrix G_th;
    G_th.numRows = totalCells;
    G_th.rowPointers.resize(totalCells + 1, 0);

    std::vector<std::set<int>> rowColSets(totalCells);
    for (int iz = 0; iz < nz; ++iz) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                int p = cellIndex(ix, iy, iz);
                rowColSets[p].insert(p);
                if (ix > 0)      rowColSets[p].insert(cellIndex(ix - 1, iy, iz));
                if (ix < nx - 1) rowColSets[p].insert(cellIndex(ix + 1, iy, iz));
                if (iy > 0)      rowColSets[p].insert(cellIndex(ix, iy - 1, iz));
                if (iy < ny - 1) rowColSets[p].insert(cellIndex(ix, iy + 1, iz));
                if (iz > 0)      rowColSets[p].insert(cellIndex(ix, iy, iz - 1));
                if (iz < nz - 1) rowColSets[p].insert(cellIndex(ix, iy, iz + 1));
            }
        }
    }

    int nnz = 0;
    for (int r = 0; r < totalCells; ++r) {
        G_th.rowPointers[r] = nnz;
        for (int c : rowColSets[r]) {
            G_th.colIndices.push_back(c);
            nnz++;
        }
    }
    G_th.rowPointers[totalCells] = nnz;
    G_th.values.assign(nnz, 0.0);

    // 2. Stamp 3D FVM Conductances (Harmonic mean across heterogeneous Si/Package interfaces)
    // State unknown is temperature rise dT = T - T_ambient [K]
    std::vector<double> rhsPower(totalCells, 0.0);

    for (int iz = 0; iz < nz; ++iz) {
        double kThis = layerConductivity(iz);
        double gx = (kThis * areaX) / dx;
        double gy = (kThis * areaY) / dy;

        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                int p = cellIndex(ix, iy, iz);
                double diagSum = 0.0;

                // -X and +X neighbors (adiabatic Neumann at die outer edges)
                if (ix > 0) {
                    int w = cellIndex(ix - 1, iy, iz);
                    G_th.addValue(p, w, -gx);
                    diagSum += gx;
                }
                if (ix < nx - 1) {
                    int e = cellIndex(ix + 1, iy, iz);
                    G_th.addValue(p, e, -gx);
                    diagSum += gx;
                }

                // -Y and +Y neighbors
                if (iy > 0) {
                    int s = cellIndex(ix, iy - 1, iz);
                    G_th.addValue(p, s, -gy);
                    diagSum += gy;
                }
                if (iy < ny - 1) {
                    int n = cellIndex(ix, iy + 1, iz);
                    G_th.addValue(p, n, -gy);
                    diagSum += gy;
                }

                // -Z neighbor (toward Silicon active layer)
                if (iz > 0) {
                    double kDown = layerConductivity(iz - 1);
                    double kHarm = (2.0 * kThis * kDown) / std::max(kThis + kDown, 1e-12);
                    double gzDown = (kHarm * areaZ) / dz;
                    int b = cellIndex(ix, iy, iz - 1);
                    G_th.addValue(p, b, -gzDown);
                    diagSum += gzDown;
                }

                // +Z neighbor (toward Package / Heat Spreader) or Convective Heatsink Boundary
                if (iz < nz - 1) {
                    double kUp = layerConductivity(iz + 1);
                    double kHarm = (2.0 * kThis * kUp) / std::max(kThis + kUp, 1e-12);
                    double gzUp = (kHarm * areaZ) / dz;
                    int t = cellIndex(ix, iy, iz + 1);
                    G_th.addValue(p, t, -gzUp);
                    diagSum += gzUp;
                } else {
                    // Top package surface Robin convective boundary: q = h_conv * A * (T - T_amb)
                    double gConv = config.hConvection * areaZ;
                    diagSum += gConv;
                }

                // Subtle lateral edge convection so even single-layer grids remain strictly SPD
                if (ix == 0 || ix == nx - 1) diagSum += (0.05 * config.hConvection * areaX);
                if (iy == 0 || iy == ny - 1) diagSum += (0.05 * config.hConvection * areaY);

                // Optional transient thermal capacitance term: C_th / dt
                if (dtTransient > 0.0) {
                    double cTh = layerVolumetricCp(iz) * voxelVol;
                    double gCap = cTh / dtTransient;
                    diagSum += gCap;
                    double prevRise = tempFieldK[p] - config.ambientTempK;
                    rhsPower[p] += gCap * prevRise;
                }

                G_th.addValue(p, p, diagSum);
            }
        }
    }

    // 3. Inject Localized Component Joule Power (bilinearly distributed onto a 2x2 footprint)
    for (const auto& src : sources) {
        if (src.powerWatts <= 0.0) continue;

        double u = std::max(0.0, std::min(src.normX, 1.0)) * (nx - 1);
        double v = std::max(0.0, std::min(src.normY, 1.0)) * (ny - 1);

        int ix0 = static_cast<int>(std::floor(u));
        int iy0 = static_cast<int>(std::floor(v));
        int ix1 = std::min(ix0 + 1, nx - 1);
        int iy1 = std::min(iy0 + 1, ny - 1);

        double fx = u - ix0;
        double fy = v - iy0;
        int iz = std::max(0, std::min(src.layerIdx, nz - 1));

        rhsPower[cellIndex(ix0, iy0, iz)] += src.powerWatts * (1.0 - fx) * (1.0 - fy);
        rhsPower[cellIndex(ix1, iy0, iz)] += src.powerWatts * fx * (1.0 - fy);
        rhsPower[cellIndex(ix0, iy1, iz)] += src.powerWatts * (1.0 - fx) * fy;
        rhsPower[cellIndex(ix1, iy1, iz)] += src.powerWatts * fx * fy;
    }

    // 4. Solve Sparse FVM System G_th * dT = P
    // Match FreyaLab's hybrid architecture: Left-Looking CSR LU for small grids (< 300 cells),
    // and BTF + AMD KLU for larger 3D FVM grids (>= 300 cells).
    std::vector<double> deltaT(totalCells, 0.0);
    bool solved = false;

    if (totalCells >= KluSolver::KLU_ACTIVATION_THRESHOLD) {
        KluSymbolicContext kluCache = KluSolver::analyzeSymbolic(G_th);
        if (KluSolver::factorizeNumeric(G_th, kluCache)) {
            deltaT = KluSolver::solveNumeric(kluCache, rhsPower);
            solved = true;
        }
    }

    if (!solved) {
        CsrMatrix L, U;
        std::vector<int> P(totalCells, 0);
        L.initializeFromPattern(G_th);
        U.initializeFromPattern(G_th);
        if (CircuitMath::luFactorize(G_th, L, U, P)) {
            deltaT = CircuitMath::luSolve(L, U, P, rhsPower);
            solved = true;
        }
    }

    // 5. Convert Temperature Rise dT back to Absolute Kelvin T(x, y, z) and track extrema
    minTempK = 1e30;
    maxTempK = -1e30;
    for (int i = 0; i < totalCells; ++i) {
        double rise = (solved && std::isfinite(deltaT[i])) ? std::max(0.0, deltaT[i]) : 0.0;
        double absK = config.ambientTempK + rise;
        tempFieldK[i] = absK;
        if (absK < minTempK) minTempK = absK;
        if (absK > maxTempK) maxTempK = absK;
    }

    if (minTempK > maxTempK) {
        minTempK = config.ambientTempK;
        maxTempK = config.ambientTempK;
    }

    return tempFieldK;
}
