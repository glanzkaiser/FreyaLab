// BSIM.h
#pragma once

#include "freyalab/CircuitComponents.h"
#include "freyalab/osdi.h"
#include "freyalab/math.h"
#include <cmath>
#include <complex>
#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

enum class MultiGateGeometry {
    FINFET_TRIGATE = 1, // GEOMOD = 1: 3-sided Tri-Gate FinFET
    GAA_NANOSHEET  = 2  // GEOMOD = 2: 4-sided Gate-All-Around RibbonFET / Nanosheet
};

struct BsimCmgParameters {
    MultiGateGeometry geomod = MultiGateGeometry::FINFET_TRIGATE;
    bool isNChannel = true;

    // 3D Multi-Gate Physical Dimensions
    double L       = 14e-9;   // Effective channel length (14 nm)
    double TFIN    = 7e-9;    // Fin width or nanosheet thickness (7 nm)
    double HFIN    = 30e-9;   // Fin height for FinFET (30 nm)
    double WSHEET  = 25e-9;   // Nanosheet width for RibbonFET / GAA-FET (25 nm)
    int    NFIN    = 4;       // Number of parallel fins or stacked nanosheets
    double EOT     = 0.8e-9;  // Equivalent Oxide Thickness (0.8 nm HKMG)

    // Threshold, Subthreshold & Short-Channel DIBL Parameters
    double VTH0    = 0.35;    // Long-channel threshold voltage (V)
    double NFACTOR = 1.2;     // Subthreshold swing ideality factor
    double ETA0    = 0.08;    // Drain-Induced Barrier Lowering (DIBL) coefficient
    double DVT0    = 0.04;    // Short-channel Vth roll-off offset (V)

    // Mobility, Velocity Saturation & Output Conductance Parameters
    double U0      = 0.035;   // Low-field carrier mobility (m^2 / V*s)
    double VSAT    = 1.0e5;   // Carrier saturation velocity (m/s)
    double PCLM    = 0.05;    // Channel-length modulation parameter
    double DELTA   = 0.01;    // Smoothing parameter for Vdseff transition (V)

    // Per-Fin / Per-Sheet Parasitic Overlap & Fringe Capacitances
    double CGS0    = 0.4e-15; // Gate-to-Source capacitance per fin/sheet (F)
    double CGD0    = 0.4e-15; // Gate-to-Drain capacitance per fin/sheet (F)

    // Raw key=value map forwarded directly to BSIM-CMG 112.1.0 OSDI parameters
    std::map<std::string, double> rawOverrides;
};

/**
 * Singleton-capable dynamic library loader for OpenVAF-compiled bsimcmg.osdi.
 */
class BsimCmgOsdiLoader {
private:
    void* libHandle = nullptr;
    const OsdiDescriptor* descriptor = nullptr;
    bool attemptedLoad = false;

public:
    BsimCmgOsdiLoader() = default;
    ~BsimCmgOsdiLoader();

    BsimCmgOsdiLoader(const BsimCmgOsdiLoader&) = delete;
    BsimCmgOsdiLoader& operator=(const BsimCmgOsdiLoader&) = delete;

    // Shared global instance (searches "models/bsimcmg.osdi" and "./bsimcmg.osdi" automatically)
    static BsimCmgOsdiLoader& getGlobalInstance();

    bool load(const std::string& osdiPath);
    void unload();

    bool isLoaded() const { return descriptor != nullptr; }
    const OsdiDescriptor* getDescriptor() const { return descriptor; }
};


/**
 * Multi-Gate Transistor Component (FinFET, RibbonFET, GAA-FET).
 * Automatically uses UC Berkeley BSIM-CMG 112.1.0 via OSDI if bsimcmg.osdi is loaded,
 * and seamlessly falls back to the built-in Phase-1 analytical model otherwise.
 */
class BsimCmgFet : public Component {
private:
    BsimCmgParameters params;
    double weff = 0.0;
    double cox  = 0.0;

    // --- OSDI Runtime State Buffers ---
    const OsdiDescriptor* osdiDesc = nullptr;
    std::vector<uint8_t> modelData;
    std::vector<uint8_t> instanceData;
    std::vector<double> stateHistN;       // Converged OSDI internal states at step n
    std::vector<double> stateHistNc;      // Converged OSDI internal states at TR midpoint n+gamma
    std::vector<double> stateCheckpoint;
    bool osdiReady = false;

    void setupOsdiBuffers();
    void setOsdiParamByName(const std::string& paramName, double val);
    MnaIndex resolveMnaRowIndex(uint32_t osdiNodeIdx, MnaIndex extraIdx) const;

    // --- Built-In Phase-1 Analytical Fallback Kernel ---
    void evaluatePhase1Physics(double vgs, double vds,
                               double& ids, double& gm, double& gds) const;
    void stampPhase1(MnaMatrixContext& ctx);
    void stampACPhase1(CsrMatrixComplex& A_c, double omega, const std::vector<double>& dcState);

public:
    BsimCmgFet(const std::string& name,
               const std::vector<NetNode>& nodes,
               ComponentType subType = ComponentType::FINFET,
               const BsimCmgParameters& modelParams = BsimCmgParameters());

    int getExtraRowRequirement() const override;
    void allocateRows(int& globalExtraRowCounter) override;

    void stamp(MnaMatrixContext& ctx) override;
    void updateHistory(const MnaMatrixContext& ctx) override;
    void getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns,
                            MnaIndex extraIdx = MnaIndex{}) const override;
    void stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                 double omega, const std::vector<double>& dcState) override;

    double calculateVoltageDrop(const std::vector<double>& state) const override;
    double calculateCurrent(const SimulationContext& simCtx,
                            MnaIndex extraIdx = MnaIndex{}) const override;

    void resetState() override;
    void saveCheckpoint() override;
    void restoreCheckpoint() override;

    bool isUsingOsdi() const { return osdiReady; }
    const BsimCmgParameters& getParameters() const { return params; }
    double getEffectiveWidth() const { return weff; }

    double getValue() const override { return params.VTH0; }
    void setValue(double newVth0) override;
};

