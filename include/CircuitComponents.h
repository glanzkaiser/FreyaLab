// CircuitComponents.h
#pragma once
#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <utility>
#include <complex>

struct CsrMatrix; 
struct CsrMatrixComplex;

// Represents a 1-based netlist node (0 = Ground)
struct NetNode {
    int id = 0;
    bool isGround() const { return id <= 0; }
    bool operator==(const NetNode& other) const { return id == other.id; }
    bool operator!=(const NetNode& other) const { return id != other.id; }
    bool operator<(const NetNode& other) const { return id < other.id; }
};

// Represents a 0-based linear solver row/column (-1 = Ground / invalid)
struct MnaIndex {
    int idx = -1;
    bool isValid() const { return idx >= 0; }
    static MnaIndex fromNetNode(NetNode n) { return MnaIndex{ n.id - 1 }; }
    static MnaIndex fromExtraRow(int rawIdx) { return MnaIndex{ rawIdx }; }
    bool operator==(const MnaIndex& other) const { return idx == other.idx; }
    bool operator!=(const MnaIndex& other) const { return idx != other.idx; }
    bool operator<(const MnaIndex& other) const { return idx < other.idx; }
};

// Bounds-checked voltage lookup (converts 1-based NetNode to 0-based MnaIndex safely)
inline double getSafeVoltage(const std::vector<double>& state, NetNode node) {
    if (node.isGround()) return 0.0;
    MnaIndex m = MnaIndex::fromNetNode(node);
    if (!m.isValid() || static_cast<size_t>(m.idx) >= state.size()) return 0.0;
    return state[m.idx];
}

inline double getSafeVoltage(const std::vector<double>* state, NetNode node) {
    if (!state) return 0.0;
    return getSafeVoltage(*state, node);
}

// Bounds-checked RHS vector injection using 0-based MnaIndex
inline void addSafeRhs(std::vector<double>& B, MnaIndex row, double val) {
    if (row.isValid() && static_cast<size_t>(row.idx) < B.size()) {
        B[row.idx] += val;
    }
}

inline void addSafeRhs(std::vector<std::complex<double>>& B_c, MnaIndex row, std::complex<double> val) {
    if (row.isValid() && static_cast<size_t>(row.idx) < B_c.size()) {
        B_c[row.idx] += val;
    }
}

// Bounds-checked RHS vector assignment using 0-based MnaIndex
inline void setSafeRhs(std::vector<double>& B, MnaIndex row, double val) {
    if (row.isValid() && static_cast<size_t>(row.idx) < B.size()) {
        B[row.idx] = val;
    }
}
// Bounds-checked complex RHS vector assignment using 0-based MnaIndex
inline void setSafeRhs(std::vector<std::complex<double>>& B_c, MnaIndex row, std::complex<double> val) {
    if (row.isValid() && static_cast<size_t>(row.idx) < B_c.size()) {
        B_c[row.idx] = val;
    }
}

// Bounds-checked read helper for 0-based extra MNA rows (e.g., branch currents)
inline double getSafeExtraValue(const std::vector<double>& state, MnaIndex row) {
    if (!row.isValid() || static_cast<size_t>(row.idx) >= state.size()) return 0.0;
    return state[row.idx];
}

inline double getSafeExtraValue(const std::vector<double>* state, MnaIndex row) {
    if (!state) return 0.0;
    return getSafeExtraValue(*state, row);
}

enum class SourceType { DC, AC, PULSE };

// High-level MNA Primitive Categories
enum class ComponentType { 
    RESISTOR, VOLTAGE_SOURCE, CURRENT_SOURCE, CAPACITOR, INDUCTOR, 
    DIODE, BJT, MOSFET, IGBT, DARLINGTON, FINFET, RIBBONFET, GAA_FET, 
    THYRISTOR_DEVICE, ANALOG_IC, DIGITAL_IC, AC_SOURCE_DEVICE, 
    PROTECTION_DEVICE, ELECTROMECHANICAL, PWM_CONTROLLER 
};

// Specialized Groupings for physics and behavioral dispatch
enum class ThyristorType { SCR, TRIAC, DIAC };
enum class AnalogIcType { OPAMP, VOLTAGE_REGULATOR, INSTRUMENTATION_AMP, TIMER_555 };
enum class DigitalIcType { LOGIC_GATE, MICROCONTROLLER, MICROPROCESSOR, DSP, FPGA };
enum class AcSourceType { ALTERNATOR, SIGNAL_GENERATOR, CRYSTAL_OSCILLATOR };
enum class ProtectionType { FUSE, CIRCUIT_BREAKER, TVS_DIODE, GDT };
enum class ElectromechanicalType { SWITCH, RELAY, CONNECTOR, TERMINAL_BLOCK };
enum class TransformerType { STEP_UP, STEP_DOWN, ISOLATION, AUTOTRANSFORMER, TOROIDAL };
enum class DiodeApplication { RECTIFIER, SIGNAL, SCHOTTKY, VARACTOR, TUNNEL, PIN };
enum class JfetPolarity { NJFET, PJFET };
enum class GateType { AND, OR, NAND, NOR, XOR, XNOR, NOT };
enum class BatteryChemistry { LEAD_ACID, LITHIUM_ION, NIMH };
enum class SignalShape { SINE, SQUARE, TRIANGLE, SAWTOOTH, HARMONIC_PULSE };


// TR-BDF2 Integrated State Matrix Context (Type-Safe MnaIndex Layer)
struct MnaMatrixContext {
    CsrMatrix& G;
    std::vector<double>& B;
    const std::vector<double>* currentIterV = nullptr;   // v_{k} (Current Newton-Raphson Guess)
    const std::vector<double>* historicStepV = nullptr;  // v_{n} (Converged Solution from Past Step)
    const std::vector<double>* midpointStepV = nullptr;  // v_{n+gamma} (Trapezoidal Midpoint Substep)
    double time = 0.0;
    double dt = 1e-6;
    int orderStep = 0;                        // 0 = Static DC, 1 = TR Midpoint, 2 = BDF2 Step Phase
    MnaIndex assignedExtraIdx{};              // Default invalid (-1)
};

struct SimulationContext {
    const std::vector<double>& state;
    double dt = 0.0;
    int orderStep = 0; // 0 = DC Operating Point, 1 = TR Midpoint, 2 = BDF2 Step Phase
};

// Base Component Interface
class Component {
protected:
    std::string name;
    ComponentType type;
    std::vector<NetNode> nodes;
    MnaIndex assignedExtraRowIdx{};

public:
    Component(const std::string& name, ComponentType type, const std::vector<NetNode>& nodes)
        : name(name), type(type), nodes(nodes) {}

    virtual ~Component() = default;

    std::string getName() const { return name; }
    ComponentType getType() const { return type; }
    const std::vector<NetNode>& getNodes() const { return nodes; }

    // Dynamic Extra MNA Row Topology Tracking
    virtual int getExtraRowRequirement() const { return 0; }
    virtual void allocateRows(int& globalExtraRowCounter) {
        int req = getExtraRowRequirement();
        if (req > 0) {
            assignedExtraRowIdx = MnaIndex::fromExtraRow(globalExtraRowCounter);
            globalExtraRowCounter += req;
        }
    }
    void setAssignedExtraIdx(MnaIndex idx) { assignedExtraRowIdx = idx; }
    MnaIndex getAssignedExtraIdx() const { return assignedExtraRowIdx; }

    // Numerical Solver Pipeline Hooks
    virtual void stamp(MnaMatrixContext& ctx) = 0;
    virtual void updateHistory(const MnaMatrixContext& ctx) {}
    virtual void getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns,
                                    MnaIndex extraIdx = MnaIndex{}) const = 0;
    virtual double calculateCurrent(const SimulationContext& simCtx,
                                    MnaIndex extraIdx = MnaIndex{}) const = 0;
    virtual double calculateVoltageDrop(const std::vector<double>& state) const = 0;

    // Phasor-Domain Small-Signal Multi-Frequency Stamping Interface
    virtual void stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                         double omega, const std::vector<double>& dcState) {
        (void)A_c; (void)B_c; (void)omega; (void)dcState;
    }

    // Checkpoint Rollback Pipeline Interfaces
    virtual void resetState() {}
    virtual void saveCheckpoint() {}
    virtual void restoreCheckpoint() {}

    virtual double getValue() const { return 0.0; }
    // Parametric Value Mutation Interface for DC Sweeps and Optimizers
    virtual void setValue(double newValue) { (void)newValue; }
};


// Primitives Declarations
class Resistor : public Component {
    double value;
public:
    Resistor(const std::string& name, const std::vector<NetNode>& nodes, double r);
    void stamp(MnaMatrixContext& ctx) override;
    void getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns,
                            MnaIndex extraIdx = MnaIndex{}) const override;
    double calculateVoltageDrop(const std::vector<double>& state) const override;
    double calculateCurrent(const SimulationContext& simCtx,
                            MnaIndex extraIdx = MnaIndex{}) const override;
    void stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                 double omega, const std::vector<double>& dcState) override;

    void setValue(double newValue) override { value = newValue; }
    double getValue() const override { return value; }
};

class Diode : public Component {
    double Is = 1e-14;
    double N = 1.0;
private:
    // Small-signal frequency domain AC depletion and carrier metrics
    double Cj0 = 2e-12;       // Zero-bias junction capacitance (2 pF)
    double Vj = 0.7;          // Junction potential (0.7 V)
    double M = 0.5;           // Grading coefficient (0.5 for abrupt junctions)
    double tt = 1e-9;         // Transit time / carrier lifetime (1 ns)
    double fc = 0.5;          // Forward-bias capacitance coefficient symmetry boundary

public:
    Diode(const std::string& name, const std::vector<NetNode>& nodes);
    void stamp(MnaMatrixContext& ctx) override;
    void getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns,
                            MnaIndex extraIdx = MnaIndex{}) const override;
    void stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                 double omega, const std::vector<double>& dcState) override;
    double calculateVoltageDrop(const std::vector<double>& state) const override;
    double calculateCurrent(const SimulationContext& simCtx,
                            MnaIndex extraIdx = MnaIndex{}) const override;
};

class Capacitor : public Component {
    double capacitance;
    double q_hist_n = 0.0;
    double q_hist_nc = 0.0;
    double q_hist_n_checkpoint = 0.0;
    double q_hist_nc_checkpoint = 0.0;
    double last_dt = 1e-3;
public:
    Capacitor(const std::string& name, const std::vector<NetNode>& nodes, double c);
    void stamp(MnaMatrixContext& ctx) override;
    void updateHistory(const MnaMatrixContext& ctx) override;
    void getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns,
                            MnaIndex extraIdx = MnaIndex{}) const override;
    double calculateVoltageDrop(const std::vector<double>& state) const override;
    double calculateCurrent(const SimulationContext& simCtx,
                            MnaIndex extraIdx = MnaIndex{}) const override;
    void stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                 double omega, const std::vector<double>& dcState) override;

    double getValue() const override { return capacitance; }
    void saveCheckpoint() override;
    void restoreCheckpoint() override;
};

class Inductor : public Component {
    double inductance;
    double i_hist_n = 0.0;
    double i_hist_nc = 0.0;
    double i_hist_n_checkpoint = 0.0;
    double i_hist_nc_checkpoint = 0.0;

public:
    Inductor(const std::string& name, const std::vector<NetNode>& nodes, double l);
    int getExtraRowRequirement() const override { return 1; }
    void allocateRows(int& globalExtraRowCounter) override;
    void stamp(MnaMatrixContext& ctx) override;
    void updateHistory(const MnaMatrixContext& ctx) override;
    void getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns,
                            MnaIndex extraIdx = MnaIndex{}) const override;
    double calculateVoltageDrop(const std::vector<double>& state) const override;
    double calculateCurrent(const SimulationContext& simCtx,
                            MnaIndex extraIdx = MnaIndex{}) const override;
    void stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                 double omega, const std::vector<double>& dcState) override;
    double getValue() const override { return inductance; }
    void saveCheckpoint() override;
    void restoreCheckpoint() override;
};

class Mosfet : public Component {
private:
    bool isNChannel;
    double vt0 = 0.7;
    double beta = 1e-3;
    double lambda = 0.02;

    // Small-signal frequency domain AC intrinsic and overlap parasitics
    double Cgs0 = 1e-12;       // Zero-bias Gate-Source overlap capacitance (1 pF)
    double Cgd0 = 0.5e-12;     // Zero-bias Gate-Drain overlap capacitance (0.5 pF)

public:
    Mosfet(const std::string& name, const std::vector<NetNode>& nodes, bool nChannel);

    // Core Numerical Solver Pipeline Hooks
    void stamp(MnaMatrixContext& ctx) override;
    void getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns,
                            MnaIndex extraIdx = MnaIndex{}) const override;

    // Multi-Frequency Small-Signal Linearized Phasor Overrides
    void stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                 double omega, const std::vector<double>& dcState) override;

    // Post-Facto Telemetry Probing and Waveform Trackers
    double calculateVoltageDrop(const std::vector<double>& state) const override;
    double calculateCurrent(const SimulationContext& simCtx,
                            MnaIndex extraIdx = MnaIndex{}) const override;
};

class VbicBjtAdvanced : public Component {
private:
    bool isNpn;
    double Is = 1e-16;
    double Nf = 1.0;
    double Nr = 1.0;
    double Ikf = 0.01;

    // Small-signal frequency domain dynamic base charge parasitics
    double Cje0 = 5e-12;       // Zero-bias Base-Emitter capacitance (5 pF)
    double Vje = 0.75;         // Base-Emitter junction potential
    double Me = 0.33;          // Base-Emitter grading coefficient
    double tf = 0.1e-9;        // Forward transit time (100 ps)

    double Cjc0 = 2e-12;       // Zero-bias Base-Collector capacitance (2 pF)
    double Vjc = 0.65;         // Base-Collector junction potential
    double Mc = 0.5;           // Base-Collector grading coefficient
    double tr = 10e-9;         // Reverse transit time (10 ns)

public:
    VbicBjtAdvanced(const std::string& name, const std::vector<NetNode>& nodes, bool npnVariant);

    // Core Numerical Solver Pipeline Hooks
    void stamp(MnaMatrixContext& ctx) override;
    void getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns,
                            MnaIndex extraIdx = MnaIndex{}) const override;

    // Multi-Frequency Small-Signal Linearized Phasor Overrides
    void stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                 double omega, const std::vector<double>& dcState) override;

    // Post-Facto Telemetry Probing and Waveform Trackers
    double calculateVoltageDrop(const std::vector<double>& state) const override;
    double calculateCurrent(const SimulationContext& simCtx,
                            MnaIndex extraIdx = MnaIndex{}) const override;
};

class VoltageSource : public Component {
    double value;
    double acMagnitude = 1.0;
    double acPhase = 0.0;
public:
    VoltageSource(const std::string& name, const std::vector<NetNode>& nodes, double v);
    int getExtraRowRequirement() const override { return 1; }
    void allocateRows(int& globalExtraRowCounter) override;
    void stamp(MnaMatrixContext& ctx) override;
    void getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns,
                            MnaIndex extraIdx = MnaIndex{}) const override;
    double calculateVoltageDrop(const std::vector<double>& state) const override;
    double calculateCurrent(const SimulationContext& simCtx,
                            MnaIndex extraIdx = MnaIndex{}) const override;
    void stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                 double omega, const std::vector<double>& dcState) override;
    double getValue() const override { return value; }
    void setAcParameters(double mag, double phase);
    void setValue(double newValue) override { value = newValue; }
};

class BehaviorLogicGateComponent : public Component {
private:
    GateType gate;
    double thresholdVoltage;
    double outputHighV;
    double outputLowV;
    double outputResistance;
public:
    BehaviorLogicGateComponent(const std::string& name, const std::vector<NetNode>& nodes, GateType gateType,
                               double vThresh = 2.5, double vHigh = 5.0, double vLow = 0.0, double rOut = 10.0);
    int getExtraRowRequirement() const override { return 1; }
    void allocateRows(int& globalExtraRowCounter) override;
    void stamp(MnaMatrixContext& ctx) override;
    void getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns,
                            MnaIndex extraIdx = MnaIndex{}) const override;
    double calculateVoltageDrop(const std::vector<double>& state) const override;
    double calculateCurrent(const SimulationContext& simCtx,
                            MnaIndex extraIdx = MnaIndex{}) const override;
    void stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                 double omega, const std::vector<double>& dcState) override;
};


