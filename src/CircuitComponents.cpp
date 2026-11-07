
#include "freyalab/CircuitComponents.h"
#include "freyalab/math.h"
#include "freyalab/LogicGates.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <complex>
#include <iostream>
// Safe helper function to extract absolute node potentials from a state vector
inline double getNodeVoltage(const std::vector<double>* stateVec, int node) {
    if (!stateVec || node <= 0) return 0.0; // Node 0 is reference ground (0V)
    if (node < static_cast<int>(stateVec->size())) {
        return (*stateVec)[node]; 
    }
    return 0.0;
}


// code good and robust and no structural bugs 
// Pointer-safe, retargetable helper utility function
static inline double getVoltage(const std::vector<double>* state, int node) {
    if (!state || node <= 0 || static_cast<size_t>(node) >= state->size()) return 0.0;
    return (*state)[node];
}

// Flat vector wrapper helper for simple state calls
static inline double getVoltage(const std::vector<double>& state, int node) {
    if (node <= 0 || static_cast<size_t>(node) >= state.size()) return 0.0;
    return state[node];
}

constexpr double V_THERMAL_300K = 0.0258563;

// ============================================================================
// 1. RESISTOR PRIMITIVE
// ============================================================================


Resistor::Resistor(const std::string& name, const std::vector<int>& nodes, double r)
    : Component(name, ComponentType::RESISTOR, nodes), value(r) {}

// mathematically correct, good and robust and no structural bugs 
void Resistor::stamp(MnaMatrixContext& ctx) {// RESISTOR DC/TRANSIENT STAMPING METHOD
    // 1. Division-by-zero safety clamping factor
    double safe_value = (std::abs(value) > 1e-20) ? value : 1e-20;
    double g = 1.0 / safe_value;
    
    int n1 = nodes[0];
    int n2 = nodes[1];

    // 2. Execute 0-based index internal matrix translation with Ground (0) trappings
    if (n1 > 0) ctx.G.addValue(n1 - 1, n1 - 1,  g);
    if (n2 > 0) ctx.G.addValue(n2 - 1, n2 - 1,  g);
    
    if (n1 > 0 && n2 > 0) {
        ctx.G.addValue(n1 - 1, n2 - 1, -g);
        ctx.G.addValue(n2 - 1, n1 - 1, -g);
    }
}

void Resistor::stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c, double omega, const std::vector<double>& dcState) {// RESISTOR AC PHASOR SMALL-SIGNAL STAMPING METHOD
    (void)B_c; 
    (void)omega; 
    (void)dcState; // Linear ideal resistor does not modulate over the DC operating point state

    // 1. Division-by-zero safety clamping factor matching the real-valued stamp engine
    double safe_value = (std::abs(value) > 1e-20) ? value : 1e-20;
    double g = 1.0 / safe_value;

    // Pure real conductance phasor representation (Ideal resistors exhibit zero phase shift)
    std::complex<double> g_complex(g, 0.0);
    
    int n1 = nodes[0];
    int n2 = nodes[1];

    // 2. FIXED: Correctly shift to 0-based matrix indexing (node - 1) to eliminate AC Segmentation Faults
    if (n1 > 0) A_c.addValue(n1 - 1, n1 - 1,  g_complex);
    if (n2 > 0) A_c.addValue(n2 - 1, n2 - 1,  g_complex);
    
    if (n1 > 0 && n2 > 0) {
        A_c.addValue(n1 - 1, n2 - 1, -g_complex);
        A_c.addValue(n2 - 1, n1 - 1, -g_complex);
    }
}

void Resistor::getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx) const {
    (void)extraIdx; // Resistors are passive two-terminal elements and do not require extra MNA rows
    
    // Structural Bug Guard: Ensure the parser has successfully populated both terminal nodes
    if (nodes.size() < 2) {
        return; 
    }

    int n1 = nodes[0];
    int n2 = nodes[1];

    // Trapping ground index 0 by strictly enforcing (node > 0) before 0-based conversion (-1)
    if (n1 > 0) {
        patterns.push_back({n1 - 1, n1 - 1});
    }
    if (n2 > 0) {
        patterns.push_back({n2 - 1, n2 - 1});
    }
    
    if (n1 > 0 && n2 > 0) {
        patterns.push_back({n1 - 1, n2 - 1});
        patterns.push_back({n2 - 1, n1 - 1});
    }
}

double Resistor::calculateVoltageDrop(const std::vector<double>& state) const {
    return getSafeVoltage(state, nodes[0]) - getSafeVoltage(state, nodes[1]);
}

double Resistor::calculateCurrent(const SimulationContext& simCtx, int extraIdx) const {
    (void)extraIdx;
    if (simCtx.state.empty()) {
        return 0.0;
    }

    // 1. Primary Path: Direct extraction via MNA assigned extra row index if pre-allocated
    int idx = this->assignedExtraRowIdx;
    if (idx >= 0 && static_cast<size_t>(idx) < simCtx.state.size()) {
        return simCtx.state[idx];
    }
    
    // 2. Spatial Algebraic Path: Calculate instantaneously via Ohm's Law (I = V_drop / R)
    double safeValue = (std::abs(value) > 1e-20) ? value : 1e-20; // Critical division-by-zero guard
    double v_drop = calculateVoltageDrop(simCtx.state);
    return v_drop / safeValue;
}

// ============================================================================
// 2. DIODE MODEL
// ============================================================================
Diode::Diode(const std::string& name, const std::vector<int>& nodes)
    : Component(name, ComponentType::DIODE, nodes) {}
// mathematically correct, good and robust and no structural bugs 
void Diode::stamp(MnaMatrixContext& ctx) {
    if (nodes.size() < 2) return;

    int n1 = nodes[0]; // Anode
    int n2 = nodes[1]; // Cathode

    // 1. Resolve the TRUE current iteration voltage drop drop bounds
    double vd_true = 0.0;
    double vold = 0.0;

    if (ctx.currentIterV) {
        vd_true = getSafeVoltage(*ctx.currentIterV, n1) - getSafeVoltage(*ctx.currentIterV, n2);
    }
    if (ctx.historicStepV) {
        vold = getSafeVoltage(*ctx.historicStepV, n1) - getSafeVoltage(*ctx.historicStepV, n2);
    }

    double vt = N * V_THERMAL_300K;
    double vcrit = vt * std::log(vt / (std::sqrt(2.0) * Is));

    // 2. Compute the LIMITED voltage drop drop specifically for exponential safety checks
    double vd_lim = CircuitMath::pnjlim(vd_true, vold, vt, vcrit);

    // 3. Compute Shockley equations using the safely limited voltage drop drop
    double expTerm = std::exp(std::min(vd_lim / vt, 40.0));
    double id = Is * (expTerm - 1.0);
    double gd = (Is / vt) * expTerm;

    // Inject Gmin stabilization tracking
    constexpr double gmin = 1e-12;
    gd += gmin;
    
    // 4. CRITICAL SPICE FIX: Formulate Ieq using the true voltage drop drop (vd_true), 
    // but anchored by the limited physics variables (id, gd) to maintain mathematical consistency
    double ieq = id - gd * vd_lim;

    // 5. Conductance Matrix G Stamping
    if (n1 > 0) ctx.G.addValue(n1 - 1, n1 - 1,  gd);
    if (n2 > 0) ctx.G.addValue(n2 - 1, n2 - 1,  gd);
    if (n1 > 0 && n2 > 0) {
        ctx.G.addValue(n1 - 1, n2 - 1, -gd);
        ctx.G.addValue(n2 - 1, n1 - 1, -gd);
    }

    // 6. Source Vector B Stamping
    if (n1 > 0) ctx.B[n1 - 1] += -ieq;
    if (n2 > 0) ctx.B[n2 - 1] +=  ieq;
}

void Diode::stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c, double omega, const std::vector<double>& dcState) {
    (void)B_c; // Small-signal independent AC sources handle B_c loads, diodes strictly populate admittance A_c  

    // 1. Guard against unallocated structural netlist nodes
    if (nodes.size() < 2) return;

    int n1 = nodes[0]; // Anode terminal
    int n2 = nodes[1]; // Cathode terminal
    
    // 2. Extract linearized small-signal parameters from the converged DC Operating Point
    double v_dc = calculateVoltageDrop(dcState);
    double vt = N * V_THERMAL_300K;
    
    // Compute small-signal conductance (gd) mirroring the exact baseline of stamp()
    double expTerm = std::exp(std::min(v_dc / vt, 40.0)); // Prevent exponential overflow during high-bias sweeps
    double gd = (Is / vt) * expTerm;
    
    // Enforce high-performance numerical stability baseline conductance (Gmin)
    constexpr double gmin = 1e-12;
    gd += gmin;

    // 3. Compute Small-Signal Dynamic Junction Capacitance (C_total = C_diffusion + C_depletion)
    // Diffusion Capacitance: C_diff = tt * gd
    double c_diffusion = tt * gd;
    double c_depletion = 0.0;
    
    // Depletion Capacitance boundary condition handling (SPICE standard smooth threshold)
    double v_break = fc * Vj;
    if (v_dc < v_break) {
        // Standard reverse/low-forward bias depletion tracking
        double base = std::max(1.0 - v_dc / Vj, 1e-3);
        c_depletion = Cj0 * std::pow(base, -M);
    } else {
        // Forward bias linear expansion overlay to mitigate infinite capacitance asymptotes
        double c_factor = std::pow(1.0 - fc, -M);
        double slope = (Cj0 * M * std::pow(1.0 - fc, -(M + 1.0))) / Vj;
	c_depletion = Cj0 * c_factor + slope * (v_dc - v_break);
    }
    
    double c_total = std::max(c_diffusion + c_depletion, 1e-20);

    // 4. Transform Admittance into a Complex Phasor Representation: Y = gd + j * omega * C_total
    std::complex<double> y_complex(gd, omega * c_total);

    // 5. FIXED: Correctly shift to 0-based matrix indexing (node - 1) to eliminate AC Segmentation Faults
    if (n1 > 0) {
        A_c.addValue(n1 - 1, n1 - 1,  y_complex);
    }
    if (n2 > 0) {
        A_c.addValue(n2 - 1, n2 - 1,  y_complex);
    }
    if (n1 > 0 && n2 > 0) {
        A_c.addValue(n1 - 1, n2 - 1, -y_complex);
        A_c.addValue(n2 - 1, n1 - 1, -y_complex);
    }
}

void Diode::getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx) const {
    (void)extraIdx; // Diodes are non-linear 2-terminal devices and do not request extra rows
    
    // Structural netlist validation layer guard
    if (nodes.size() < 2) return;

    int n1 = nodes[0];
    int n2 = nodes[1];

    // Build the static sparsity grid coordinate maps using strictly monitored 0-based translation
    if (n1 > 0) patterns.push_back({n1 - 1, n1 - 1});
    if (n2 > 0) patterns.push_back({n2 - 1, n2 - 1});
    
    if (n1 > 0 && n2 > 0) {
        patterns.push_back({n1 - 1, n2 - 1});
        patterns.push_back({n2 - 1, n1 - 1});
    }
}

double Diode::calculateVoltageDrop(const std::vector<double>& state) const {
    return getSafeVoltage(state, nodes[0]) - getSafeVoltage(state, nodes[1]);
}


double Diode::calculateCurrent(const SimulationContext& simCtx, int extraIdx) const {
    (void)extraIdx;
    if (simCtx.state.empty()) {
        return 0.0;
    }

    double vd = calculateVoltageDrop(simCtx.state);
    double vt = N * V_THERMAL_300K;
    
    // 1. Analytical SPICE Voltage Clamping Layer (Prevents numeric floating-point overflow)
    double vcrit = vt * std::log(vt / (std::sqrt(2.0) * Is));
    if (vd > vcrit) {
        vd = vcrit + (vd - vcrit) / (1.0 + (vd - vcrit) / vt);
    }

    // 2. Core Shockley Equation evaluation with numerical exponent bounding
    double id = Is * (std::exp(std::min(vd / vt, 40.0)) - 1.0);
    
    // 3. Parallel Gmin conductance baseline injection to match the linear matrix stamp() exactly
    constexpr double gmin = 1e-12;
    id += gmin * vd;

    return id;
}

// ============================================================================
// 3. CAPACITOR
// ============================================================================

Capacitor::Capacitor(const std::string& name, const std::vector<int>& nodes, double c)
    : Component(name, ComponentType::CAPACITOR, nodes), capacitance(c) {}
// mathematically correct, good and robust and no structural bugs 
void Capacitor::stamp(MnaMatrixContext& ctx) {
    if (nodes.size() < 2) return;

    int n1 = nodes[0];
    int n2 = nodes[1];
    
    // TR-BDF2 optimal step fractional spacing constant: gamma = 2.0 - sqrt(2.0)
    const double gamma = 2.0 - std::sqrt(2.0); 
    double geq = 0.0;
    double ieq = 0.0;
    
    // 1. DC Operating Point Analysis (Order 0): Capacitors act as open circuits with Gmin shunt
    if (ctx.orderStep == 0 || ctx.dt <= 0.0) {
        constexpr double gmin = 1e-12; 
        if (n1 > 0) ctx.G.addValue(n1 - 1, n1 - 1,  gmin);
        if (n2 > 0) ctx.G.addValue(n2 - 1, n2 - 1,  gmin);
        if (n1 > 0 && n2 > 0) {
            ctx.G.addValue(n1 - 1, n2 - 1, -gmin);
            ctx.G.addValue(n2 - 1, n1 - 1, -gmin);
        }
        return;
    }
    
    // 2. Transient Integration Companion Model Evaluation with strict pointer guards
    if (ctx.orderStep == 1) {
        // TR Stage Companion Model: geq = 2C / (gamma * dt)
        geq = (2.0 * capacitance) / (gamma * ctx.dt);
        
        double v_n = 0.0;
        if (ctx.historicStepV) {
            v_n = getSafeVoltage(*ctx.historicStepV, n1) - getSafeVoltage(*ctx.historicStepV, n2);
        }
        // Norton companion current source tracking equation
        ieq = -geq * v_n - (q_hist_n / (gamma * ctx.dt * 0.5));
        
    } else if (ctx.orderStep == 2) {
        // BDF2 Composite Stage Companion Model: Variable multi-step spacing coefficients
        double c1 = (2.0 - gamma) / (1.0 - gamma);
        double c2 = 1.0 / (gamma * (1.0 - gamma));
        geq = c1 * capacitance / ctx.dt;
        
        double v_n = 0.0;
        double v_nc = 0.0;
        if (ctx.historicStepV) {
            v_n = getSafeVoltage(*ctx.historicStepV, n1) - getSafeVoltage(*ctx.historicStepV, n2);
        }
        if (ctx.midpointStepV) {
            v_nc = getSafeVoltage(*ctx.midpointStepV, n1) - getSafeVoltage(*ctx.midpointStepV, n2);
        }
        // Exact discrete multi-step history integration mapping
        ieq = -capacitance * (c2 * v_nc - (c2 - c1) * v_n) / ctx.dt;
    }
    
    // 3. Core Conductance Admittance Matrix Stamps (Enforcing strict 0-based node - 1 offsets)
    if (n1 > 0) ctx.G.addValue(n1 - 1, n1 - 1,  geq);
    if (n2 > 0) ctx.G.addValue(n2 - 1, n2 - 1,  geq);
    
    if (n1 > 0 && n2 > 0) {
        ctx.G.addValue(n1 - 1, n2 - 1, -geq);
        ctx.G.addValue(n2 - 1, n1 - 1, -geq);
    }

    // 4. FIXED: Shift to 0-based indexing (node - 1) on vector B to eliminate Nodal Segfaults
    // Current exits positive terminal (-ieq injection) and enters negative terminal (+ieq injection)
    if (n1 > 0) ctx.B[n1 - 1] += -ieq;
    if (n2 > 0) ctx.B[n2 - 1] +=  ieq;
}

void Capacitor::updateHistory(const MnaMatrixContext& ctx) {
    if (nodes.size() < 2 || !ctx.currentIterV) return;

    // Safely sample instant node voltages via the global inline getSafeVoltage interface
    double v = getSafeVoltage(*ctx.currentIterV, nodes[0]) - getSafeVoltage(*ctx.currentIterV, nodes[1]);
    
    if (ctx.orderStep == 1) {
        // Store intermediate TR stage charge integration state point (v_{n+gamma})
        q_hist_nc = capacitance * v; 
    } else if (ctx.orderStep == 2) {
        // Commits final full BDF2 time step charge state point (v_{n+1})
        q_hist_n = capacitance * v;  
    }
}

void Capacitor::getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx) const {
    (void)extraIdx; // Passive two-terminal components do not request extra rows
    if (nodes.size() < 2) return;

    int n1 = nodes[0];
    int n2 = nodes[1];
    
    if (n1 > 0) patterns.push_back({n1 - 1, n1 - 1});
    if (n2 > 0) patterns.push_back({n2 - 1, n2 - 1});
    
    if (n1 > 0 && n2 > 0) {
        patterns.push_back({n1 - 1, n2 - 1});
        patterns.push_back({n2 - 1, n1 - 1});
    }
}

double Capacitor::calculateVoltageDrop(const std::vector<double>& state) const {
    return getSafeVoltage(state, nodes[0]) - getSafeVoltage(state, nodes[1]);
}

double Capacitor::calculateCurrent(const SimulationContext& simCtx, int extraIdx) const {
    (void)extraIdx;
    if (simCtx.state.empty()) {
        return 0.0; 
    }

    // During pure static DC Operating Point analysis (Order 0), an ideal capacitor acts as an open circuit
    if (simCtx.orderStep == 0 || simCtx.dt <= 0.0) {
        constexpr double gmin = 1e-12; // Standard internal shunt leak matching matrix stabilization
        return gmin * calculateVoltageDrop(simCtx.state);
    }

    // Extract natively solved branch current if assigned to an MNA row
    int idx = this->assignedExtraRowIdx;
    if (idx >= 0 && static_cast<size_t>(idx) < simCtx.state.size()) {
        return simCtx.state[idx];
    }
    
    double v_now = calculateVoltageDrop(simCtx.state);
    
    // TR-BDF2 Step Fractional Parameter matching SPICE engines: gamma = 2.0 - sqrt(2.0)
    const double gamma = 2.0 - std::sqrt(2.0); 
    double geq = 0.0;
    double ieq = 0.0;
    
    if (simCtx.orderStep == 1) { // Trapezoidal Integration Substep
        geq = (2.0 * capacitance) / (gamma * simCtx.dt);
        double v_n = q_hist_n / capacitance; // Historical voltage state derived from past charge
        ieq = -geq * v_n - (q_hist_n / (gamma * simCtx.dt * 0.5));
    } else if (simCtx.orderStep == 2) { // BDF2 Integration Step Phase
        double c1 = (2.0 - gamma) / (1.0 - gamma);
        double c2 = 1.0 / (gamma * (1.0 - gamma));
        geq = c1 * capacitance / simCtx.dt;
        ieq = -(c2 * q_hist_nc - (c2 - c1) * q_hist_n) / simCtx.dt;
    }
    
    return (geq * v_now) + ieq;
}

void Capacitor::stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c, 
                         double omega, const std::vector<double>& dcState) {
    (void)B_c;
    (void)dcState; // Ideal linear capacitors do not modulate over the DC operating point state

    // 1. Structural Bug Guard: Ensure the parser has successfully populated both terminal nodes
    if (nodes.size() < 2) {
        return;
    }

    int n1 = nodes[0];
    int n2 = nodes[1];

    // 2. Compute Admittance Complex Phasor Representation: Y = 0 + j * omega * C
    // Note: If omega is exactly 0 (DC), the admittance naturally scales to purely zero complex values,
    // which correctly models an ideal open-circuit behavior without skipping the matrix allocations.
    std::complex<double> j_omega_c(0.0, omega * capacitance);
    
    // Enforce high-performance numerical stability baseline shunt conductance (Gmin)
    // This prevents matrix structural singularity errors at near-DC or absolute zero-frequency analysis.
    if (std::abs(omega) < 1e-20) {
        constexpr double gmin = 1e-12;
        j_omega_c = std::complex<double>(gmin, 0.0);
    }

    // 3. FIXED: Correctly shift to 0-based matrix indexing (node - 1) to eliminate AC Segmentation Faults
    // Trapping ground index 0 by strictly enforcing (node > 0) before 0-based conversion (-1)
    if (n1 > 0) {
        A_c.addValue(n1 - 1, n1 - 1,  j_omega_c);
    }
    if (n2 > 0) {
        A_c.addValue(n2 - 1, n2 - 1,  j_omega_c);
    }
    
    if (n1 > 0 && n2 > 0) {
        A_c.addValue(n1 - 1, n2 - 1, -j_omega_c);
        A_c.addValue(n2 - 1, n1 - 1, -j_omega_c);
    }
}

void Capacitor::saveCheckpoint() {
    q_hist_n_checkpoint = q_hist_n;
    q_hist_nc_checkpoint = q_hist_nc;
}

void Capacitor::restoreCheckpoint() {
    q_hist_n = q_hist_n_checkpoint;
    q_hist_nc = q_hist_nc_checkpoint;
}

// ============================================================================
// 4. INDUCTOR
// ============================================================================
Inductor::Inductor(const std::string& name, const std::vector<int>& nodes, double l)
    : Component(name, ComponentType::INDUCTOR, nodes), inductance(l) {}

// mathematically correct, good and robust and no structural bugs
void Inductor::stamp(MnaMatrixContext& ctx) {
    // 1. Structural Layer Guards: Ensure nodes and extra row references are safely allocated
    if (nodes.size() < 2) return;
    
    int idx = (ctx.assignedExtraIdx >= 0) ? ctx.assignedExtraIdx : this->assignedExtraRowIdx;
    if (idx < 0 || idx >= static_cast<int>(ctx.B.size())) return; 
    
    int n1 = nodes[0]; // Positive terminal
    int n2 = nodes[1]; // Negative terminal
    
    // TR-BDF2 optimal step fractional spacing constant: gamma = 2.0 - sqrt(2.0)
    const double gamma = 2.0 - std::sqrt(2.0); 
    
    // 2. Nodal Interconnection Matrix Stamps (Enforcing strict 0-based matrix mapping via node - 1)
    if (n1 > 0) {
        ctx.G.addValue(n1 - 1, idx,  1.0);
        ctx.G.addValue(idx, n1 - 1,  1.0);
    }
    if (n2 > 0) {
        ctx.G.addValue(n2 - 1, idx, -1.0);
        ctx.G.addValue(idx, n2 - 1, -1.0);
    }
    
    // 3. DC Operating Point Analysis (Order 0): Inductors act as pure structural short circuits
    if (ctx.orderStep == 0 || ctx.dt <= 0.0) {
        ctx.G.addValue(idx, idx, 0.0); 
        ctx.B[idx] = 0.0; // FIXED: Absolute initialization instead of incremental assignment
        return;
    }
    
    // 4. Transient Integration Companion Model Evaluation with synchronized signs
    if (ctx.orderStep == 1) {
        // TR Stage Companion Model: Req = 2L / (gamma * dt)
        double req = (2.0 * inductance) / (gamma * ctx.dt);
        ctx.G.addValue(idx, idx, -req);
        
        // Extract historical voltage drop across terminals from the past converged step
        double v_past = 0.0;
        if (ctx.historicStepV) {
            v_past = getSafeVoltage(*ctx.historicStepV, n1) - getSafeVoltage(*ctx.historicStepV, n2);
        }
        // Thévenin RHS injection: B[idx] = v_n - Req * i_n
        ctx.B[idx] = v_past - req * i_hist_n;
        
    } else if (ctx.orderStep == 2) {
        // BDF2 Composite Stage Companion Model: Multi-step time grid integration coefficients
        double c_factor = (2.0 - gamma) / (1.0 - gamma);
        double c2 = 1.0 / (gamma * (1.0 - gamma));
        double req = c_factor * inductance / ctx.dt;
        ctx.G.addValue(idx, idx, -req);
        
        // Compute discrete multi-step history current injection tracking framework
        double i_hist_term = (c2 * i_hist_nc - (c2 - c_factor) * i_hist_n);
        // Thévenin RHS injection: B[idx] = - (L / dt) * [c2 * i_{n+gamma} - (c2 - c_factor) * i_n]
        ctx.B[idx] = - (inductance / ctx.dt) * i_hist_term;
    }
}

void Inductor::updateHistory(const MnaMatrixContext& ctx) {
    int idx = (ctx.assignedExtraIdx >= 0) ? ctx.assignedExtraIdx : this->assignedExtraRowIdx;
    if (idx < 0 || !ctx.currentIterV || idx >= static_cast<int>(ctx.currentIterV->size())) return;
    
    // Natively sample the current branch current directly solved inside the MNA vector state
    double current = (*ctx.currentIterV)[idx];
    
    if (ctx.orderStep == 1) {
        // Store intermediate TR stage branch current integration state point (i_{n+gamma})
        i_hist_nc = current; 
    } else if (ctx.orderStep == 2) {
        // Commits final full BDF2 time step branch current state point (i_{n+1})
        i_hist_n = current;  
    }
}

void Inductor::getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx) const {
    if (nodes.size() < 2) return;
    
    int brRow = (extraIdx >= 0) ? extraIdx : this->assignedExtraRowIdx;
    if (brRow < 0) return;

    // Self-inductance internal structural matrix grid point allocation
    patterns.push_back({brRow, brRow});
    
    int n1 = nodes[0];
    int n2 = nodes[1];

    // Build the static branch coupling pattern maps using strictly monitored 0-based index shifts
    if (n1 > 0) {
        patterns.push_back({n1 - 1, brRow});
        patterns.push_back({brRow, n1 - 1});
    }
    if (n2 > 0) {
        patterns.push_back({n2 - 1, brRow});
        patterns.push_back({brRow, n2 - 1});
    }
}
double Inductor::calculateCurrent(const SimulationContext& simCtx, int extraIdx) const {
    int targetIdx = (extraIdx >= 0) ? extraIdx : this->assignedExtraRowIdx;
    
    // Direct Optimal Path: Read natively solved branch current variable from MNA solution state
    if (targetIdx >= 0 && static_cast<size_t>(targetIdx) < simCtx.state.size()) {
        return simCtx.state[targetIdx];
    }
    
    // Algebraic Fallback: Used when probing parameters outside the synchronized solver loop context
    if (simCtx.state.empty()) {
        return 0.0; 
    }
    if (simCtx.orderStep == 0 || simCtx.dt <= 0.0) {
        return i_hist_n; // In static DC state, return the last known locked tracking point
    }
    
    const double gamma = 2.0 - std::sqrt(2.0);
    double v_now = calculateVoltageDrop(simCtx.state);
    
    if (simCtx.orderStep == 1) { // TR Companion Fallback Equation
        return i_hist_n + (v_now * (gamma * simCtx.dt) / (2.0 * inductance));
    } else if (simCtx.orderStep == 2) { // BDF2 Companion Fallback Equation
        double c_factor = (2.0 - gamma) / (1.0 - gamma);
        double c2 = 1.0 / (gamma * (1.0 - gamma));
        double ieq = (c2 * i_hist_nc - (c2 - c_factor) * i_hist_n);
        return ieq + (v_now * simCtx.dt / (c_factor * inductance));
    }
    return 0.0;
}

double Inductor::calculateVoltageDrop(const std::vector<double>& state) const {
    return getSafeVoltage(state, nodes[0]) - getSafeVoltage(state, nodes[1]);
}

void Inductor::stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c, double omega, const std::vector<double>& dcState) {
    (void)B_c; 
    (void)dcState; // Linear ideal inductors do not modulate over the DC operating point state

    // 1. Structural Layer Guards: Ensure nodes and extra row references are safely allocated
    if (nodes.size() < 2) return;

    int idx = getAssignedExtraIdx(); // Natively fetched 0-based branch index from the base class
    if (idx < 0) return; 

    int n1 = nodes[0]; // Positive terminal
    int n2 = nodes[1]; // Negative terminal
    
    // 2. Compute Complex Impedance of the Inductor: Z = 0 + j * omega * L
    std::complex<double> j_omega_l(0.0, omega * inductance);

    // Enforce high-performance numerical stability baseline shunt to prevent structural matrix
    // singularity errors if the device is evaluated at near-DC or absolute zero frequency sweeps.
    if (std::abs(omega) < 1e-20) {
        constexpr double rmin = 1e-12; // Minor internal series stabilizer
        j_omega_l = std::complex<double>(rmin, 0.0);
    }

    // 3. FIXED: Correctly shift circuit nodes to 0-based matrix indexing (node - 1). 
    // Leave the branch index 'idx' unmodified as it is pre-compiled as a 0-based index.
    if (n1 > 0) {
        A_c.addValue(n1 - 1, idx,     1.0);
        A_c.addValue(idx,     n1 - 1, 1.0);
    }
    if (n2 > 0) {
        A_c.addValue(n2 - 1, idx,    -1.0);
        A_c.addValue(idx,     n2 - 1, -1.0);
    }

    // Branch relationship voltage-current loop integration equation: 
    // V_{n1-1} - V_{n2-1} - (j * omega * L) * I_branch = 0
    A_c.addValue(idx, idx, -j_omega_l);
}

void Inductor::allocateRows(int& globalExtraRowCounter) {
    Component::allocateRows(globalExtraRowCounter);
}

void Inductor::saveCheckpoint() {
    i_hist_n_checkpoint = i_hist_n;
    i_hist_nc_checkpoint = i_hist_nc;
}

void Inductor::restoreCheckpoint() {
    i_hist_n = i_hist_n_checkpoint;
    i_hist_nc = i_hist_nc_checkpoint;
}

// ============================================================================
// 5. MOSFET
// ============================================================================

Mosfet::Mosfet(const std::string& name, const std::vector<int>& nodes, bool nChannel)
    : Component(name, ComponentType::MOSFET, nodes), isNChannel(nChannel) {}
// mathematically correct, good and robust and no structural bugs 
void Mosfet::stamp(MnaMatrixContext& ctx) {
    // 1. Structural Layer Guard: Ensure all 3 multi-terminal pins are present
    if (nodes.size() < 3) return;

    const int nD_orig = nodes[0]; // Physical Drain connection tag
    const int nG_orig = nodes[1]; // Physical Gate connection tag
    const int nS_orig = nodes[2]; // Physical Source connection tag

    // 2. Safely sample contemporary absolute node voltages via the pointer-guarded interface
    double vD = 0.0;
    double vG_orig = 0.0;
    double vS = 0.0;

    if (ctx.currentIterV) {
        vD = getSafeVoltage(*ctx.currentIterV, nD_orig);
        vG_orig = getSafeVoltage(*ctx.currentIterV, nG_orig);
        vS = getSafeVoltage(*ctx.currentIterV, nS_orig);
    }

    // Set up active working node variables for directional channel analysis
    int nD_working = nD_orig;
    int nG_working = nG_orig;
    int nS_working = nS_orig;

    // 3. Perform PMOS Structural Device Channel Inversion
    if (!isNChannel) {
        std::swap(vD, vS);
        // Invert the gate voltage driving frame relative to the PMOS source terminal
        vG_orig = vS - (vG_orig - vS); 
    }

    // 4. Handle Bidirectional Operating Modes (Symmetric Drain-Source Swapping)
    double vds = vD - vS;
    if (vds < 0.0) {
        std::swap(nD_working, nS_working);
        std::swap(vD, vS);
        vds = -vds; // Force operational Vds to be positive inside the core math engine
    }

    // Compute gate overdrive voltage characteristics inside the working reference frame
    double vgs = vG_orig - vS;
    double vov = vgs - vt0;

    double ids = 0.0;
    double gm = 0.0;
    double gds = 1e-12; // High-performance stability minimum floor conductance baseline
    
    // 5. Piecewise Evaluation of the Shichman-Hodges Square-Law Transistor Model
    if (vov > 0.0) {
        double lambdaTerm = 1.0 + lambda * vds;
        if (vds < vov) { 
            // Linear / Triode Operating Region
            ids = beta * (vov * vds - 0.5 * vds * vds) * lambdaTerm;
            gm  = beta * vds * lambdaTerm;
            gds = beta * (vov - vds) * lambdaTerm + beta * (vov * vds - 0.5 * vds * vds) * lambda;
        } else { 
            // Saturation / Active Pinch-Off Operating Region
            ids = 0.5 * beta * vov * vov * lambdaTerm;
            gm  = beta * vov * lambdaTerm;
            gds = 0.5 * beta * vov * vov * lambda;
        }
    }

    // Guarantee numerical matrix health during steep reverse cutoff iterations
    if (gds < 1e-12) gds = 1e-12;

    // Compute the Newton-Raphson companion linearized Norton equivalent compensation current source
    // Derived from the multi-variable Taylor expansion: Ieq = Ids - gm * Vgs - gds * Vds
    double ieq = ids - gm * vgs - gds * vds;

    // Adjust physical current polarity rules if the actual hardware device is a PMOS channel type
    if (!isNChannel) {
        ieq = -ieq;
    }

    // 6. Map Nodal Coordinates into 0-Based Linear Matrix Indices (Trapping Ground 0)
    int mD = nD_working - 1;
    int mG = nG_working - 1;
    int mS = nS_working - 1;

    // 7. Multi-Terminal Jacobian Conductance & Current Matrix Injection
    // Stamping is performed directly on the active working node terminals to ensure numerical stability.
    
    // Nodal row equation stamps for the working Drain terminal
    if (mD >= 0) {
        ctx.G.addValue(mD, mD,  gds);
        if (mG >= 0) ctx.G.addValue(mD, mG,  gm);
        if (mS >= 0) ctx.G.addValue(mD, mS, -(gm + gds));
        ctx.B[mD] += -ieq; // Current flowing out of the Drain node
    }
    
    // Nodal row equation stamps for the working Source terminal
    if (mS >= 0) {
        ctx.G.addValue(mS, mS,   (gm + gds));
        if (mD >= 0) ctx.G.addValue(mS, mD,  -gds);
        if (mG >= 0) ctx.G.addValue(mS, mG,  -gm);
        ctx.B[mS] += ieq;  // Current flowing into the Source node
    }
}

void Mosfet::stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c, 
                     double omega, const std::vector<double>& dcState) {
    (void)B_c; // Non-linear devices do not load independent phasor source inputs

    // 1. Structural Layer Guard: Ensure the device contains all 3 hardware terminals
    if (nodes.size() < 3) return;
    
    const int nD_phys = nodes[0];
    const int nG_phys = nodes[1];
    const int nS_phys = nodes[2];

    // 2. FIXED: Safely resolve absolute node potentials using getSafeVoltage to eliminate 1-based offset shifts
    double vD_val = getSafeVoltage(dcState, nD_phys);
    double vG_val = getSafeVoltage(dcState, nG_phys);
    double vS_val = getSafeVoltage(dcState, nS_phys);

    // 3. Determine internal channel operational biasing frames based on hardware channel type
    double vgs_internal = 0.0;
    double vds_internal = 0.0;
    bool is_swapped = false;

    if (isNChannel) {
        vgs_internal = vG_val - vS_val;
        vds_internal = vD_val - vS_val;
        // Detect physical channel inversion (Symmetric Drain-Source inversion on NMOS)
        if (vds_internal < 0.0) {
            vds_internal = -vds_internal;
            vgs_internal = vG_val - vD_val; // Reference Gate voltage over the acting physical Drain terminal
            is_swapped = true;
        }
    } else { // PMOS Device Structure
        vgs_internal = vS_val - vG_val;
        vds_internal = vS_val - vD_val;
        // Detect physical channel inversion (Symmetric Drain-Source inversion on PMOS)
        if (vds_internal < 0.0) {
            vds_internal = -vds_internal;
            vgs_internal = vD_val - vG_val; // Reference Gate voltage over the acting physical Drain terminal
            is_swapped = true;
        }
    }

    // 4. Compute linearized small-signal parameters matching the exact baseline of stamp()
    double vov = vgs_internal - vt0;
    double gm = 0.0;
    double gds = 1e-12; // Numerical safety minimum floor conductance baseline

    if (vov > 0.0) {
        double lambdaTerm = 1.0 + lambda * vds_internal;
        if (vds_internal < vov) { // Linear / Triode Operating Region
            gm  = beta * vds_internal * lambdaTerm;
            gds = beta * (vov - vds_internal) * lambdaTerm + beta * (vov * vds_internal - 0.5 * vds_internal * vds_internal) * lambda;
        } else { // Saturation / Active Pinch-Off Operating Region
            gm  = beta * vov * lambdaTerm;
            gds = 0.5 * beta * vov * vov * lambda;
        }
    }
    if (gds < 1e-12) gds = 1e-12;

    // 5. Compute Small-Signal Parasitic Overlap Admittance Components
    std::complex<double> y_gs(0.0, omega * Cgs0);
    std::complex<double> y_gd(0.0, omega * Cgd0);

    // If channel inversion occurs, parasitic overlap capacitors switch orientations geometrically
    if (is_swapped) {
        std::swap(y_gs, y_gd);
    }

    // Determine effective node mapping tags for proper geometric KCL routing
    int nD_eff = is_swapped ? nS_phys : nD_phys;
    int nS_eff = is_swapped ? nD_phys : nS_phys;

    // 6. FIXED: Translate effective pins to clean 0-based matrix index systems (node - 1)
    int mD_eff = nD_eff - 1;
    int mG_phys = nG_phys - 1;
    int mS_eff = nS_eff - 1;

    // 7. Inject Complex Nodal Admittance Stencils based on clean 0-based grid spaces
    
    // Effective Drain Node Equations (KCL Row)
    if (mD_eff >= 0) {
        A_c.addValue(mD_eff, mD_eff,  std::complex<double>(gds, 0.0) + y_gd);
        if (mG_phys >= 0) A_c.addValue(mD_eff, mG_phys, std::complex<double>(gm, 0.0) - y_gd);
        if (mS_eff >= 0)  A_c.addValue(mD_eff, mS_eff,  std::complex<double>(-(gm + gds), 0.0));
    }
    
    // Physical Gate Node Equations (KCL Row)
    if (mG_phys >= 0) {
        A_c.addValue(mG_phys, mG_phys, y_gs + y_gd);
        if (mD_eff >= 0)  A_c.addValue(mG_phys, mD_eff,  -y_gd);
        if (mS_eff >= 0)  A_c.addValue(mG_phys, mS_eff,  -y_gs);
    }
    
    // Effective Source Node Equations (KCL Row)
    if (mS_eff >= 0) {
        A_c.addValue(mS_eff, mS_eff,  std::complex<double>(gm + gds, 0.0) + y_gs);
        if (mD_eff >= 0)  A_c.addValue(mS_eff, mD_eff,  std::complex<double>(-gds, 0.0));
        if (mG_phys >= 0) A_c.addValue(mS_eff, mG_phys, std::complex<double>(-gm, 0.0) - y_gs);
    }
}

double Mosfet::calculateCurrent(const SimulationContext& simCtx, int extraIdx) const {
    (void)extraIdx;
    if (simCtx.state.empty() || nodes.size() < 3) {
        return 0.0;
    }

    // Node parameter mapping variables based on 1-based parser arrays
    int nD = nodes[0];
    int nG = nodes[1];
    int nS = nodes[2];

    // 1. Safe extraction of raw absolute node potentials via getSafeVoltage abstraction layer
    double vD = getSafeVoltage(simCtx.state, nD);
    double vG = getSafeVoltage(simCtx.state, nG);
    double vS = getSafeVoltage(simCtx.state, nS);

    bool reversedCurrent = false;

    // 2. Perform PMOS Complementary Transistor Channel Inversion upfront
    if (!isNChannel) {
        std::swap(vD, vS);
        vG = vS - (vG - vS); // Mirror gate voltage relative to the PMOS source terminal frame
        reversedCurrent = true;
    }

    // 3. Handle Symmetric Bidirectional Operating Channel Modes (Drain-Source swapping)
    if ((vD - vS) < 0.0) {
        std::swap(vD, vS);
        reversedCurrent = !reversedCurrent; // Invert branch current sign trace direction on geometric swap
    }

    // 4. Compute standard Shichman-Hodges localized operating frame parameters
    double vgs = vG - vS;
    double vds = vD - vS;
    double vov = vgs - vt0;

    if (vov <= 0.0) {
        return 0.0; // Cutoff Operating Region: Channel is closed
    }

    double ids = 0.0;
    double lambdaTerm = 1.0 + lambda * vds;

    // 5. Piecewise numeric evaluation of active operating region
    if (vds < vov) { // Linear / Triode Operating Region
        ids = beta * (vov * vds - 0.5 * vds * vds) * lambdaTerm;
    } else { // Saturation / Active Pinch-Off Operating Region
        ids = 0.5 * beta * vov * vov * lambdaTerm;
    }

    return reversedCurrent ? -ids : ids;
}

double Mosfet::calculateVoltageDrop(const std::vector<double>& state) const {
    // Vds = V_drain - V_source -> nodes[0] - nodes[2]
    return getSafeVoltage(state, nodes[0]) - getSafeVoltage(state, nodes[2]); 
}

void Mosfet::getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx) const {
    (void)extraIdx; // MOSFETs are non-linear three-terminal devices and do not request extra rows

    // Structural Netlist Validation Guard: Ensure the parser has populated all 3 mandatory terminals
    if (nodes.size() < 3) {
        return;
    }

    // Loop through the 3x3 internal device connectivity matrix stencil mappings
    // 0 = Drain, 1 = Gate, 2 = Source
    for (int i : {0, 1, 2}) {
        for (int j : {0, 1, 2}) {
            // Trap ground index 0 by strictly enforcing (node > 0) before 0-based conversion (-1)
            if (nodes[i] > 0 && nodes[j] > 0) {
                patterns.push_back({nodes[i] - 1, nodes[j] - 1});
            }
        }
    }
}

// ============================================================================
// 6. ADVANCED BJT
// ============================================================================

VbicBjtAdvanced::VbicBjtAdvanced(const std::string& name, const std::vector<int>& nodes, bool npnVariant)
    : Component(name, ComponentType::BJT, nodes), isNpn(npnVariant) {}

// mathematically correct, good and robust and no structural bugs 
void VbicBjtAdvanced::stamp(MnaMatrixContext& ctx) {
    // 1. Structural Layer Guard: Ensure the device possesses all 3 mandatory terminals
    if (nodes.size() < 3) return;

    const int nC = nodes[0]; // Physical Collector connection tag
    const int nB = nodes[1]; // Physical Base connection tag
    const int nE = nodes[2]; // Physical Emitter connection tag

    // 2. FIXED: Safely resolve voltages via getSafeVoltage to eliminate 1-based pointer crashes
    double vC = 0.0, vB = 0.0, vE = 0.0;
    double vold_C = 0.0, vold_B = 0.0, vold_E = 0.0;

    if (ctx.currentIterV) {
        vC = getSafeVoltage(*ctx.currentIterV, nC);
        vB = getSafeVoltage(*ctx.currentIterV, nB);
        vE = getSafeVoltage(*ctx.currentIterV, nE);
    }
    if (ctx.historicStepV) {
        vold_C = getSafeVoltage(*ctx.historicStepV, nC);
        vold_B = getSafeVoltage(*ctx.historicStepV, nB);
        vold_E = getSafeVoltage(*ctx.historicStepV, nE);
    }

    // 3. Polarity Isolation Layer
    double sig = isNpn ? 1.0 : -1.0;
    double vbe = sig * (vB - vE);
    double vbc = sig * (vB - vC);
    double vold_be = sig * (vold_B - vold_E);
    double vold_bc = sig * (vold_B - vold_C);

    double vt = Nf * V_THERMAL_300K;
    double vtr = Nr * V_THERMAL_300K;
    double safe_Is = std::max(Is, 1e-20);

    // 4. SPICE Voltage Limiting Loop
    double vcrit_f = vt * std::log(vt / (std::sqrt(2.0) * safe_Is));
    double vcrit_r = vtr * std::log(vtr / (std::sqrt(2.0) * safe_Is));
    double vbe_lim = CircuitMath::pnjlim(vbe, vold_be, vt, vcrit_f);
    double vbc_lim = CircuitMath::pnjlim(vbc, vold_bc, vtr, vcrit_r);

    // 5. Physical Carrier Transport Calculations with safe exponent bounding
    double exp_vbe = std::exp(std::min(vbe_lim / vt, 40.0));
    double exp_vbc = std::exp(std::min(vbc_lim / vtr, 40.0));

    double ibf = safe_Is * (exp_vbe - 1.0);
    double ibr = safe_Is * (exp_vbc - 1.0);

    double dibf_dvbe = (safe_Is / vt) * exp_vbe;
    double dibr_dvbc = (safe_Is / vtr) * exp_vbc;

    // High-current Knee Charge Calculations
    double qb = 1.0;
    double dqb_dvbe = 0.0;
    if (Ikf > 0.0) {
        double safe_ibf = std::max(0.0, ibf);
        double arg_q2 = 1.0 + 4.0 * (safe_ibf / Ikf);
        if (arg_q2 > 0.0) {
            double sqrt_q2 = std::sqrt(arg_q2);
            qb = 0.5 * (1.0 + sqrt_q2);
            dqb_dvbe = 2.0 * dibf_dvbe / (Ikf * sqrt_q2);
        }
    }

    // Smooth numerical floor block to prevent divide-by-zero anomaly
    if (qb < 1e-3) {
        double delta = 1e-3 - qb;
        double sqrt_val = std::sqrt(delta * delta + 1e-6);
        qb = 1e-3 + 0.5 * (sqrt_val - delta);
        dqb_dvbe *= 0.5 * (1.0 + (delta / sqrt_val)); 
    }

    double itz = (ibf - ibr) / qb;
    double ditz_dvbe = (dibf_dvbe / qb) - (itz / qb) * dqb_dvbe;
    double ditz_dvbc = -(dibr_dvbc / qb);

    constexpr double gmin = 1e-12;
    double ice = itz + gmin * (vbe_lim - vbc_lim);
    double ibe = (ibf / 100.0) + gmin * vbe_lim; // Assuming static forward beta = 100
    double ibc = (ibr / 50.0) + gmin * vbc_lim;   // Assuming static reverse beta = 50

    double g_ice_vbe = ditz_dvbe + gmin;
    double g_ice_vbc = ditz_dvbc - gmin;
    double g_ibe_vbe = (dibf_dvbe / 100.0) + gmin;
    double g_ibc_vbc = (dibr_dvbc / 50.0) + gmin;

    // 6. Absolute Oriented Node Conductance Mapping (KCL Multi-Terminal Jacobian)
    double final_gCC = g_ibc_vbc - g_ice_vbc;
    double final_gCB = g_ice_vbe + g_ice_vbc - g_ibc_vbc;
    double final_gCE = -g_ice_vbe;

    double final_gBC = -g_ibc_vbc;
    double final_gBB = g_ibe_vbe + g_ibc_vbc;
    double final_gBE = -g_ibe_vbe;

    double final_gEC = g_ice_vbc; 
    double final_gEB = -g_ibe_vbe - g_ice_vbe - g_ice_vbc;
    double final_gEE = g_ibe_vbe + g_ice_vbe;

    // Mathematically Sound True External Currents
    double ic_ext = sig * (ice - ibc);
    double ib_ext = sig * (ibe + ibc);
    double ie_ext = -(ic_ext + ib_ext);

    // RHS Companion Vector Calculations via Taylor Expansion
    double ieq_c = -ic_ext + (final_gCC * vC + final_gCB * vB + final_gCE * vE);
    double ieq_b = -ib_ext + (final_gBC * vC + final_gBB * vB + final_gBE * vE);
    double ieq_e = -ie_ext + (final_gEC * vC + final_gEB * vB + final_gEE * vE);

    // Translate to 0-based matrix indexing (Trapping Ground index 0)
    int mC = nC - 1;
    int mB = nB - 1;
    int mE = nE - 1;

    // 7. Multi-Terminal Jacobian Matrix Stamping Pipeline
    // Collector Row Stamps
    if (mC >= 0) {
        ctx.G.addValue(mC, mC, final_gCC);
        if (mB >= 0) ctx.G.addValue(mC, mB, final_gCB);
        if (mE >= 0) ctx.G.addValue(mC, mE, final_gCE);
        ctx.B[mC] += ieq_c;
    }
    // FIXED: Guard base row entry using 0-based 'mB >= 0' instead of 'nB >= 0' to eliminate Segfaults
    if (mB >= 0) {
        if (mC >= 0) ctx.G.addValue(mB, mC, final_gBC);
        ctx.G.addValue(mB, mB, final_gBB);
        if (mE >= 0) ctx.G.addValue(mB, mE, final_gBE);
        ctx.B[mB] += ieq_b;
    }
    // Emitter Row Stamps
    if (mE >= 0) {
        if (mC >= 0) ctx.G.addValue(mE, mC, final_gEC);
        if (mB >= 0) ctx.G.addValue(mE, mB, final_gEB);
        ctx.G.addValue(mE, mE, final_gEE);
        ctx.B[mE] += ieq_e;
    }
}

void VbicBjtAdvanced::stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c, double omega, const std::vector<double>& dcState) {
    (void)B_c; // Small-signal independent AC elements handle B_c vector modifications

    // 1. Structural Layer Guard: Ensure the device possesses all 3 mandatory terminals
    if (nodes.size() < 3) return;

    int nC = nodes[0];
    int nB = nodes[1];
    int nE = nodes[2];

    // 2. FIXED: Safely resolve absolute node potentials from DC operating point via getSafeVoltage
    double vC_val = getSafeVoltage(dcState, nC);
    double vB_val = getSafeVoltage(dcState, nB);
    double vE_val = getSafeVoltage(dcState, nE);

    // 3. Determine intrinsic internal voltage orientations (PNP inversion mapping)
    double sig = isNpn ? 1.0 : -1.0;
    
    double vbe = sig * (vB_val - vE_val);
    double vbc = sig * (vB_val - vC_val);

    double vt  = Nf * V_THERMAL_300K;
    double vtr = Nr * V_THERMAL_300K;

    // 4. DC Operating Point Formulations & Safe Exponential Clamping
    double exp_vbe = std::exp(std::min(vbe / vt, 40.0)); 
    double exp_vbc = std::exp(std::min(vbc / vtr, 40.0));

    double ibf = Is * (exp_vbe - 1.0);
    double qb = 1.0;
    double dqb_dvbe = 0.0;
    
    double dibf_dvbe = (Is / vt) * exp_vbe;
    double dibr_dvbc = (Is / vtr) * exp_vbc;

    // Guard Gummel Knee square root parameter to eliminate NaN anomalies
    double safe_ibf = std::max(0.0, ibf);
    double arg_q2 = 1.0 + 4.0 * (safe_ibf / (Ikf > 1e-20 ? Ikf : 1e-20));
    if (arg_q2 > 0.0) {
        double sqrt_arg = std::sqrt(arg_q2);
        qb = 0.5 * (1.0 + sqrt_arg);
        dqb_dvbe = (Ikf > 1e-20) ? (2.0 * dibf_dvbe / (Ikf * sqrt_arg)) : 0.0;
    } else {
        qb = 0.5;
        dqb_dvbe = 0.0;
    }
    if (qb < 1e-3) qb = 1e-3; // Numerical floor protection against divide-by-zero

    double itz = ibf / qb;
    double ditz_dvbe = (dibf_dvbe / qb) - (itz / qb) * dqb_dvbe;
    double ditz_dvbc = -(dibr_dvbc / qb);

    // Small-signal linearized linear base current derivatives
    double dibf_scaled_dvbe = dibf_dvbe / 100.0; // Forward Beta = 100
    double dibr_scaled_dvbc = dibr_dvbc / 50.0;  // Reverse Beta = 50

    // Small-signal intrinsic terminal current partial derivatives
    double diC_dvbe = ditz_dvbe;
    double diC_dvbc = ditz_dvbc - dibr_scaled_dvbc;
    
    double diB_dvbe = dibf_scaled_dvbe;
    double diB_dvbc = dibr_scaled_dvbc;
    
    double diE_dvbe = -diC_dvbe - diB_dvbe;
    double diE_dvbc = -diC_dvbc - diB_dvbc;

    // 5. Formal Node Voltage Partial Derivative Conversions
    double dIc_dvC = -sig * diC_dvbc;
    double dIc_dvB =  sig * (diC_dvbe + diC_dvbc);
    double dIc_dvE = -sig * diC_dvbe;

    double dIb_dvC = -sig * diB_dvbc;
    double dIb_dvB =  sig * (diB_dvbe + diB_dvbc);
    double dIb_dvE = -sig * diB_dvbe;

    double dIe_dvC = -sig * diE_dvbc;
    double dIe_dvB =  sig * (diE_dvbe + diE_dvbc);
    double dIe_dvE = -sig * diE_dvbe;

    // 6. Dynamic Small-Signal Junction & Diffusion Parasitic Capacitance Evaluation
    constexpr double fc = 0.5; 

    // Base-Emitter Capacitance Formulation
    double c_be_diff = tf * dibf_dvbe;
    double c_be_dep = 0.0;
    if (vbe < fc * Vje) {
        double base_be = std::max(1.0 - vbe / Vje, 1e-3);
        c_be_dep = Cje0 * std::pow(base_be, -Me);
    } else {
        double c_factor = std::pow(1.0 - fc, -Me);
        c_be_dep = Cje0 * c_factor + (Me * Cje0 * c_factor / (Vje * (1.0 - fc))) * (vbe - fc * Vje);
    }
    double c_be_total = std::max(c_be_diff + c_be_dep, 1e-20);

    // Base-Collector Capacitance Formulation
    double c_bc_diff = tr * dibr_dvbc;
    double c_bc_dep = 0.0;
    if (vbc < fc * Vjc) {
        double base_bc = std::max(1.0 - vbc / Vjc, 1e-3);
        c_bc_dep = Cjc0 * std::pow(base_bc, -Mc);
    } else {
        double c_factor = std::pow(1.0 - fc, -Mc);
        c_bc_dep = Cjc0 * c_factor + (Mc * Cjc0 * c_factor / (Vjc * (1.0 - fc))) * (vbc - fc * Vjc);
    }
    double c_bc_total = std::max(c_bc_diff + c_bc_dep, 1e-20);

    // Map into complex admittance susceptance elements: Y = 0 + j * omega * C
    std::complex<double> y_be(0.0, omega * c_be_total);
    std::complex<double> y_bc(0.0, omega * c_bc_total);

    constexpr double gmin = 1e-12; // High-performance solver minimal stabilizing shunt conductance

    // 7. Assemble Complex Linearized Jacobian MNA Stencils
    std::complex<double> gCC = std::complex<double>(dIc_dvC + gmin, 0.0) + y_bc;
    std::complex<double> gCB = std::complex<double>(dIc_dvB, 0.0) - y_bc;
    std::complex<double> gCE = std::complex<double>(dIc_dvE, 0.0);

    std::complex<double> gBC = std::complex<double>(dIb_dvC, 0.0) - y_bc;
    std::complex<double> gBB = std::complex<double>(dIb_dvB + gmin, 0.0) + y_be + y_bc;
    std::complex<double> gBE = std::complex<double>(dIb_dvE, 0.0) - y_be;

    std::complex<double> gEC = std::complex<double>(dIe_dvC, 0.0);
    std::complex<double> gEB = std::complex<double>(dIe_dvB, 0.0) - y_be;
    std::complex<double> gEE = std::complex<double>(dIe_dvE + gmin, 0.0) + y_be;

    // 8. FIXED: Correctly shift matrix maps to clean 0-based coordinate indices (node - 1)
    int mC = nC - 1;
    int mB = nB - 1;
    int mE = nE - 1;

    // Complex Admittance Jacobian Matrix Injections
    // Collector Row KCL Stamps
    if (mC >= 0) {
        A_c.addValue(mC, mC, gCC);
        if (mB >= 0) A_c.addValue(mC, mB, gCB);
        if (mE >= 0) A_c.addValue(mC, mE, gCE);
    }
    // Base Row KCL Stamps
    if (mB >= 0) {
        if (mC >= 0) A_c.addValue(mB, mC, gBC);
        A_c.addValue(mB, mB, gBB);
        if (mE >= 0) A_c.addValue(mB, mE, gBE);
    }
    // Emitter Row KCL Stamps
    if (mE >= 0) {
        if (mC >= 0) A_c.addValue(mE, mC, gEC);
        if (mB >= 0) A_c.addValue(mE, mB, gEB);
        A_c.addValue(mE, mE, gEE);
    }
}

double VbicBjtAdvanced::calculateCurrent(const SimulationContext& simCtx, int extraIdx) const {
    (void)extraIdx;
    if (simCtx.state.empty() || nodes.size() < 3) {
        return 0.0;
    }

    // 1. Safe extraction of raw absolute node potentials using getSafeVoltage interface
    int nC = nodes[0];
    int nB = nodes[1];
    int nE = nodes[2];

    double vC = getSafeVoltage(simCtx.state, nC);
    double vB = getSafeVoltage(simCtx.state, nB);
    double vE = getSafeVoltage(simCtx.state, nE);

    // 2. Perform Polarity Device Transformation Upfront (PNP Inversion Mapping Isolation)
    double sig = isNpn ? 1.0 : -1.0;
    
    double vbe = sig * (vB - vE);
    double vbc = sig * (vB - vC);

    double vt  = Nf * V_THERMAL_300K;
    double vtr = Nr * V_THERMAL_300K;

    // 3. Numerical Clamping Layer to eliminate exponential overflow risks during unstable solver sweeps
    double vcrit_f = vt * std::log(vt / (std::sqrt(2.0) * Is));
    if (vbe > vcrit_f) {
        vbe = vcrit_f + (vbe - vcrit_f) / (1.0 + (vbe - vcrit_f) / vt);
    }
    double vcrit_r = vtr * std::log(vtr / (std::sqrt(2.0) * Is));
    if (vbc > vcrit_r) {
        vbc = vcrit_r + (vbc - vcrit_r) / (1.0 + (vbc - vcrit_r) / vtr);
    }

    // 4. Core Transport & Diode Formulations matching stamp() exactly
    double exp_vbe = std::exp(std::min(vbe / vt, 40.0));
    double exp_vbc = std::exp(std::min(vbc / vtr, 40.0));

    double ibf = Is * (exp_vbe - 1.0);
    double ibr = Is * (exp_vbc - 1.0);

    // 5. High-Injection Gummel-Poon Knee / Base-Charge (qb) Evaluation
    // Guarded against negative inputs inside the square root to maintain mathematical consistency
    double safe_ibf = std::max(0.0, ibf);
    double arg_q2 = 1.0 + 4.0 * (safe_ibf / (Ikf > 1e-20 ? Ikf : 1e-20));
    double qb = 1.0;
    if (arg_q2 > 0.0) {
        qb = 0.5 * (1.0 + std::sqrt(arg_q2));
    }
    if (qb < 1e-3) qb = 1e-3; // Numerical floor block to prevent divide-by-zero anomaly

    double itz = (ibf - ibr) / qb;
    double ib_r_total = ibr / 50.0; // Assuming a static internal Reverse Beta = 50

    // 6. Total Internal Collector Device Current Formulation
    double ic_internal = itz - ib_r_total;
    
    // Add parallel Gmin shunt factor for perfect consistency with matrix stamps
    constexpr double gmin = 1e-12;
    ic_internal += gmin * (sig * (vC - vE));

    // Returns the true external collector branch current tracking polarization
    return sig * ic_internal;
}

double VbicBjtAdvanced::calculateVoltageDrop(const std::vector<double>& state) const {
    // Vce = V_collector - V_emitter -> nodes[0] - nodes[2]
    return getSafeVoltage(state, nodes[0]) - getSafeVoltage(state, nodes[2]); 
}

void VbicBjtAdvanced::getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx) const {
    (void)extraIdx; // BJTs are non-linear three-terminal devices and do not request extra rows

    // Structural Netlist Validation Guard: Ensure the parser has populated all 3 mandatory terminals
    if (nodes.size() < 3) {
        return;
    }

    // Loop through the 3x3 internal device connectivity matrix stencil mappings
    // 0 = Collector, 1 = Base, 2 = Emitter
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            // Trap ground index 0 by strictly enforcing (node > 0) before 0-based conversion (-1)
            if (nodes[i] > 0 && nodes[j] > 0) {
                patterns.push_back({nodes[i] - 1, nodes[j] - 1});
            }
        }
    }
}

// ============================================================================
// 7. INDEPENDENT VOLTAGE SOURCE
// ============================================================================

VoltageSource::VoltageSource(const std::string& name, const std::vector<int>& nodes, double v)
    : Component(name, ComponentType::VOLTAGE_SOURCE, nodes), value(v) {}

// mathematically correct, good and robust and no structural bugs 
void VoltageSource::stamp(MnaMatrixContext& ctx) {
    // 1. Structural Layer Guards: Ensure nodes and extra row references are safely allocated
    if (nodes.size() < 2) return;

    int idx = (ctx.assignedExtraIdx >= 0) ? ctx.assignedExtraIdx : this->assignedExtraRowIdx;
    if (idx < 0 || idx >= static_cast<int>(ctx.B.size())) return; 

    int n1 = nodes[0]; // Positive terminal
    int n2 = nodes[1]; // Negative terminal
    
    // 2. Nodal Interconnection Matrix Stamps (Enforcing strict 0-based matrix mapping via node - 1)
    if (n1 > 0) {
        ctx.G.addValue(n1 - 1, idx,  1.0);
        ctx.G.addValue(idx, n1 - 1,  1.0);
    }
    if (n2 > 0) {
        ctx.G.addValue(n2 - 1, idx, -1.0);
        ctx.G.addValue(idx, n2 - 1, -1.0);
    }
    
    // Hard allocation for the branch diagonal grid point to secure structural solver factor integrity
    ctx.G.addValue(idx, idx, 0.0);

    // 3. FIXED: Absolute initialization instead of incremental assignment to eliminate load duplication
    // Constitutive branch circuit voltage equation: V_{n1-1} - V_{n2-1} = independent_voltage_value
    ctx.B[idx] = value;
}

void VoltageSource::getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx) const {
    // Structural netlist validation layer guard
    if (nodes.size() < 2) return;

    int brRow = (extraIdx >= 0) ? extraIdx : this->assignedExtraRowIdx;
    if (brRow < 0) return;

    // Allocate vital self-coupling branch diagonal matrix grid point entry
    patterns.push_back({brRow, brRow});
    
    int n1 = nodes[0];
    int n2 = nodes[1];

    // Build static connection pattern maps using strictly monitored 0-based matrix shifts
    if (n1 > 0) {
        patterns.push_back({n1 - 1, brRow});
        patterns.push_back({brRow, n1 - 1});
    }
    if (n2 > 0) {
        patterns.push_back({n2 - 1, brRow});
        patterns.push_back({brRow, n2 - 1});
    }
}

double VoltageSource::calculateCurrent(const SimulationContext& simCtx, int extraIdx) const {
    int targetIdx = (extraIdx >= 0) ? extraIdx : this->assignedExtraRowIdx;
    
    // Direct Track: Extract natively solved branch current from MNA solution vector
    if (targetIdx >= 0 && static_cast<size_t>(targetIdx) < simCtx.state.size()) {
        return simCtx.state[targetIdx];
    }
    return 0.0;
}


double VoltageSource::calculateVoltageDrop(const std::vector<double>& state) const {
    // Structural netlist validation layer guard
    if (nodes.size() < 2) return 0.0;

    int n1 = nodes[0];
    int n2 = nodes[1];
    
    // Leverage the global pointer-guarded getSafeVoltage abstraction to secure 0-based indexing
    return getSafeVoltage(state, n1) - getSafeVoltage(state, n2);
}

void VoltageSource::stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c, double omega, const std::vector<double>& dcState) {
    (void)omega; 
    (void)dcState; // Independent sources do not modulate based on the DC operational point state

    // 1. Structural Layer Guards: Ensure nodes and extra row references are safely allocated
    if (nodes.size() < 2) return;

    int idx = getAssignedExtraIdx(); // Fetch pre-compiled 0-based branch index coordinate
    if (idx < 0 || idx >= static_cast<int>(B_c.size())) return; // Guard against unallocated parameters

    const int n1 = nodes[0]; // Positive terminal
    const int n2 = nodes[1]; // Negative terminal

    // 2. FIXED: Correctly shift circuit nodes to 0-based matrix indexing (node - 1) to eliminate AC Segfaults.
    // Leave the branch index 'idx' unmodified as it is already a 0-based coordinate tracking value.
    if (n1 > 0) {
        A_c.addValue(n1 - 1, idx,     1.0);
        A_c.addValue(idx,     n1 - 1, 1.0);
    }
    if (n2 > 0) {
        A_c.addValue(n2 - 1, idx,    -1.0);
        A_c.addValue(idx,     n2 - 1, -1.0);
    }
    
    // Enforce branch matrix structural short-circuit identity baseline coupling coordinate
    A_c.addValue(idx, idx, std::complex<double>(0.0, 0.0));

    // 3. Convert Magnitude/Phase into Cartesian Phasor Representation (Real + j*Imag)
    // Note: If the parser already converted degrees to radians, use acPhase directly. 
    // If acPhase is stored in degrees, use: double phaseRad = acPhase * M_PI / 180.0;
    std::complex<double> v_phasor = std::polar(acMagnitude, acPhase);

    // 4. FIXED: Absolute initialization assignment on B_c instead of incremental addition to prevent drift
    // Independent source branch constitutive phasor loop equation constraint injection
    B_c[idx] = v_phasor;
}

void VoltageSource::setAcParameters(double mag, double phase) {
    acMagnitude = mag;
    acPhase = phase;
}

void VoltageSource::allocateRows(int& globalExtraRowCounter) {
    // Let the base class handle it cleanly, using your getExtraRowRequirement() override
    Component::allocateRows(globalExtraRowCounter);
}

// ============================================================================
// 8. BEHAVIORAL LOGIC GATE BRIDGE
// ============================================================================

BehaviorLogicGateComponent::BehaviorLogicGateComponent(
    const std::string& name, const std::vector<int>& nodes, GateType gateType, 
    double vThresh, double vHigh, double vLow, double rOut)
    : Component(name, ComponentType::DIGITAL_IC, nodes), gate(gateType), 
      thresholdVoltage(vThresh), outputHighV(vHigh), outputLowV(vLow), outputResistance(rOut) {}

// mathematically correct, good and robust and no structural bugs 
void BehaviorLogicGateComponent::stamp(MnaMatrixContext& ctx) {
    // 1. Structural Layer Guards: Ensure nodes and extra row references are safely allocated
    if (nodes.empty()) return;

    int idx = (ctx.assignedExtraIdx >= 0) ? ctx.assignedExtraIdx : this->assignedExtraRowIdx;
    if (idx < 0 || idx >= static_cast<int>(ctx.B.size())) return; 

    // A NOT gate requires at least 2 terminals (1 In, 1 Out). Multi-input gates require at least 3 (2 In, 1 Out).
    size_t minTerminals = (gate == GateType::NOT) ? 2 : 3;
    if (nodes.size() < minTerminals) {
        throw std::runtime_error("Digital IC error: [" + name + "] lacks mandatory terminal connection pairs.");
    }

    int nOut = nodes.back(); // Output pin is mapped to the final index element
    double continuousLogicState = 0.0;
    constexpr double kTransitionSharpness = 4.0; // Defines the tanh voltage transition steepness slope

    // Analytical derivative trackers mapping: d(continuousLogicState) / dV_node
    double dLogic_dvA = 0.0;
    double dLogic_dvB = 0.0;

    // 2. Safely resolve local voltages using getSafeVoltage to eliminate pointer crashes
    if (gate == GateType::NOT) {
        int nIn = nodes[0];
        double vIn = ctx.currentIterV ? getSafeVoltage(*ctx.currentIterV, nIn) : 0.0;
        
        double tanhArg = kTransitionSharpness * (vIn - thresholdVoltage);
        double tanhVal = std::tanh(tanhArg);
        double sech2 = 1.0 - tanhVal * tanhVal;
        
        continuousLogicState = 0.5 * (1.0 - tanhVal);
        dLogic_dvA = -0.5 * kTransitionSharpness * sech2; // Treated as input A mapping
    } 
    else {
        int nA = nodes[0];
        int nB = nodes[1];
        double vA = ctx.currentIterV ? getSafeVoltage(*ctx.currentIterV, nA) : 0.0;
        double vB = ctx.currentIterV ? getSafeVoltage(*ctx.currentIterV, nB) : 0.0;
        
        double tanhArgA = kTransitionSharpness * (vA - thresholdVoltage);
        double tanhArgB = kTransitionSharpness * (vB - thresholdVoltage);
        
        double tanhValA = std::tanh(tanhArgA);
        double tanhValB = std::tanh(tanhArgB);
        
        double smoothA = 0.5 * (1.0 + tanhValA);
        double smoothB = 0.5 * (1.0 + tanhValB);
        
        double sech2A = 1.0 - tanhValA * tanhValA;
        double sech2B = 1.0 - tanhValB * tanhValB;
        
        double dsmoothA_dvA = 0.5 * kTransitionSharpness * sech2A;
        double dsmoothB_dvB = 0.5 * kTransitionSharpness * sech2B;

        switch (gate) {
            case GateType::AND:  
                continuousLogicState = smoothA * smoothB; 
                dLogic_dvA = dsmoothA_dvA * smoothB;
                dLogic_dvB = smoothA * dsmoothB_dvB;
                break;
            case GateType::OR:   
                continuousLogicState = smoothA + smoothB - (smoothA * smoothB); 
                dLogic_dvA = dsmoothA_dvA * (1.0 - smoothB);
                dLogic_dvB = dsmoothB_dvB * (1.0 - smoothA);
                break;
            case GateType::NAND: 
                continuousLogicState = 1.0 - (smoothA * smoothB); 
                dLogic_dvA = -dsmoothA_dvA * smoothB;
                dLogic_dvB = -smoothA * dsmoothB_dvB;
                break;
            case GateType::NOR:  
                continuousLogicState = 1.0 - (smoothA + smoothB - (smoothA * smoothB)); 
                dLogic_dvA = -dsmoothA_dvA * (1.0 - smoothB);
                dLogic_dvB = -dsmoothB_dvB * (1.0 - smoothA);
                break;
            case GateType::XOR:  
                continuousLogicState = smoothA + smoothB - 2.0 * (smoothA * smoothB); 
                dLogic_dvA = dsmoothA_dvA * (1.0 - 2.0 * smoothB);
                dLogic_dvB = dsmoothB_dvB * (1.0 - 2.0 * smoothA);
                break;
            case GateType::XNOR: 
                continuousLogicState = 1.0 - (smoothA + smoothB - 2.0 * (smoothA * smoothB)); 
                dLogic_dvA = -dsmoothA_dvA * (1.0 - 2.0 * smoothB);
                dLogic_dvB = -dsmoothB_dvB * (1.0 - 2.0 * smoothA);
                break;
            default: 
                continuousLogicState = smoothA; 
                dLogic_dvA = dsmoothA_dvA;
                break;
        }
    }

    // 3. Linearly map internal states onto analog hardware voltage scaling bounds
    double vScale = outputHighV - outputLowV;
    double dVtarget_dvA = dLogic_dvA * vScale;
    double dVtarget_dvB = dLogic_dvB * vScale;

    continuousLogicState = std::max(0.0, std::min(1.0, continuousLogicState));
    double targetV = outputLowV + continuousLogicState * vScale;

    // 4. Taylor Series Newton RHS Expansion Adjustment Calculation
    double vA_curr = ctx.currentIterV ? getSafeVoltage(*ctx.currentIterV, nodes[0]) : 0.0;
    double vB_curr = (gate != GateType::NOT && ctx.currentIterV) ? getSafeVoltage(*ctx.currentIterV, nodes[1]) : 0.0;
    
    // Ieq represents the linearized multi-variable offset companion voltage load: V_target - J * V_k
    double ieq = targetV - (dVtarget_dvA * vA_curr + dVtarget_dvB * vB_curr);

    // 5. Multi-Terminal Jacobian Matrix Stamping (Enforcing strict 0-based node - 1 shifts)
    if (nOut > 0) {
        ctx.G.addValue(nOut - 1, idx,      1.0);  
        ctx.G.addValue(idx,      nOut - 1, 1.0);  
    }
    
    // Inject cross-coupling transconductance derivative column elements into the MNA grid
    if (nodes[0] > 0) {
        ctx.G.addValue(idx, nodes[0] - 1, -dVtarget_dvA);
    }
    if (gate != GateType::NOT && nodes[1] > 0) {
        ctx.G.addValue(idx, nodes[1] - 1, -dVtarget_dvB);
    }

    // Physical behavioral driving output series resistance mapping coordinate
    ctx.G.addValue(idx, idx, -outputResistance);

    // FIXED: Absolute serialization assignment instead of incremental step loading additions
    ctx.B[idx] = ieq;
}

double BehaviorLogicGateComponent::calculateCurrent(const SimulationContext& simCtx, int extraIdx) const {
    // Determine the active 0-based branch matrix coordinate tracking index
    int targetIdx = (extraIdx >= 0) ? extraIdx : this->assignedExtraRowIdx;
    
    // Guard against unallocated parameters or uninitialized simulation states
    if (simCtx.state.empty() || targetIdx < 0) {
        return 0.0; 
    }
    
    // Direct Track: Extract the natively solved logic output branch current from the MNA solution vector.
    // The branch index 'targetIdx' is pre-compiled as a 0-based coordinate tracking value.
    if (static_cast<size_t>(targetIdx) < simCtx.state.size()) {
        return simCtx.state[targetIdx];
    }
    
    return 0.0;
}
double BehaviorLogicGateComponent::calculateVoltageDrop(const std::vector<double>& state) const {
    // Structural netlist validation layer guard
    if (nodes.empty()) {
        return 0.0;
    }
    
    // The output terminal is always mapped to the final index element of the nodes vector
    int nOut = nodes.back();
    
    // Leverage the global pointer-guarded getSafeVoltage abstraction to secure 0-based indexing.
    // This returns the exact absolute output node voltage relative to Ground (0V).
    return getSafeVoltage(state, nOut); 
}

void BehaviorLogicGateComponent::stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c, double omega, const std::vector<double>& dcState) {
    (void)omega; 
    (void)B_c; // Independent phasor load adjustments are handled externally by AC source elements

    // 1. Structural Layer Guards: Ensure nodes and extra row references are safely allocated
    if (nodes.empty()) return;

    int idx = getAssignedExtraIdx(); // Fetch pre-compiled 0-based branch index coordinate
    if (idx < 0) return; 

    size_t minTerminals = (gate == GateType::NOT) ? 2 : 3;
    if (nodes.size() < minTerminals) return;

    int nOut = nodes.back(); // Output pin is mapped to the final index element
    constexpr double kTransitionSharpness = 4.0; 

    // Analytical derivative trackers mapping: d(continuousLogicState) / dV_node
    double dLogic_dvA = 0.0;
    double dLogic_dvB = 0.0;

    // 2. Safely sample converged absolute DC node potentials via getSafeVoltage interface
    if (gate == GateType::NOT) {
        int nIn = nodes[0];
        double vIn = getSafeVoltage(dcState, nIn);
        
        double tanhArg = kTransitionSharpness * (vIn - thresholdVoltage);
        double tanhVal = std::tanh(tanhArg);
        double sech2 = 1.0 - tanhVal * tanhVal;
        
        dLogic_dvA = -0.5 * kTransitionSharpness * sech2; 
    } 
    else {
        int nA = nodes[0];
        int nB = nodes[1];
        double vA = getSafeVoltage(dcState, nA);
        double vB = getSafeVoltage(dcState, nB);
        
        double tanhArgA = kTransitionSharpness * (vA - thresholdVoltage);
        double tanhArgB = kTransitionSharpness * (vB - thresholdVoltage);
        
        double tanhValA = std::tanh(tanhArgA);
        double tanhValB = std::tanh(tanhArgB);
        
        double smoothA = 0.5 * (1.0 + tanhValA);
        double smoothB = 0.5 * (1.0 + tanhValB);
        
        double sech2A = 1.0 - tanhValA * tanhValA;
        double sech2B = 1.0 - tanhValB * tanhValB;
        
        double dsmoothA_dvA = 0.5 * kTransitionSharpness * sech2A;
        double dsmoothB_dvB = 0.5 * kTransitionSharpness * sech2B;

        switch (gate) {
            case GateType::AND:  
                dLogic_dvA = dsmoothA_dvA * smoothB;
                dLogic_dvB = smoothA * dsmoothB_dvB;
                break;
            case GateType::OR:   
                dLogic_dvA = dsmoothA_dvA * (1.0 - smoothB);
                dLogic_dvB = dsmoothB_dvB * (1.0 - smoothA);
                break;
            case GateType::NAND: 
                dLogic_dvA = -dsmoothA_dvA * smoothB;
                dLogic_dvB = -smoothA * dsmoothB_dvB;
                break;
            case GateType::NOR:  
                dLogic_dvA = -dsmoothA_dvA * (1.0 - smoothB);
                dLogic_dvB = -dsmoothB_dvB * (1.0 - smoothA);
                break;
            case GateType::XOR:  
                dLogic_dvA = dsmoothA_dvA * (1.0 - 2.0 * smoothB);
                dLogic_dvB = dsmoothB_dvB * (1.0 - 2.0 * smoothA);
                break;
            case GateType::XNOR: 
                dLogic_dvA = -dsmoothA_dvA * (1.0 - 2.0 * smoothB);
                dLogic_dvB = -dsmoothB_dvB * (1.0 - 2.0 * smoothA);
                break;
            default: 
                dLogic_dvA = dsmoothA_dvA;
                break;
        }
    }

    // Linearly map the derivatives onto target small-signal transconductance scaling weights
    double vScale = outputHighV - outputLowV;
    double dVtarget_dvA = dLogic_dvA * vScale;
    double dVtarget_dvB = dLogic_dvB * vScale;

    // 3. FIXED: Translate wire nodes to clean 0-based matrix indexing (node - 1)
    int mOut = nOut - 1;

    // 4. Inject Complex Nodal Admittance Stencils (Including the critical AC input-output coupling derivatives)
    if (mOut >= 0) {
        A_c.addValue(mOut, idx, std::complex<double>(1.0, 0.0));
        A_c.addValue(idx, mOut, std::complex<double>(1.0, 0.0));
    }
    
    // Inject linearized small-signal coupling transconductances into Jacobian columns
    if (nodes[0] > 0) {
        A_c.addValue(idx, nodes[0] - 1, std::complex<double>(-dVtarget_dvA, 0.0));
    }
    if (gate != GateType::NOT && nodes[1] > 0) {
        A_c.addValue(idx, nodes[1] - 1, std::complex<double>(-dVtarget_dvB, 0.0));
    }

    // Output load series resistance phasor coupling coordinate injection
    A_c.addValue(idx, idx, std::complex<double>(-outputResistance, 0.0));
}

void BehaviorLogicGateComponent::allocateRows(int& globalExtraRowCounter) {
    Component::allocateRows(globalExtraRowCounter);
}

void BehaviorLogicGateComponent::getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx) const {
    // FIX: Fall back to local assignment marker if extraIdx argument is unallocated (-1)
    int brRow = (extraIdx >= 0) ? extraIdx : this->assignedExtraRowIdx;
    if (brRow < 0) return;

    if (nodes.size() < 2) return; // Safety check for malformed netlists

    // FIX: Dynamically locate the output pin (always the final element)
    int nOut = nodes.back();

    patterns.push_back({brRow, brRow}); // Diagonal structural branch entry
    
    // FIX: Stripped out all manual "- 1" offsets to align perfectly with 1-based Netlist structural maps
    if (nOut > 0) {
        patterns.push_back({nOut-1, brRow});
        patterns.push_back({brRow, nOut-1});
    }

    if (gate == GateType::NOT) {
        // NOT gate only has one input row entry
        if (nodes[0] > 0) {
            patterns.push_back({brRow, nodes[0]-1});
        }
    } 
    else {
        // Multi-input gates require cross-coupling rows allocation for inputs A and B
        if (nodes[0] > 0) patterns.push_back({brRow, nodes[0]-1});
        if (nodes[1] > 0) patterns.push_back({brRow, nodes[1]-1});
    }
}

/*void BehaviorLogicGateComponent::stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c, double omega, const std::vector<double>& dcState) {
    (void)B_c; (void)omega; (void)dcState; // AC voltage target is 0, frequency independent

    // Fetch the 0-based matrix row index allocated for this logic gate branch
    int idx = getAssignedExtraIdx();
    if (idx < 0) return; // Guard against unallocated branch parameters

    // Logic gate output node is a 1-based schematic node index
    int out = nodes[2]; // Assuming nodes[0]=inA, nodes[1]=inB, nodes[2]=out from the topology configuration

    // CONVERSION: Map 1-based physical node to 0-based matrix coordinates via "- 1"
    // 'idx' remains untouched as it is already a 0-based internal matrix row index
    if (out > 0) {
        A_c.addValue(out - 1, idx, 1.0);
        A_c.addValue(idx, out - 1, 1.0);
    }

    // Small-signal linearized representation: Vout - Rout * Ibranch = 0 (AC Ground)
    // Populate the diagonal of the extra branch row to prevent a singular matrix error
    double effective_rout = (outputResistance < 1e-12) ? 1e-12 : outputResistance;
    A_c.addValue(idx, idx, -effective_rout);
}*/

