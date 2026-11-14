
#include "freyalab/CircuitComponents.h"
#include "freyalab/math.h"
#include "freyalab/LogicGates.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <complex>
#include <iostream>


constexpr double V_THERMAL_300K = 0.0258563;

// ============================================================================
// 1. RESISTOR PRIMITIVE
// ============================================================================

Resistor::Resistor(const std::string& name, const std::vector<NetNode>& nodes, double r)
    : Component(name, ComponentType::RESISTOR, nodes), value(r) {}

void Resistor::stamp(MnaMatrixContext& ctx) {
    if (nodes.size() < 2) return;

    // 1. Division-by-zero safety clamping factor
    double safe_value = (std::abs(value) > 1e-20) ? value : 1e-20;
    double g = 1.0 / safe_value;

    // 2. Translate 1-based NetNodes into 0-based MnaIndices (Ground automatically becomes invalid)
    MnaIndex m1 = MnaIndex::fromNetNode(nodes[0]);
    MnaIndex m2 = MnaIndex::fromNetNode(nodes[1]);

    // 3. Inject conductance stencil with compile-time type safety
    if (m1.isValid()) ctx.G.addValue(m1.idx, m1.idx,  g);
    if (m2.isValid()) ctx.G.addValue(m2.idx, m2.idx,  g);

    if (m1.isValid() && m2.isValid()) {
        ctx.G.addValue(m1.idx, m2.idx, -g);
        ctx.G.addValue(m2.idx, m1.idx, -g);
    }
}

void Resistor::stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                       double omega, const std::vector<double>& dcState) {
    (void)B_c;
    (void)omega;
    (void)dcState;

    if (nodes.size() < 2) return;

    double safe_value = (std::abs(value) > 1e-20) ? value : 1e-20;
    double g = 1.0 / safe_value;
    std::complex<double> g_complex(g, 0.0);

    MnaIndex m1 = MnaIndex::fromNetNode(nodes[0]);
    MnaIndex m2 = MnaIndex::fromNetNode(nodes[1]);

    if (m1.isValid()) A_c.addValue(m1.idx, m1.idx,  g_complex);
    if (m2.isValid()) A_c.addValue(m2.idx, m2.idx,  g_complex);

    if (m1.isValid() && m2.isValid()) {
        A_c.addValue(m1.idx, m2.idx, -g_complex);
        A_c.addValue(m2.idx, m1.idx, -g_complex);
    }
}

void Resistor::getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns, MnaIndex extraIdx) const {
    (void)extraIdx;

    if (nodes.size() < 2) return;

    MnaIndex m1 = MnaIndex::fromNetNode(nodes[0]);
    MnaIndex m2 = MnaIndex::fromNetNode(nodes[1]);

    if (m1.isValid()) {
        patterns.push_back({m1, m1});
    }
    if (m2.isValid()) {
        patterns.push_back({m2, m2});
    }
    if (m1.isValid() && m2.isValid()) {
        patterns.push_back({m1, m2});
        patterns.push_back({m2, m1});
    }
}

double Resistor::calculateVoltageDrop(const std::vector<double>& state) const {
    if (nodes.size() < 2) return 0.0;
    return getSafeVoltage(state, nodes[0]) - getSafeVoltage(state, nodes[1]);
}

double Resistor::calculateCurrent(const SimulationContext& simCtx, MnaIndex extraIdx) const {
    (void)extraIdx;
    if (simCtx.state.empty()) {
        return 0.0;
    }

    // 1. Primary Path: Direct extraction via MNA assigned extra row index if pre-allocated
    MnaIndex idx = this->assignedExtraRowIdx;
    if (idx.isValid() && static_cast<size_t>(idx.idx) < simCtx.state.size()) {
        return simCtx.state[idx.idx];
    }

    // 2. Spatial Algebraic Path: Calculate instantaneously via Ohm's Law (I = V_drop / R)
    double safeValue = (std::abs(value) > 1e-20) ? value : 1e-20;
    double v_drop = calculateVoltageDrop(simCtx.state);
    return v_drop / safeValue;
}


// ============================================================================
// 2. DIODE MODEL
// ============================================================================
Diode::Diode(const std::string& name, const std::vector<NetNode>& nodes)
    : Component(name, ComponentType::DIODE, nodes) {}

void Diode::stamp(MnaMatrixContext& ctx) {
    if (nodes.size() < 2) return;

    NetNode n1 = nodes[0]; // Anode
    NetNode n2 = nodes[1]; // Cathode

    // 1. Resolve current iteration and historic voltage drops via bounds-checked NetNode lookups
    double vd_true = 0.0;
    double vold = 0.0;

    if (ctx.currentIterV) {
        vd_true = getSafeVoltage(ctx.currentIterV, n1) - getSafeVoltage(ctx.currentIterV, n2);
    }
    if (ctx.historicStepV) {
        vold = getSafeVoltage(ctx.historicStepV, n1) - getSafeVoltage(ctx.historicStepV, n2);
    }

    double safeIs = std::max(Is, 1e-20);
    double vt = N * V_THERMAL_300K;
    double vcrit = vt * std::log(vt / (std::sqrt(2.0) * safeIs));

    // 2. Compute limited voltage drop for exponential overflow safety
    double vd_lim = CircuitMath::pnjlim(vd_true, vold, vt, vcrit);

    // 3. Compute Shockley equations using the safely limited voltage drop
    double expTerm = std::exp(std::min(vd_lim / vt, 40.0));
    double id = safeIs * (expTerm - 1.0);
    double gd = (safeIs / vt) * expTerm;

    constexpr double gmin = 1e-12;
    gd += gmin;

    // 4. Formulate Norton equivalent current anchored at vd_lim
    double ieq = id - gd * vd_lim;

    // 5. Translate 1-based NetNodes to 0-based MnaIndices (trapping Ground automatically)
    MnaIndex m1 = MnaIndex::fromNetNode(n1);
    MnaIndex m2 = MnaIndex::fromNetNode(n2);

    // 6. Conductance Matrix G Stamping
    if (m1.isValid()) ctx.G.addValue(m1.idx, m1.idx,  gd);
    if (m2.isValid()) ctx.G.addValue(m2.idx, m2.idx,  gd);
    if (m1.isValid() && m2.isValid()) {
        ctx.G.addValue(m1.idx, m2.idx, -gd);
        ctx.G.addValue(m2.idx, m1.idx, -gd);
    }

    // 7. Bounds-Checked Source Vector B Stamping
    addSafeRhs(ctx.B, m1, -ieq);
    addSafeRhs(ctx.B, m2,  ieq);
}

void Diode::stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                    double omega, const std::vector<double>& dcState) {
    (void)B_c;

    if (nodes.size() < 2) return;

    NetNode n1 = nodes[0]; // Anode terminal
    NetNode n2 = nodes[1]; // Cathode terminal

    // 1. Extract linearized small-signal parameters from the converged DC Operating Point
    double v_dc = calculateVoltageDrop(dcState);
    double safeIs = std::max(Is, 1e-20);
    double vt = N * V_THERMAL_300K;

    double expTerm = std::exp(std::min(v_dc / vt, 40.0));
    double gd = (safeIs / vt) * expTerm;

    constexpr double gmin = 1e-12;
    gd += gmin;

    // 2. Compute Small-Signal Dynamic Junction Capacitance (C_total = C_diffusion + C_depletion)
    double c_diffusion = tt * gd;
    double c_depletion = 0.0;

    double v_break = fc * Vj;
    if (v_dc < v_break) {
        double base = std::max(1.0 - v_dc / Vj, 1e-3);
        c_depletion = Cj0 * std::pow(base, -M);
    } else {
        double c_factor = std::pow(1.0 - fc, -M);
        double slope = (Cj0 * M * std::pow(1.0 - fc, -(M + 1.0))) / Vj;
        c_depletion = Cj0 * c_factor + slope * (v_dc - v_break);
    }

    double c_total = std::max(c_diffusion + c_depletion, 1e-20);

    // 3. Transform Admittance into Complex Phasor Representation: Y = gd + j * omega * C_total
    std::complex<double> y_complex(gd, omega * c_total);

    // 4. Translate 1-based NetNodes to 0-based MnaIndices
    MnaIndex m1 = MnaIndex::fromNetNode(n1);
    MnaIndex m2 = MnaIndex::fromNetNode(n2);

    if (m1.isValid()) A_c.addValue(m1.idx, m1.idx,  y_complex);
    if (m2.isValid()) A_c.addValue(m2.idx, m2.idx,  y_complex);

    if (m1.isValid() && m2.isValid()) {
        A_c.addValue(m1.idx, m2.idx, -y_complex);
        A_c.addValue(m2.idx, m1.idx, -y_complex);
    }
}

void Diode::getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns, MnaIndex extraIdx) const {
    (void)extraIdx; // Diodes are non-linear 2-terminal devices and do not request extra rows

    if (nodes.size() < 2) return;

    MnaIndex m1 = MnaIndex::fromNetNode(nodes[0]);
    MnaIndex m2 = MnaIndex::fromNetNode(nodes[1]);

    if (m1.isValid()) patterns.push_back({m1, m1});
    if (m2.isValid()) patterns.push_back({m2, m2});

    if (m1.isValid() && m2.isValid()) {
        patterns.push_back({m1, m2});
        patterns.push_back({m2, m1});
    }
}

double Diode::calculateVoltageDrop(const std::vector<double>& state) const {
    if (nodes.size() < 2) return 0.0;
    return getSafeVoltage(state, nodes[0]) - getSafeVoltage(state, nodes[1]);
}

double Diode::calculateCurrent(const SimulationContext& simCtx, MnaIndex extraIdx) const {
    (void)extraIdx;
    if (simCtx.state.empty() || nodes.size() < 2) {
        return 0.0;
    }

    double vd = calculateVoltageDrop(simCtx.state);
    double safeIs = std::max(Is, 1e-20);
    double vt = N * V_THERMAL_300K;

    // 1. Analytical SPICE Voltage Clamping Layer (Prevents numeric floating-point overflow)
    double vcrit = vt * std::log(vt / (std::sqrt(2.0) * safeIs));
    if (vd > vcrit) {
        vd = vcrit + (vd - vcrit) / (1.0 + (vd - vcrit) / vt);
    }

    // 2. Core Shockley Equation evaluation with numerical exponent bounding
    double id = safeIs * (std::exp(std::min(vd / vt, 40.0)) - 1.0);

    // 3. Parallel Gmin conductance baseline injection to match the linear matrix stamp() exactly
    constexpr double gmin = 1e-12;
    id += gmin * vd;

    return id;
}


// ============================================================================
// 3. CAPACITOR
// ============================================================================

Capacitor::Capacitor(const std::string& name, const std::vector<NetNode>& nodes, double c)
    : Component(name, ComponentType::CAPACITOR, nodes), capacitance(c) {}
// mathematically correct, good and robust and no structural bugs 
void Capacitor::stamp(MnaMatrixContext& ctx) {
    if (nodes.size() < 2) return;

    NetNode n1 = nodes[0];
    NetNode n2 = nodes[1];

    MnaIndex m1 = MnaIndex::fromNetNode(n1);
    MnaIndex m2 = MnaIndex::fromNetNode(n2);

    // TR-BDF2 optimal step fractional spacing constant: gamma = 2.0 - sqrt(2.0)
    const double gamma = 2.0 - std::sqrt(2.0);
    double geq = 0.0;
    double ieq = 0.0;

    // 1. DC Operating Point Analysis (Order 0): Capacitors act as open circuits with Gmin shunt
    if (ctx.orderStep == 0 || ctx.dt <= 0.0) {
        constexpr double gmin = 1e-12;
        if (m1.isValid()) ctx.G.addValue(m1.idx, m1.idx,  gmin);
        if (m2.isValid()) ctx.G.addValue(m2.idx, m2.idx,  gmin);
        if (m1.isValid() && m2.isValid()) {
            ctx.G.addValue(m1.idx, m2.idx, -gmin);
            ctx.G.addValue(m2.idx, m1.idx, -gmin);
        }
        return;
    }

    // 2. Transient Integration Companion Model Evaluation with bounds-checked NetNode lookups
    if (ctx.orderStep == 1) {
        // TR Stage Companion Model: geq = 2C / (gamma * dt)
        geq = (2.0 * capacitance) / (gamma * ctx.dt);

        double v_n = 0.0;
        if (ctx.historicStepV) {
            v_n = getSafeVoltage(ctx.historicStepV, n1) - getSafeVoltage(ctx.historicStepV, n2);
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
            v_n = getSafeVoltage(ctx.historicStepV, n1) - getSafeVoltage(ctx.historicStepV, n2);
        }
        if (ctx.midpointStepV) {
            v_nc = getSafeVoltage(ctx.midpointStepV, n1) - getSafeVoltage(ctx.midpointStepV, n2);
        }
        // Exact discrete multi-step history integration mapping
        ieq = -capacitance * (c2 * v_nc - (c2 - c1) * v_n) / ctx.dt;
    }

    // 3. Conductance Admittance Matrix Stamps using validated 0-based MnaIndices
    if (m1.isValid()) ctx.G.addValue(m1.idx, m1.idx,  geq);
    if (m2.isValid()) ctx.G.addValue(m2.idx, m2.idx,  geq);

    if (m1.isValid() && m2.isValid()) {
        ctx.G.addValue(m1.idx, m2.idx, -geq);
        ctx.G.addValue(m2.idx, m1.idx, -geq);
    }

    // 4. Bounds-checked RHS injection (current exits positive terminal, enters negative terminal)
    addSafeRhs(ctx.B, m1, -ieq);
    addSafeRhs(ctx.B, m2,  ieq);
}

void Capacitor::updateHistory(const MnaMatrixContext& ctx) {
    if (nodes.size() < 2 || !ctx.currentIterV) return;

    // Safely sample instant node voltages via the bounds-checked NetNode helper
    double v = getSafeVoltage(ctx.currentIterV, nodes[0]) - getSafeVoltage(ctx.currentIterV, nodes[1]);

    if (ctx.orderStep == 1) {
        // Store intermediate TR stage charge integration state point (v_{n+gamma})
        q_hist_nc = capacitance * v;
    } else if (ctx.orderStep == 2) {
        // Commit final full BDF2 time step charge state point (v_{n+1})
        q_hist_n = capacitance * v;
    }
}

void Capacitor::getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns, MnaIndex extraIdx) const {
    (void)extraIdx; // Passive two-terminal components do not request extra rows
    if (nodes.size() < 2) return;

    MnaIndex m1 = MnaIndex::fromNetNode(nodes[0]);
    MnaIndex m2 = MnaIndex::fromNetNode(nodes[1]);

    if (m1.isValid()) patterns.push_back({m1, m1});
    if (m2.isValid()) patterns.push_back({m2, m2});

    if (m1.isValid() && m2.isValid()) {
        patterns.push_back({m1, m2});
        patterns.push_back({m2, m1});
    }
}

double Capacitor::calculateVoltageDrop(const std::vector<double>& state) const {
    if (nodes.size() < 2) return 0.0;
    return getSafeVoltage(state, nodes[0]) - getSafeVoltage(state, nodes[1]);
}

double Capacitor::calculateCurrent(const SimulationContext& simCtx, MnaIndex extraIdx) const {
    (void)extraIdx;
    if (simCtx.state.empty() || nodes.size() < 2) {
        return 0.0;
    }

    // During pure static DC Operating Point analysis (Order 0), an ideal capacitor acts as an open circuit
    if (simCtx.orderStep == 0 || simCtx.dt <= 0.0) {
        constexpr double gmin = 1e-12;
        return gmin * calculateVoltageDrop(simCtx.state);
    }

    // Extract natively solved branch current if assigned to an MNA extra row
    MnaIndex idx = this->assignedExtraRowIdx;
    if (idx.isValid() && static_cast<size_t>(idx.idx) < simCtx.state.size()) {
        return simCtx.state[idx.idx];
    }

    double v_now = calculateVoltageDrop(simCtx.state);
    double safeCap = (std::abs(capacitance) > 1e-24) ? capacitance : 1e-24;

    // TR-BDF2 Step Fractional Parameter: gamma = 2.0 - sqrt(2.0)
    const double gamma = 2.0 - std::sqrt(2.0);
    double geq = 0.0;
    double ieq = 0.0;

    if (simCtx.orderStep == 1) {
        // Trapezoidal Integration Substep
        geq = (2.0 * safeCap) / (gamma * simCtx.dt);
        double v_n = q_hist_n / safeCap;
        ieq = -geq * v_n - (q_hist_n / (gamma * simCtx.dt * 0.5));
    } else if (simCtx.orderStep == 2) {
        // BDF2 Integration Step Phase
        double c1 = (2.0 - gamma) / (1.0 - gamma);
        double c2 = 1.0 / (gamma * (1.0 - gamma));
        geq = c1 * safeCap / simCtx.dt;
        ieq = -(c2 * q_hist_nc - (c2 - c1) * q_hist_n) / simCtx.dt;
    }

    return (geq * v_now) + ieq;
}

void Capacitor::stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                        double omega, const std::vector<double>& dcState) {
    (void)B_c;
    (void)dcState; // Ideal linear capacitors do not modulate over the DC operating point state

    if (nodes.size() < 2) {
        return;
    }

    // 1. Translate 1-based NetNodes into 0-based MnaIndices (Ground automatically becomes invalid)
    MnaIndex m1 = MnaIndex::fromNetNode(nodes[0]);
    MnaIndex m2 = MnaIndex::fromNetNode(nodes[1]);

    // 2. Compute Admittance Complex Phasor Representation: Y = 0 + j * omega * C
    std::complex<double> j_omega_c(0.0, omega * capacitance);

    // Enforce Gmin shunt baseline at near-DC or zero frequency to prevent singular matrices
    if (std::abs(omega) < 1e-20) {
        constexpr double gmin = 1e-12;
        j_omega_c = std::complex<double>(gmin, 0.0);
    }

    // 3. Inject complex admittance stencil using validated 0-based MnaIndices
    if (m1.isValid()) {
        A_c.addValue(m1.idx, m1.idx,  j_omega_c);
    }
    if (m2.isValid()) {
        A_c.addValue(m2.idx, m2.idx,  j_omega_c);
    }

    if (m1.isValid() && m2.isValid()) {
        A_c.addValue(m1.idx, m2.idx, -j_omega_c);
        A_c.addValue(m2.idx, m1.idx, -j_omega_c);
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
Inductor::Inductor(const std::string& name, const std::vector<NetNode>& nodes, double l)
    : Component(name, ComponentType::INDUCTOR, nodes), inductance(l) {}

void Inductor::stamp(MnaMatrixContext& ctx) {
    // 1. Structural Layer Guards: Ensure nodes and extra row references are safely allocated
    if (nodes.size() < 2) return;

    MnaIndex idx = ctx.assignedExtraIdx.isValid() ? ctx.assignedExtraIdx : this->assignedExtraRowIdx;
    if (!idx.isValid() || static_cast<size_t>(idx.idx) >= ctx.B.size()) return;

    NetNode n1 = nodes[0]; // Positive terminal
    NetNode n2 = nodes[1]; // Negative terminal

    MnaIndex m1 = MnaIndex::fromNetNode(n1);
    MnaIndex m2 = MnaIndex::fromNetNode(n2);

    // TR-BDF2 optimal step fractional spacing constant: gamma = 2.0 - sqrt(2.0)
    const double gamma = 2.0 - std::sqrt(2.0);

    // 2. Nodal Interconnection Matrix Stamps using validated 0-based MnaIndices
    if (m1.isValid()) {
        ctx.G.addValue(m1.idx, idx.idx,  1.0);
        ctx.G.addValue(idx.idx, m1.idx,  1.0);
    }
    if (m2.isValid()) {
        ctx.G.addValue(m2.idx, idx.idx, -1.0);
        ctx.G.addValue(idx.idx, m2.idx, -1.0);
    }

    // 3. DC Operating Point Analysis (Order 0): Inductors act as pure structural short circuits
    if (ctx.orderStep == 0 || ctx.dt <= 0.0) {
        ctx.G.addValue(idx.idx, idx.idx, 0.0);
        setSafeRhs(ctx.B, idx, 0.0);
        return;
    }

    // 4. Transient Integration Companion Model Evaluation with synchronized signs
    if (ctx.orderStep == 1) {
        // TR Stage Companion Model: Req = 2L / (gamma * dt)
        double req = (2.0 * inductance) / (gamma * ctx.dt);
        ctx.G.addValue(idx.idx, idx.idx, -req);

        // Extract historical voltage drop across terminals via bounds-checked NetNode helper
        double v_past = 0.0;
        if (ctx.historicStepV) {
            v_past = getSafeVoltage(ctx.historicStepV, n1) - getSafeVoltage(ctx.historicStepV, n2);
        }
        // Thévenin RHS injection: B[idx] = v_n - Req * i_n
        setSafeRhs(ctx.B, idx, v_past - req * i_hist_n);

    } else if (ctx.orderStep == 2) {
        // BDF2 Composite Stage Companion Model: Multi-step time grid integration coefficients
        double c_factor = (2.0 - gamma) / (1.0 - gamma);
        double c2 = 1.0 / (gamma * (1.0 - gamma));
        double req = c_factor * inductance / ctx.dt;
        ctx.G.addValue(idx.idx, idx.idx, -req);

        // Compute discrete multi-step history current injection tracking framework
        double i_hist_term = (c2 * i_hist_nc - (c2 - c_factor) * i_hist_n);
        // Thévenin RHS injection: B[idx] = - (L / dt) * [c2 * i_{n+gamma} - (c2 - c_factor) * i_n]
        setSafeRhs(ctx.B, idx, -(inductance / ctx.dt) * i_hist_term);
    }
}

void Inductor::updateHistory(const MnaMatrixContext& ctx) {
    MnaIndex idx = ctx.assignedExtraIdx.isValid() ? ctx.assignedExtraIdx : this->assignedExtraRowIdx;
    if (!idx.isValid() || !ctx.currentIterV || static_cast<size_t>(idx.idx) >= ctx.currentIterV->size()) {
        return;
    }

    // Safely sample the solved branch current from the MNA state vector
    double current = getSafeExtraValue(ctx.currentIterV, idx);

    if (ctx.orderStep == 1) {
        // Store intermediate TR stage branch current integration state point (i_{n+gamma})
        i_hist_nc = current;
    } else if (ctx.orderStep == 2) {
        // Commit final full BDF2 time step branch current state point (i_{n+1})
        i_hist_n = current;
    }
}

void Inductor::getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns, MnaIndex extraIdx) const {
    if (nodes.size() < 2) return;

    MnaIndex brRow = extraIdx.isValid() ? extraIdx : this->assignedExtraRowIdx;
    if (!brRow.isValid()) return;

    // Self-inductance internal structural matrix grid point allocation
    patterns.push_back({brRow, brRow});

    MnaIndex m1 = MnaIndex::fromNetNode(nodes[0]);
    MnaIndex m2 = MnaIndex::fromNetNode(nodes[1]);

    // Build the static branch coupling pattern maps using validated 0-based MnaIndices
    if (m1.isValid()) {
        patterns.push_back({m1, brRow});
        patterns.push_back({brRow, m1});
    }
    if (m2.isValid()) {
        patterns.push_back({m2, brRow});
        patterns.push_back({brRow, m2});
    }
}

double Inductor::calculateCurrent(const SimulationContext& simCtx, MnaIndex extraIdx) const {
    MnaIndex targetIdx = extraIdx.isValid() ? extraIdx : this->assignedExtraRowIdx;

    // Direct Optimal Path: Read natively solved branch current variable from MNA solution state
    if (targetIdx.isValid() && static_cast<size_t>(targetIdx.idx) < simCtx.state.size()) {
        return getSafeExtraValue(simCtx.state, targetIdx);
    }

    // Algebraic Fallback: Used when probing parameters outside the synchronized solver loop context
    if (simCtx.state.empty()) {
        return 0.0;
    }
    if (simCtx.orderStep == 0 || simCtx.dt <= 0.0) {
        return i_hist_n;
    }

    double safeInductance = (std::abs(inductance) > 1e-24) ? inductance : 1e-24;
    const double gamma = 2.0 - std::sqrt(2.0);
    double v_now = calculateVoltageDrop(simCtx.state);

    if (simCtx.orderStep == 1) { // TR Companion Fallback Equation
        return i_hist_n + (v_now * (gamma * simCtx.dt) / (2.0 * safeInductance));
    } else if (simCtx.orderStep == 2) { // BDF2 Companion Fallback Equation
        double c_factor = (2.0 - gamma) / (1.0 - gamma);
        double c2 = 1.0 / (gamma * (1.0 - gamma));
        double ieq = (c2 * i_hist_nc - (c2 - c_factor) * i_hist_n);
        return ieq + (v_now * simCtx.dt / (c_factor * safeInductance));
    }
    return 0.0;
}

double Inductor::calculateVoltageDrop(const std::vector<double>& state) const {
    if (nodes.size() < 2) return 0.0;
    return getSafeVoltage(state, nodes[0]) - getSafeVoltage(state, nodes[1]);
}

void Inductor::stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                       double omega, const std::vector<double>& dcState) {
    (void)B_c;
    (void)dcState; // Linear ideal inductors do not modulate over the DC operating point state

    // 1. Structural Layer Guards: Ensure nodes and extra row references are safely allocated
    if (nodes.size() < 2) return;

    MnaIndex idx = getAssignedExtraIdx(); // 0-based branch current MnaIndex
    if (!idx.isValid()) return;

    MnaIndex m1 = MnaIndex::fromNetNode(nodes[0]); // Positive terminal
    MnaIndex m2 = MnaIndex::fromNetNode(nodes[1]); // Negative terminal

    // 2. Compute Complex Impedance of the Inductor: Z = 0 + j * omega * L
    std::complex<double> j_omega_l(0.0, omega * inductance);

    // Enforce numerical stability baseline series resistance at near-DC or zero frequency
    if (std::abs(omega) < 1e-20) {
        constexpr double rmin = 1e-12;
        j_omega_l = std::complex<double>(rmin, 0.0);
    }

    // 3. Inject branch-current coupling entries using validated 0-based MnaIndices
    if (m1.isValid()) {
        A_c.addValue(m1.idx,  idx.idx,  1.0);
        A_c.addValue(idx.idx, m1.idx,   1.0);
    }
    if (m2.isValid()) {
        A_c.addValue(m2.idx,  idx.idx, -1.0);
        A_c.addValue(idx.idx, m2.idx,  -1.0);
    }

    // Branch constitutive equation: V_{m1} - V_{m2} - (j * omega * L) * I_branch = 0
    A_c.addValue(idx.idx, idx.idx, -j_omega_l);
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
Mosfet::Mosfet(const std::string& name, const std::vector<NetNode>& nodes, bool nChannel)
    : Component(name, ComponentType::MOSFET, nodes), isNChannel(nChannel) {}

void Mosfet::stamp(MnaMatrixContext& ctx) {
    // 1. Structural Layer Guard: Ensure all 3 multi-terminal pins are present
    if (nodes.size() < 3) return;

    NetNode nD_working = nodes[0]; // Physical Drain NetNode
    NetNode nG_working = nodes[1]; // Physical Gate NetNode
    NetNode nS_working = nodes[2]; // Physical Source NetNode

    // 2. Safely sample contemporary absolute node voltages via bounds-checked NetNode helper
    double vD = 0.0;
    double vG = 0.0;
    double vS = 0.0;

    if (ctx.currentIterV) {
        vD = getSafeVoltage(ctx.currentIterV, nD_working);
        vG = getSafeVoltage(ctx.currentIterV, nG_working);
        vS = getSafeVoltage(ctx.currentIterV, nS_working);
    }

    // 3. Polarity Normalization (+1 for NMOS, -1 for PMOS)
    const double pol = isNChannel ? 1.0 : -1.0;
    vD *= pol;
    vG *= pol;
    vS *= pol;

    // 4. Handle Bidirectional Operating Modes (Symmetric Drain-Source Swapping)
    double vds = vD - vS;
    if (vds < 0.0) {
        std::swap(nD_working, nS_working);
        std::swap(vD, vS);
        vds = -vds;
    }

    // Compute gate overdrive voltage in the normalized working frame
    double vgs = vG - vS;
    double vov = vgs - vt0;

    double ids = 0.0;
    double gm  = 0.0;
    double gds = 1e-12; // Minimum floor conductance for matrix conditioning

    // 5. Piecewise Evaluation of the Shichman-Hodges Square-Law Model
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

    if (gds < 1e-12) gds = 1e-12;

    // Norton equivalent compensation current converted back to physical terminal polarity
    double ieq = pol * (ids - gm * vgs - gds * vds);

    // 6. Translate 1-based NetNodes into 0-based MnaIndices (trapping Ground automatically)
    MnaIndex mD = MnaIndex::fromNetNode(nD_working);
    MnaIndex mG = MnaIndex::fromNetNode(nG_working);
    MnaIndex mS = MnaIndex::fromNetNode(nS_working);

    // 7. Multi-Terminal Jacobian Conductance & Bounds-Checked RHS Injection
    if (mD.isValid()) {
        ctx.G.addValue(mD.idx, mD.idx,  gds);
        if (mG.isValid()) ctx.G.addValue(mD.idx, mG.idx,  gm);
        if (mS.isValid()) ctx.G.addValue(mD.idx, mS.idx, -(gm + gds));
        addSafeRhs(ctx.B, mD, -ieq); // Current leaving working Drain
    }

    if (mS.isValid()) {
        ctx.G.addValue(mS.idx, mS.idx,  (gm + gds));
        if (mD.isValid()) ctx.G.addValue(mS.idx, mD.idx, -gds);
        if (mG.isValid()) ctx.G.addValue(mS.idx, mG.idx, -gm);
        addSafeRhs(ctx.B, mS, ieq);   // Current entering working Source
    }
}
void Mosfet::stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                     double omega, const std::vector<double>& dcState) {
    (void)B_c; // Nonlinear devices do not load independent phasor source inputs

    // 1. Structural Layer Guard: Ensure the device contains all 3 hardware terminals
    if (nodes.size() < 3) return;

    const NetNode nD_phys = nodes[0];
    const NetNode nG_phys = nodes[1];
    const NetNode nS_phys = nodes[2];

    // 2. Safely resolve absolute node potentials from DC operating point via NetNode helper
    double vD_val = getSafeVoltage(dcState, nD_phys);
    double vG_val = getSafeVoltage(dcState, nG_phys);
    double vS_val = getSafeVoltage(dcState, nS_phys);

    // 3. Polarity normalization (+1 for NMOS, -1 for PMOS) & channel inversion check
    const double pol = isNChannel ? 1.0 : -1.0;
    double vD = pol * vD_val;
    double vG = pol * vG_val;
    double vS = pol * vS_val;

    bool is_swapped = false;
    double vds_internal = vD - vS;
    if (vds_internal < 0.0) {
        std::swap(vD, vS);
        vds_internal = -vds_internal;
        is_swapped = true;
    }
    double vgs_internal = vG - vS;

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

    // Determine effective NetNode mapping tags for proper geometric KCL routing
    NetNode nD_eff = is_swapped ? nS_phys : nD_phys;
    NetNode nS_eff = is_swapped ? nD_phys : nS_phys;

    // 6. Translate effective NetNodes into 0-based MnaIndices (trapping Ground automatically)
    MnaIndex mD_eff  = MnaIndex::fromNetNode(nD_eff);
    MnaIndex mG_phys = MnaIndex::fromNetNode(nG_phys);
    MnaIndex mS_eff  = MnaIndex::fromNetNode(nS_eff);

    // 7. Inject Complex Nodal Admittance Stencils using validated 0-based MnaIndices
    if (mD_eff.isValid()) {
        A_c.addValue(mD_eff.idx, mD_eff.idx, std::complex<double>(gds, 0.0) + y_gd);
        if (mG_phys.isValid()) A_c.addValue(mD_eff.idx, mG_phys.idx, std::complex<double>(gm, 0.0) - y_gd);
        if (mS_eff.isValid())  A_c.addValue(mD_eff.idx, mS_eff.idx,  std::complex<double>(-(gm + gds), 0.0));
    }

    if (mG_phys.isValid()) {
        A_c.addValue(mG_phys.idx, mG_phys.idx, y_gs + y_gd);
        if (mD_eff.isValid())  A_c.addValue(mG_phys.idx, mD_eff.idx,  -y_gd);
        if (mS_eff.isValid())  A_c.addValue(mG_phys.idx, mS_eff.idx,  -y_gs);
    }

    if (mS_eff.isValid()) {
        A_c.addValue(mS_eff.idx, mS_eff.idx, std::complex<double>(gm + gds, 0.0) + y_gs);
        if (mD_eff.isValid())  A_c.addValue(mS_eff.idx, mD_eff.idx,  std::complex<double>(-gds, 0.0));
        if (mG_phys.isValid()) A_c.addValue(mS_eff.idx, mG_phys.idx, std::complex<double>(-gm, 0.0) - y_gs);
    }
}

double Mosfet::calculateCurrent(const SimulationContext& simCtx, MnaIndex extraIdx) const {
    (void)extraIdx;
    if (simCtx.state.empty() || nodes.size() < 3) {
        return 0.0;
    }

    // 1. Safe extraction of raw absolute node potentials via NetNode helper
    double vD = getSafeVoltage(simCtx.state, nodes[0]);
    double vG = getSafeVoltage(simCtx.state, nodes[1]);
    double vS = getSafeVoltage(simCtx.state, nodes[2]);

    // 2. Polarity normalization (+1 for NMOS, -1 for PMOS)
    const double pol = isNChannel ? 1.0 : -1.0;
    vD *= pol;
    vG *= pol;
    vS *= pol;

    bool reversedCurrent = !isNChannel;

    // 3. Handle symmetric bidirectional drain-source swapping
    if ((vD - vS) < 0.0) {
        std::swap(vD, vS);
        reversedCurrent = !reversedCurrent;
    }

    // 4. Compute localized Shichman-Hodges operating parameters
    double vgs = vG - vS;
    double vds = vD - vS;
    double vov = vgs - vt0;

    if (vov <= 0.0) {
        return 0.0;
    }

    double ids = 0.0;
    double lambdaTerm = 1.0 + lambda * vds;

    if (vds < vov) {
        ids = beta * (vov * vds - 0.5 * vds * vds) * lambdaTerm;
    } else {
        ids = 0.5 * beta * vov * vov * lambdaTerm;
    }

    return reversedCurrent ? -ids : ids;
}

double Mosfet::calculateVoltageDrop(const std::vector<double>& state) const {
    if (nodes.size() < 3) return 0.0;
    // Vds = V_drain - V_source -> nodes[0] - nodes[2]
    return getSafeVoltage(state, nodes[0]) - getSafeVoltage(state, nodes[2]);
}

void Mosfet::getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns, MnaIndex extraIdx) const {
    (void)extraIdx; // MOSFETs are non-linear three-terminal devices and do not request extra rows

    // Structural Netlist Validation Guard: Ensure the parser has populated all 3 mandatory terminals
    if (nodes.size() < 3) {
        return;
    }

    // Loop through the 3x3 internal device connectivity matrix stencil mappings
    // 0 = Drain, 1 = Gate, 2 = Source
    for (int i : {0, 1, 2}) {
        MnaIndex mi = MnaIndex::fromNetNode(nodes[i]);
        if (!mi.isValid()) continue;

        for (int j : {0, 1, 2}) {
            MnaIndex mj = MnaIndex::fromNetNode(nodes[j]);
            if (mj.isValid()) {
                patterns.push_back({mi, mj});
            }
        }
    }
}

// ============================================================================
// 6. ADVANCED BJT
// ============================================================================

VbicBjtAdvanced::VbicBjtAdvanced(const std::string& name, const std::vector<NetNode>& nodes, bool npnVariant)
    : Component(name, ComponentType::BJT, nodes), isNpn(npnVariant) {}

void VbicBjtAdvanced::stamp(MnaMatrixContext& ctx) {
    // 1. Structural Layer Guard: Ensure the device possesses all 3 mandatory terminals
    if (nodes.size() < 3) return;

    const NetNode nC = nodes[0]; // Physical Collector NetNode
    const NetNode nB = nodes[1]; // Physical Base NetNode
    const NetNode nE = nodes[2]; // Physical Emitter NetNode

    // 2. Safely resolve node potentials via bounds-checked NetNode helpers
    double vC = 0.0, vB = 0.0, vE = 0.0;
    double vold_C = 0.0, vold_B = 0.0, vold_E = 0.0;

    if (ctx.currentIterV) {
        vC = getSafeVoltage(ctx.currentIterV, nC);
        vB = getSafeVoltage(ctx.currentIterV, nB);
        vE = getSafeVoltage(ctx.currentIterV, nE);
    }
    if (ctx.historicStepV) {
        vold_C = getSafeVoltage(ctx.historicStepV, nC);
        vold_B = getSafeVoltage(ctx.historicStepV, nB);
        vold_E = getSafeVoltage(ctx.historicStepV, nE);
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
    double ibe = (ibf / 100.0) + gmin * vbe_lim; // Static forward beta = 100
    double ibc = (ibr / 50.0) + gmin * vbc_lim;   // Static reverse beta = 50

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

    // True External Currents
    double ic_ext = sig * (ice - ibc);
    double ib_ext = sig * (ibe + ibc);
    double ie_ext = -(ic_ext + ib_ext);

    // RHS Companion Vector Calculations via Taylor Expansion
    double ieq_c = -ic_ext + (final_gCC * vC + final_gCB * vB + final_gCE * vE);
    double ieq_b = -ib_ext + (final_gBC * vC + final_gBB * vB + final_gBE * vE);
    double ieq_e = -ie_ext + (final_gEC * vC + final_gEB * vB + final_gEE * vE);

    // 7. Translate 1-based NetNodes into 0-based MnaIndices (trapping Ground automatically)
    MnaIndex mC = MnaIndex::fromNetNode(nC);
    MnaIndex mB = MnaIndex::fromNetNode(nB);
    MnaIndex mE = MnaIndex::fromNetNode(nE);

    // 8. Multi-Terminal Jacobian Matrix & Bounds-Checked RHS Stamping Pipeline
    if (mC.isValid()) {
        ctx.G.addValue(mC.idx, mC.idx, final_gCC);
        if (mB.isValid()) ctx.G.addValue(mC.idx, mB.idx, final_gCB);
        if (mE.isValid()) ctx.G.addValue(mC.idx, mE.idx, final_gCE);
        addSafeRhs(ctx.B, mC, ieq_c);
    }
    if (mB.isValid()) {
        if (mC.isValid()) ctx.G.addValue(mB.idx, mC.idx, final_gBC);
        ctx.G.addValue(mB.idx, mB.idx, final_gBB);
        if (mE.isValid()) ctx.G.addValue(mB.idx, mE.idx, final_gBE);
        addSafeRhs(ctx.B, mB, ieq_b);
    }
    if (mE.isValid()) {
        if (mC.isValid()) ctx.G.addValue(mE.idx, mC.idx, final_gEC);
        if (mB.isValid()) ctx.G.addValue(mE.idx, mB.idx, final_gEB);
        ctx.G.addValue(mE.idx, mE.idx, final_gEE);
        addSafeRhs(ctx.B, mE, ieq_e);
    }
}

void VbicBjtAdvanced::stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                              double omega, const std::vector<double>& dcState) {
    (void)B_c; // Small-signal independent AC elements handle B_c vector modifications

    // 1. Structural Layer Guard: Ensure the device possesses all 3 mandatory terminals
    if (nodes.size() < 3) return;

    const NetNode nC = nodes[0];
    const NetNode nB = nodes[1];
    const NetNode nE = nodes[2];

    // 2. Safely resolve absolute node potentials from DC operating point via NetNode helper
    double vC_val = getSafeVoltage(dcState, nC);
    double vB_val = getSafeVoltage(dcState, nB);
    double vE_val = getSafeVoltage(dcState, nE);

    // 3. Determine intrinsic internal voltage orientations (PNP inversion mapping)
    double sig = isNpn ? 1.0 : -1.0;

    double vbe = sig * (vB_val - vE_val);
    double vbc = sig * (vB_val - vC_val);

    double vt  = Nf * V_THERMAL_300K;
    double vtr = Nr * V_THERMAL_300K;
    double safe_Is = std::max(Is, 1e-20);

    // 4. DC Operating Point Formulations & Safe Exponential Clamping
    double exp_vbe = std::exp(std::min(vbe / vt, 40.0));
    double exp_vbc = std::exp(std::min(vbc / vtr, 40.0));

    double ibf = safe_Is * (exp_vbe - 1.0);
    double ibr = safe_Is * (exp_vbc - 1.0);
    double qb = 1.0;
    double dqb_dvbe = 0.0;

    double dibf_dvbe = (safe_Is / vt) * exp_vbe;
    double dibr_dvbc = (safe_Is / vtr) * exp_vbc;

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

    double itz = (ibf - ibr) / qb;
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

    // 5. Formal Node Voltage Partial Derivative Conversions (sig^2 = +1 for both NPN and PNP)
    double dIc_dvC = -diC_dvbc;
    double dIc_dvB =  (diC_dvbe + diC_dvbc);
    double dIc_dvE = -diC_dvbe;

    double dIb_dvC = -diB_dvbc;
    double dIb_dvB =  (diB_dvbe + diB_dvbc);
    double dIb_dvE = -diB_dvbe;

    double dIe_dvC = -diE_dvbc;
    double dIe_dvB =  (diE_dvbe + diE_dvbc);
    double dIe_dvE = -diE_dvbe;

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

    constexpr double gmin = 1e-12; // Stabilizing shunt conductance

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

    // 8. Translate 1-based NetNodes into 0-based MnaIndices (trapping Ground automatically)
    MnaIndex mC = MnaIndex::fromNetNode(nC);
    MnaIndex mB = MnaIndex::fromNetNode(nB);
    MnaIndex mE = MnaIndex::fromNetNode(nE);

    // Complex Admittance Jacobian Matrix Injections
    if (mC.isValid()) {
        A_c.addValue(mC.idx, mC.idx, gCC);
        if (mB.isValid()) A_c.addValue(mC.idx, mB.idx, gCB);
        if (mE.isValid()) A_c.addValue(mC.idx, mE.idx, gCE);
    }
    if (mB.isValid()) {
        if (mC.isValid()) A_c.addValue(mB.idx, mC.idx, gBC);
        A_c.addValue(mB.idx, mB.idx, gBB);
        if (mE.isValid()) A_c.addValue(mB.idx, mE.idx, gBE);
    }
    if (mE.isValid()) {
        if (mC.isValid()) A_c.addValue(mE.idx, mC.idx, gEC);
        if (mB.isValid()) A_c.addValue(mE.idx, mB.idx, gEB);
        A_c.addValue(mE.idx, mE.idx, gEE);
    }
}

double VbicBjtAdvanced::calculateCurrent(const SimulationContext& simCtx, MnaIndex extraIdx) const {
    (void)extraIdx;
    if (simCtx.state.empty() || nodes.size() < 3) {
        return 0.0;
    }

    // 1. Safe extraction of raw absolute node potentials using NetNode getSafeVoltage helper
    const NetNode nC = nodes[0];
    const NetNode nB = nodes[1];
    const NetNode nE = nodes[2];

    double vC = getSafeVoltage(simCtx.state, nC);
    double vB = getSafeVoltage(simCtx.state, nB);
    double vE = getSafeVoltage(simCtx.state, nE);

    // 2. Perform Polarity Device Transformation Upfront (PNP Inversion Mapping Isolation)
    double sig = isNpn ? 1.0 : -1.0;

    double vbe = sig * (vB - vE);
    double vbc = sig * (vB - vC);

    double vt  = Nf * V_THERMAL_300K;
    double vtr = Nr * V_THERMAL_300K;
    double safe_Is = std::max(Is, 1e-20);

    // 3. Numerical Clamping Layer to eliminate exponential overflow risks during unstable solver sweeps
    double vcrit_f = vt * std::log(vt / (std::sqrt(2.0) * safe_Is));
    if (vbe > vcrit_f) {
        vbe = vcrit_f + (vbe - vcrit_f) / (1.0 + (vbe - vcrit_f) / vt);
    }
    double vcrit_r = vtr * std::log(vtr / (std::sqrt(2.0) * safe_Is));
    if (vbc > vcrit_r) {
        vbc = vcrit_r + (vbc - vcrit_r) / (1.0 + (vbc - vcrit_r) / vtr);
    }

    // 4. Core Transport & Diode Formulations matching stamp() exactly
    double exp_vbe = std::exp(std::min(vbe / vt, 40.0));
    double exp_vbc = std::exp(std::min(vbc / vtr, 40.0));

    double ibf = safe_Is * (exp_vbe - 1.0);
    double ibr = safe_Is * (exp_vbc - 1.0);

    // 5. High-Injection Gummel-Poon Knee / Base-Charge (qb) Evaluation
    double safe_ibf = std::max(0.0, ibf);
    double arg_q2 = 1.0 + 4.0 * (safe_ibf / (Ikf > 1e-20 ? Ikf : 1e-20));
    double qb = 1.0;
    if (arg_q2 > 0.0) {
        qb = 0.5 * (1.0 + std::sqrt(arg_q2));
    }
    if (qb < 1e-3) qb = 1e-3; // Numerical floor block to prevent divide-by-zero anomaly

    double itz = (ibf - ibr) / qb;
    double ib_r_total = ibr / 50.0; // Static reverse beta = 50

    // 6. Total Internal Collector Device Current Formulation
    double ic_internal = itz - ib_r_total;

    // Add parallel Gmin shunt factor for consistency with matrix stamps
    constexpr double gmin = 1e-12;
    ic_internal += gmin * (sig * (vC - vE));

    // Returns the true external collector branch current tracking polarization
    return sig * ic_internal;
}

double VbicBjtAdvanced::calculateVoltageDrop(const std::vector<double>& state) const {
    if (nodes.size() < 3) return 0.0;
    // Vce = V_collector - V_emitter -> nodes[0] - nodes[2]
    return getSafeVoltage(state, nodes[0]) - getSafeVoltage(state, nodes[2]);
}

void VbicBjtAdvanced::getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns, MnaIndex extraIdx) const {
    (void)extraIdx; // BJTs are nonlinear three-terminal devices and do not request extra rows

    // Structural Netlist Validation Guard: Ensure all 3 mandatory terminals exist
    if (nodes.size() < 3) {
        return;
    }

    // Loop through the 3x3 internal device connectivity matrix stencil mappings
    // 0 = Collector, 1 = Base, 2 = Emitter
    for (int i = 0; i < 3; ++i) {
        MnaIndex mi = MnaIndex::fromNetNode(nodes[i]);
        if (!mi.isValid()) continue;

        for (int j = 0; j < 3; ++j) {
            MnaIndex mj = MnaIndex::fromNetNode(nodes[j]);
            if (mj.isValid()) {
                patterns.push_back({mi, mj});
            }
        }
    }
}

// ============================================================================
// 7. INDEPENDENT VOLTAGE SOURCE
// ============================================================================

VoltageSource::VoltageSource(const std::string& name, const std::vector<NetNode>& nodes, double v)
    : Component(name, ComponentType::VOLTAGE_SOURCE, nodes), value(v) {}

void VoltageSource::stamp(MnaMatrixContext& ctx) {
    // 1. Structural Layer Guards: Ensure nodes and extra row references are safely allocated
    if (nodes.size() < 2) return;

    MnaIndex idx = ctx.assignedExtraIdx.isValid() ? ctx.assignedExtraIdx : this->assignedExtraRowIdx;
    if (!idx.isValid() || static_cast<size_t>(idx.idx) >= ctx.B.size()) return;

    MnaIndex m1 = MnaIndex::fromNetNode(nodes[0]); // Positive terminal
    MnaIndex m2 = MnaIndex::fromNetNode(nodes[1]); // Negative terminal

    // 2. Nodal Interconnection Matrix Stamps using validated 0-based MnaIndices
    if (m1.isValid()) {
        ctx.G.addValue(m1.idx,  idx.idx,  1.0);
        ctx.G.addValue(idx.idx, m1.idx,   1.0);
    }
    if (m2.isValid()) {
        ctx.G.addValue(m2.idx,  idx.idx, -1.0);
        ctx.G.addValue(idx.idx, m2.idx,  -1.0);
    }

    // Hard allocation for the branch diagonal grid point to secure solver structural integrity
    ctx.G.addValue(idx.idx, idx.idx, 0.0);

    // 3. Bounds-checked absolute assignment to the branch constitutive RHS row:
    // V_{m1} - V_{m2} = value
    setSafeRhs(ctx.B, idx, value);
}

void VoltageSource::getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns, MnaIndex extraIdx) const {
    if (nodes.size() < 2) return;

    MnaIndex brRow = extraIdx.isValid() ? extraIdx : this->assignedExtraRowIdx;
    if (!brRow.isValid()) return;

    // Allocate self-coupling branch diagonal matrix grid point entry
    patterns.push_back({brRow, brRow});

    MnaIndex m1 = MnaIndex::fromNetNode(nodes[0]);
    MnaIndex m2 = MnaIndex::fromNetNode(nodes[1]);

    // Build static connection pattern maps using validated 0-based MnaIndices
    if (m1.isValid()) {
        patterns.push_back({m1, brRow});
        patterns.push_back({brRow, m1});
    }
    if (m2.isValid()) {
        patterns.push_back({m2, brRow});
        patterns.push_back({brRow, m2});
    }
}

double VoltageSource::calculateCurrent(const SimulationContext& simCtx, MnaIndex extraIdx) const {
    MnaIndex targetIdx = extraIdx.isValid() ? extraIdx : this->assignedExtraRowIdx;
    return getSafeExtraValue(simCtx.state, targetIdx);
}

double VoltageSource::calculateVoltageDrop(const std::vector<double>& state) const {
    if (nodes.size() < 2) return 0.0;
    return getSafeVoltage(state, nodes[0]) - getSafeVoltage(state, nodes[1]);
}

void VoltageSource::stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                            double omega, const std::vector<double>& dcState) {
    (void)omega;
    (void)dcState; // Independent sources do not modulate based on the DC operating point state

    // 1. Structural Layer Guards: Ensure nodes and extra row references are safely allocated
    if (nodes.size() < 2) return;

    MnaIndex idx = getAssignedExtraIdx(); // 0-based branch current MnaIndex
    if (!idx.isValid() || static_cast<size_t>(idx.idx) >= B_c.size()) return;

    MnaIndex m1 = MnaIndex::fromNetNode(nodes[0]); // Positive terminal
    MnaIndex m2 = MnaIndex::fromNetNode(nodes[1]); // Negative terminal

    // 2. Inject branch-current coupling entries using validated 0-based MnaIndices
    if (m1.isValid()) {
        A_c.addValue(m1.idx,  idx.idx,  1.0);
        A_c.addValue(idx.idx, m1.idx,   1.0);
    }
    if (m2.isValid()) {
        A_c.addValue(m2.idx,  idx.idx, -1.0);
        A_c.addValue(idx.idx, m2.idx,  -1.0);
    }

    // Enforce branch matrix structural short-circuit identity baseline
    A_c.addValue(idx.idx, idx.idx, std::complex<double>(0.0, 0.0));

    // 3. Convert Magnitude/Phase into Cartesian Phasor Representation (Real + j*Imag)
    std::complex<double> v_phasor = std::polar(acMagnitude, acPhase);

    // 4. Bounds-checked absolute assignment on B_c: V_{m1} - V_{m2} = v_phasor
    setSafeRhs(B_c, idx, v_phasor);
}

void VoltageSource::setAcParameters(double mag, double phase) {
    acMagnitude = mag;
    acPhase = phase;
}

void VoltageSource::allocateRows(int& globalExtraRowCounter) {
    Component::allocateRows(globalExtraRowCounter);
}

// ============================================================================
// 8. BEHAVIORAL LOGIC GATE BRIDGE
// ============================================================================

BehaviorLogicGateComponent::BehaviorLogicGateComponent(
    const std::string& name, const std::vector<NetNode>& nodes, GateType gateType,
    double vThresh, double vHigh, double vLow, double rOut)
    : Component(name, ComponentType::DIGITAL_IC, nodes), gate(gateType),
      thresholdVoltage(vThresh), outputHighV(vHigh), outputLowV(vLow), outputResistance(rOut) {}

void BehaviorLogicGateComponent::stamp(MnaMatrixContext& ctx) {
    // 1. Structural Layer Guards: Ensure nodes and extra row references are safely allocated
    if (nodes.empty()) return;

    MnaIndex idx = ctx.assignedExtraIdx.isValid() ? ctx.assignedExtraIdx : this->assignedExtraRowIdx;
    if (!idx.isValid() || static_cast<size_t>(idx.idx) >= ctx.B.size()) return;

    // A NOT gate requires at least 2 terminals (1 In, 1 Out). Multi-input gates require at least 3 (2 In, 1 Out).
    size_t minTerminals = (gate == GateType::NOT) ? 2 : 3;
    if (nodes.size() < minTerminals) {
        throw std::runtime_error("Digital IC error: [" + name + "] lacks mandatory terminal connection pairs.");
    }

    NetNode nOut = nodes.back(); // Output pin is mapped to the final index element
    double continuousLogicState = 0.0;
    constexpr double kTransitionSharpness = 4.0; // Defines the tanh voltage transition steepness slope

    // Analytical derivative trackers mapping: d(continuousLogicState) / dV_node
    double dLogic_dvA = 0.0;
    double dLogic_dvB = 0.0;

    // 2. Safely resolve local voltages via bounds-checked NetNode helper
    if (gate == GateType::NOT) {
        NetNode nIn = nodes[0];
        double vIn = ctx.currentIterV ? getSafeVoltage(ctx.currentIterV, nIn) : 0.0;

        double tanhArg = kTransitionSharpness * (vIn - thresholdVoltage);
        double tanhVal = std::tanh(tanhArg);
        double sech2 = 1.0 - tanhVal * tanhVal;

        continuousLogicState = 0.5 * (1.0 - tanhVal);
        dLogic_dvA = -0.5 * kTransitionSharpness * sech2;
    }
    else {
        NetNode nA = nodes[0];
        NetNode nB = nodes[1];
        double vA = ctx.currentIterV ? getSafeVoltage(ctx.currentIterV, nA) : 0.0;
        double vB = ctx.currentIterV ? getSafeVoltage(ctx.currentIterV, nB) : 0.0;

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
    double vA_curr = ctx.currentIterV ? getSafeVoltage(ctx.currentIterV, nodes[0]) : 0.0;
    double vB_curr = (gate != GateType::NOT && ctx.currentIterV) ? getSafeVoltage(ctx.currentIterV, nodes[1]) : 0.0;

    double ieq = targetV - (dVtarget_dvA * vA_curr + dVtarget_dvB * vB_curr);

    // 5. Translate 1-based NetNodes into 0-based MnaIndices (trapping Ground automatically)
    MnaIndex mOut = MnaIndex::fromNetNode(nOut);
    MnaIndex mInA = MnaIndex::fromNetNode(nodes[0]);

    if (mOut.isValid()) {
        ctx.G.addValue(mOut.idx, idx.idx,  1.0);
        ctx.G.addValue(idx.idx,  mOut.idx, 1.0);
    }

    // Inject cross-coupling transconductance derivative column elements into the MNA grid
    if (mInA.isValid()) {
        ctx.G.addValue(idx.idx, mInA.idx, -dVtarget_dvA);
    }
    if (gate != GateType::NOT) {
        MnaIndex mInB = MnaIndex::fromNetNode(nodes[1]);
        if (mInB.isValid()) {
            ctx.G.addValue(idx.idx, mInB.idx, -dVtarget_dvB);
        }
    }

    // Physical behavioral driving output series resistance mapping coordinate
    ctx.G.addValue(idx.idx, idx.idx, -outputResistance);

    // Bounds-checked absolute assignment to the branch constitutive RHS row
    setSafeRhs(ctx.B, idx, ieq);
}

double BehaviorLogicGateComponent::calculateCurrent(const SimulationContext& simCtx, MnaIndex extraIdx) const {
    // Determine the active 0-based branch matrix coordinate tracking MnaIndex
    MnaIndex targetIdx = extraIdx.isValid() ? extraIdx : this->assignedExtraRowIdx;

    // Extract the solved logic output branch current via bounds-checked MnaIndex helper
    return getSafeExtraValue(simCtx.state, targetIdx);
}

double BehaviorLogicGateComponent::calculateVoltageDrop(const std::vector<double>& state) const {
    if (nodes.empty()) {
        return 0.0;
    }

    // The output terminal is always mapped to the final NetNode element of the nodes vector
    NetNode nOut = nodes.back();

    // Returns the exact absolute output node voltage relative to Ground (0V)
    return getSafeVoltage(state, nOut);
}

void BehaviorLogicGateComponent::stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                                         double omega, const std::vector<double>& dcState) {
    (void)omega;
    (void)B_c; // Independent phasor load adjustments are handled externally by AC source elements

    // 1. Structural Layer Guards: Ensure nodes and extra row references are safely allocated
    if (nodes.empty()) return;

    MnaIndex idx = getAssignedExtraIdx(); // 0-based branch MnaIndex
    if (!idx.isValid()) return;

    size_t minTerminals = (gate == GateType::NOT) ? 2 : 3;
    if (nodes.size() < minTerminals) return;

    NetNode nOut = nodes.back(); // Output pin is mapped to the final NetNode element
    constexpr double kTransitionSharpness = 4.0;

    // Analytical derivative trackers mapping: d(continuousLogicState) / dV_node
    double dLogic_dvA = 0.0;
    double dLogic_dvB = 0.0;

    // 2. Safely sample converged absolute DC node potentials via NetNode getSafeVoltage helper
    if (gate == GateType::NOT) {
        NetNode nIn = nodes[0];
        double vIn = getSafeVoltage(dcState, nIn);

        double tanhArg = kTransitionSharpness * (vIn - thresholdVoltage);
        double tanhVal = std::tanh(tanhArg);
        double sech2 = 1.0 - tanhVal * tanhVal;

        dLogic_dvA = -0.5 * kTransitionSharpness * sech2;
    }
    else {
        NetNode nA = nodes[0];
        NetNode nB = nodes[1];
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

    // 3. Translate 1-based NetNodes into 0-based MnaIndices (trapping Ground automatically)
    MnaIndex mOut = MnaIndex::fromNetNode(nOut);
    MnaIndex mInA = MnaIndex::fromNetNode(nodes[0]);

    // 4. Inject Complex Nodal Admittance Stencils using validated 0-based MnaIndices
    if (mOut.isValid()) {
        A_c.addValue(mOut.idx, idx.idx,  std::complex<double>(1.0, 0.0));
        A_c.addValue(idx.idx,  mOut.idx, std::complex<double>(1.0, 0.0));
    }

    // Inject linearized small-signal coupling transconductances into Jacobian columns
    if (mInA.isValid()) {
        A_c.addValue(idx.idx, mInA.idx, std::complex<double>(-dVtarget_dvA, 0.0));
    }
    if (gate != GateType::NOT) {
        MnaIndex mInB = MnaIndex::fromNetNode(nodes[1]);
        if (mInB.isValid()) {
            A_c.addValue(idx.idx, mInB.idx, std::complex<double>(-dVtarget_dvB, 0.0));
        }
    }

    // Output load series resistance phasor coupling coordinate injection
    A_c.addValue(idx.idx, idx.idx, std::complex<double>(-outputResistance, 0.0));
}

void BehaviorLogicGateComponent::allocateRows(int& globalExtraRowCounter) {
    Component::allocateRows(globalExtraRowCounter);
}

void BehaviorLogicGateComponent::getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns, MnaIndex extraIdx) const {
    MnaIndex brRow = extraIdx.isValid() ? extraIdx : this->assignedExtraRowIdx;
    if (!brRow.isValid()) return;

    size_t minTerminals = (gate == GateType::NOT) ? 2 : 3;
    if (nodes.size() < minTerminals) return;

    // Dynamically locate the output pin (always the final NetNode element)
    NetNode nOut = nodes.back();
    MnaIndex mOut = MnaIndex::fromNetNode(nOut);

    patterns.push_back({brRow, brRow}); // Diagonal structural branch entry

    if (mOut.isValid()) {
        patterns.push_back({mOut, brRow});
        patterns.push_back({brRow, mOut});
    }

    MnaIndex mInA = MnaIndex::fromNetNode(nodes[0]);
    if (mInA.isValid()) {
        patterns.push_back({brRow, mInA});
    }

    if (gate != GateType::NOT) {
        MnaIndex mInB = MnaIndex::fromNetNode(nodes[1]);
        if (mInB.isValid()) {
            patterns.push_back({brRow, mInB});
        }
    }
}

