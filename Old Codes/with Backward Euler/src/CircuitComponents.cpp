#include "freyalab/CircuitComponents.h"
#include "freyalab/math.h"
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Component::Component(const std::string& name, ComponentType type, const std::vector<int>& nodes) 
    : name(name), type(type), nodes(nodes) {}

// --- Resistor ---
Resistor::Resistor(const std::string& name, int nodeA, int nodeB, double resistance) 
    : Component(name, ComponentType::RESISTOR, {nodeA, nodeB}), resistance(resistance) {}

void Resistor::stamp(MnaMatrixContext& ctx) 
{
	int nA = nodes[0], nB = nodes[1];
	double g = 1.0 / resistance;
	if (nA > 0) MathSolver::accumulateCsr(ctx.G, nA - 1, nA - 1, g);
	if (nB > 0) MathSolver::accumulateCsr(ctx.G, nB - 1, nB - 1, g);
	if (nA > 0 && nB > 0) 
	{
		MathSolver::accumulateCsr(ctx.G, nA - 1, nB - 1, -g);
		MathSolver::accumulateCsr(ctx.G, nB - 1, nA - 1, -g);
	}
}

void Resistor::getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx) const 
{
	int nA = nodes[0], nB = nodes[1];
	if (nA > 0) patterns.push_back({nA - 1, nA - 1});
	if (nB > 0) patterns.push_back({nB - 1, nB - 1});
	if (nA > 0 && nB > 0) 
	{
		patterns.push_back({nA - 1, nB - 1});
		patterns.push_back({nB - 1, nA - 1});
	}
}

double Resistor::calculateCurrent(const std::vector<double>& state, int extraIdx) const 
{
	return calculateVoltageDrop(state) / resistance;
}

double Resistor::calculateVoltageDrop(const std::vector<double>& state) const 
{
	double vA = (nodes[0] > 0) ? state[nodes[0] - 1] : 0.0;
	double vB = (nodes[1] > 0) ? state[nodes[1] - 1] : 0.0;
	return vA - vB;
}

// --- Voltage Source ---
VoltageSource::VoltageSource(const std::string& name, int nodePos, int nodeNeg, double dcVal) 
    : Component(name, ComponentType::VOLTAGE_SOURCE, {nodePos, nodeNeg}), 
      srcType(SourceType::DC), value1(dcVal), value2(0), delay(0), rise(0), fall(0), width(0), period(0) {}

VoltageSource::VoltageSource(const std::string& name, int nodePos, int nodeNeg, SourceType type, 
                             double v1, double v2, double td, double tr, double tf, double tw, double tp) 
    : Component(name, ComponentType::VOLTAGE_SOURCE, {nodePos, nodeNeg}), 
      srcType(type), value1(v1), value2(v2), delay(td), rise(tr), fall(tf), width(tw), period(tp) {}

double VoltageSource::getCurrentValue(double time) const 
{
	if (srcType == SourceType::DC) return value1;
	if (srcType == SourceType::AC) 
	{
		double transientFreq = 1000.0; 
		return value1 + value2 * std::sin(2.0 * M_PI * transientFreq * time);
	}
	if (srcType == SourceType::PULSE) 
	{
		if (time < delay) return value1;
		double tRel = std::fmod(time - delay, period);
		if (tRel < 0) tRel += period;
		if (tRel < rise) return value1 + (value2 - value1) * (tRel / rise);
		if (tRel < rise + width) return value2;
		if (tRel < rise + width + fall) return value2 - (value2 - value1) * ((tRel - rise - width) / fall);
		return value1;
	}
	return 0.0;
}

void VoltageSource::stamp(MnaMatrixContext& ctx) 
{
	int nPos = nodes[0], nNeg = nodes[1], k = ctx.extraIdx;
	if (nPos > 0) 
	{
		MathSolver::accumulateCsr(ctx.G, nPos - 1, k, 1.0);
		MathSolver::accumulateCsr(ctx.G, k, nPos - 1, 1.0);
	}
	if (nNeg > 0) 
	{
		MathSolver::accumulateCsr(ctx.G, nNeg - 1, k, -1.0);
		MathSolver::accumulateCsr(ctx.G, k, nNeg - 1, -1.0);
	}
	ctx.B[k] += getCurrentValue(ctx.time); 
}

void VoltageSource::getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx) const 
{
	int nPos = nodes[0], nNeg = nodes[1], k = extraIdx;
	if (nPos > 0) 
	{
		patterns.push_back({nPos - 1, k});
		patterns.push_back({k, nPos - 1});
	}
	if (nNeg > 0) 
	{
		patterns.push_back({nNeg - 1, k});
		patterns.push_back({k, nNeg - 1});
	}
}

double VoltageSource::calculateCurrent(const std::vector<double>& state, int extraIdx) const 
{
	if (extraIdx >= 0 && extraIdx < (int)state.size()) return state[extraIdx];
	return 0.0;
}

double VoltageSource::calculateVoltageDrop(const std::vector<double>& state) const 
{
	double vP = (nodes[0] > 0) ? state[nodes[0] - 1] : 0.0;
	double vN = (nodes[1] > 0) ? state[nodes[1] - 1] : 0.0;
	return vP - vN;
}

// --- Capacitor ---
Capacitor::Capacitor(const std::string& name, int nodeA, int nodeB, double capacitance) 
    : Component(name, ComponentType::CAPACITOR, {nodeA, nodeB}), capacitance(capacitance) {}

void Capacitor::stamp(MnaMatrixContext& ctx) 
{
	int nA = nodes[0], nB = nodes[1];

	// 1. Guard against static DC operating points (Ideal Open Circuit)
	if (ctx.dt <= 1e-20) 
	{
		return; 
	}
	double geq = capacitance / ctx.dt;
	double v_old = 0.0;
	if (nA > 0) v_old += ctx.historicStepV[nA - 1];
	if (nB > 0) v_old -= ctx.historicStepV[nB - 1];
    
	lastVdrop = v_old; 
	double ieq = geq * v_old;

	if (nA > 0) MathSolver::accumulateCsr(ctx.G, nA - 1, nA - 1, geq);
	if (nB > 0) MathSolver::accumulateCsr(ctx.G, nB - 1, nB - 1, geq);
	if (nA > 0 && nB > 0) 
	{
		MathSolver::accumulateCsr(ctx.G, nA - 1, nB - 1, -geq);
		MathSolver::accumulateCsr(ctx.G, nB - 1, nA - 1, -geq);
	}
	if (nA > 0) ctx.B[nA - 1] += ieq;
	if (nB > 0) ctx.B[nB - 1] -= ieq;
}

void Capacitor::getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx) const 
{
	int nA = nodes[0], nB = nodes[1];
	if (nA > 0) patterns.push_back({nA - 1, nA - 1});
	if (nB > 0) patterns.push_back({nB - 1, nB - 1});
	if (nA > 0 && nB > 0) 
	{
		patterns.push_back({nA - 1, nB - 1});
		patterns.push_back({nB - 1, nA - 1});
	}
}

double Capacitor::calculateCurrent(const std::vector<double>& state, int extraIdx) const 
{
	return capacitance * (calculateVoltageDrop(state) - lastVdrop);
}

double Capacitor::calculateVoltageDrop(const std::vector<double>& state) const 
{
	double vA = (nodes[0] > 0) ? state[nodes[0] - 1] : 0.0;
	double vB = (nodes[1] > 0) ? state[nodes[1] - 1] : 0.0;
	return vA - vB;
}

// --- Inductor ---
Inductor::Inductor(const std::string& name, int nodeA, int nodeB, double inductance) 
    : Component(name, ComponentType::INDUCTOR, {nodeA, nodeB}), inductance(inductance) {}

// --- Inductor Transient-Ready Stamp ---
void Inductor::stamp(MnaMatrixContext& ctx) 
{
	int nA = nodes[0], nB = nodes[1], k = ctx.extraIdx;

	// 1. Guard against static DC operating points (Ideal Short Circuit / 0V Source)
	if (ctx.dt <= 1e-20) 
	{
		if (nA > 0) 
		{
			MathSolver::accumulateCsr(ctx.G, nA - 1, k, 1.0);
			MathSolver::accumulateCsr(ctx.G, k, nA - 1, 1.0);
		}
		if (nB > 0) 
		{
			MathSolver::accumulateCsr(ctx.G, nB - 1, k, -1.0);
			MathSolver::accumulateCsr(ctx.G, k, nB - 1, -1.0);
		}
		// RHS vector keeps a 0.0 force value for short circuits, do not accumulate garbage.
		return;
	}

	// 2. Compute Backward Euler transient companion values
	double req = inductance / ctx.dt;
    
	// Extract past time-step branch current flowing through the inductor extra row index k
	double i_old = ctx.historicStepV[k];

	// 3. Load Conductance Matrix G
	if (nA > 0) 
	{
		MathSolver::accumulateCsr(ctx.G, nA - 1, k, 1.0);
		MathSolver::accumulateCsr(ctx.G, k, nA - 1, 1.0);
	}
	if (nB > 0) 
	{
		MathSolver::accumulateCsr(ctx.G, nB - 1, k, -1.0);
		MathSolver::accumulateCsr(ctx.G, k, nB - 1, -1.0);
	}
    
	// Internal companion branch resistance contribution
	MathSolver::accumulateCsr(ctx.G, k, k, -req);

	// 4. Load RHS Vector B
	 // Pushes equivalent memory history onto RHS vector safely
	ctx.B[k] = - req * i_old;
}

void Inductor::getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx) const 
{
	int nA = nodes[0], nB = nodes[1], k = extraIdx;
	if (nA > 0) 
	{
		patterns.push_back({nA - 1, k});
		patterns.push_back({k, nA - 1});
	}
	if (nB > 0) 
	{
		patterns.push_back({nB - 1, k});
		patterns.push_back({k, nB - 1});
	}
	patterns.push_back({k, k});
}

double Inductor::calculateCurrent(const std::vector<double>& state, int extraIdx) const 
{
	if (extraIdx >= 0 && extraIdx < (int)state.size()) return state[extraIdx];
	return 0.0;
}

double Inductor::calculateVoltageDrop(const std::vector<double>& state) const 
{
	double vA = (nodes[0] > 0) ? state[nodes[0] - 1] : 0.0;
	double vB = (nodes[1] > 0) ? state[nodes[1] - 1] : 0.0;
	return vA - vB;
}

// --- Diode ---
Diode::Diode(const std::string& name, int anode, int cathode, double Is, double Vt) 
    : Component(name, ComponentType::DIODE, {anode, cathode}), Is(Is), Vt(Vt) {}

void Diode::stamp(MnaMatrixContext& ctx) 
{
	int nA = nodes[0], nB = nodes[1];
    
	// 1. FIX: Extract node voltages strictly from the active NR iteration workspace (Option B Alignment)
	double vA = (nA > 0) ? ctx.currentIterV[nA - 1] : 0.0;
	double vB = (nB > 0) ? ctx.currentIterV[nB - 1] : 0.0;
	double vd_target = vA - vB;

	// 2. Bounding junction step variations via standard SPICE pnjlim parameters
	constexpr double Vcrit = 0.58; 
	double vd_lim = MathSolver::pnjlim(vd_target, vd_prev_iter, Vt, Vcrit);

	// 3. Evaluate Exponential Diode Core Mechanics
	double expVd = std::exp(vd_lim / Vt);
	double gd = (Is / Vt) * expVd;
	double id = Is * (expVd - 1.0);
    
	// 4. Formulate pristine companion source currents (Ieq = Iop - J * V)
	double ieq = id - gd * vd_lim;

	// 5. Stamp Polarized Conductances to CsrMatrix G
	if (nA > 0) MathSolver::accumulateCsr(ctx.G, nA - 1, nA - 1, gd);
	if (nB > 0) MathSolver::accumulateCsr(ctx.G, nB - 1, nB - 1, gd);
	if (nA > 0 && nB > 0) 
	{
		MathSolver::accumulateCsr(ctx.G, nA - 1, nB - 1, -gd);
		MathSolver::accumulateCsr(ctx.G, nB - 1, nA - 1, -gd);
	}

	// 6. Stamp companion current injections to RHS Vector B
	if (nA > 0) ctx.B[nA - 1] -= ieq;
	if (nB > 0) ctx.B[nB - 1] += ieq;
}


void Diode::getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx) const 
{
	int nA = nodes[0], nB = nodes[1];
	if (nA > 0) patterns.push_back({nA - 1, nA - 1});
	if (nB > 0) patterns.push_back({nB - 1, nB - 1});
	if (nA > 0 && nB > 0) 
	{
		patterns.push_back({nA - 1, nB - 1});
		patterns.push_back({nB - 1, nA - 1});
	}
}

double Diode::calculateCurrent(const std::vector<double>& state, int extraIdx) const 
{
	return Is * (std::exp(calculateVoltageDrop(state) / Vt) - 1.0);
}

double Diode::calculateVoltageDrop(const std::vector<double>& state) const 
{
	double vA = (nodes[0] > 0) ? state[nodes[0] - 1] : 0.0;
	double vB = (nodes[1] > 0) ? state[nodes[1] - 1] : 0.0;
	return vA - vB;
}

// --- BJT ---
BJT::BJT(const std::string& name, int collector, int base, int emitter, bool isNpn, double Bf, double Br, double Is) 
    : Component(name, ComponentType::BJT, {collector, base, emitter}), isNpn(isNpn), Bf(Bf), Br(Br), Is(Is) {}

void BJT::stamp(MnaMatrixContext& ctx) 
{
	int nC = nodes[0], nB = nodes[1], nE = nodes[2];

	// 1. Extract absolute node voltages from the active Newton-Raphson approximation workspace (Option B)
	double vC = (nC > 0) ? ctx.currentIterV[nC - 1] : 0.0;
	double vB = (nB > 0) ? ctx.currentIterV[nB - 1] : 0.0;
	double vE = (nE > 0) ? ctx.currentIterV[nE - 1] : 0.0;

	// 2. Adjust internal physical junction voltages natively based on structural polarity
	double polarity = isNpn ? 1.0 : -1.0;
	double vbe_target = polarity * (vB - vE);
	double vbc_target = polarity * (vB - vC);
    
	constexpr double VT = 0.026;
	double Vcrit = VT * std::log(VT / (1.41421356 * Is));

	// 3. Apply standard SPICE parallel PN-junction voltage bounding loops
	double vbe_lim = MathSolver::pnjlim(vbe_target, vbe_prev, VT, Vcrit);
	double vbc_lim = MathSolver::pnjlim(vbc_target, vbc_prev, VT, Vcrit);

	// 4. Evaluate Core Ebers-Moll Exponential Elements
	double expVbe = std::exp(vbe_lim / VT);
	double expVbc = std::exp(vbc_lim / VT);

	double icc = Is * (expVbe - 1.0);
	double iec = Is * (expVbc - 1.0);

	// 5. Compute Clean Noise-Free Raw Slopes for Compact Formulation Workspace
	double gbe_raw = (Is / (Bf * VT)) * expVbe;
	double gbc_raw = (Is / (Br * VT)) * expVbc;
	double gcc_raw = (Is / VT) * expVbe;
	double gec_raw = (Is / VT) * expVbc;

	// 6. Compute Physical Device Terminal Operating-Point Currents Natively
	double Ie_op, Ic_op;
	if (isNpn) 
	{
		Ie_op = -icc + iec - (icc / Bf);
		Ic_op = icc - iec - (iec / Br);
	} 
	else 
	{ 
		Ie_op = icc - iec + (icc / Bf);
		Ic_op = -icc + iec - (iec / Br);
	}

	// 7. Compact Branch-Drop Current Formulation (Eliminates Microvolt/Picowatt Floating Point Instabilities)
	double Ieeq = Ie_op - (gbe_raw + gcc_raw) * vbe_lim + gcc_raw * vbc_lim;
	double Iceq = Ic_op + gec_raw * vbe_lim - (gbc_raw + gec_raw) * vbc_lim;
    
	// Invert source vectors if running in a active PNP structural region
	if (!isNpn) 
	{
		Ieeq = -Ieeq;
		Iceq = -Iceq;
	}
	double Ibeq = -(Ieeq + Iceq);

	// 8. Scale Small-Signal Conductances by Polarity to Stabilize Jacobian Traps
	double gbe = gbe_raw * polarity;
	double gbc = gbc_raw * polarity;
	double gcc = gcc_raw * polarity;
	double gec = gec_raw * polarity;

	// 9. Load Polarized Conductance Elements Into CsrMatrix G
	if (nB > 0) 
	{
		MathSolver::accumulateCsr(ctx.G, nB - 1, nB - 1, gbe + gbc);
		if (nC > 0) MathSolver::accumulateCsr(ctx.G, nB - 1, nC - 1, -gbc);
		if (nE > 0) MathSolver::accumulateCsr(ctx.G, nB - 1, nE - 1, -gbe);
	}
	if (nC > 0) 
	{
		MathSolver::accumulateCsr(ctx.G, nC - 1, nC - 1, gbc + gec);
		if (nB > 0) MathSolver::accumulateCsr(ctx.G, nC - 1, nB - 1, gcc - gbc);
		if (nE > 0) MathSolver::accumulateCsr(ctx.G, nC - 1, nE - 1, -gcc - gec);
	}
	if (nE > 0) 
	{
		MathSolver::accumulateCsr(ctx.G, nE - 1, nE - 1, gbe + gcc);
		if (nB > 0) MathSolver::accumulateCsr(ctx.G, nE - 1, nB - 1, gec - gbe);
		if (nC > 0) MathSolver::accumulateCsr(ctx.G, nE - 1, nC - 1, -gec - gcc);
	}

	// 10. Load Balanced Linear RHS Values into Vector B
	if (nC > 0) ctx.B[nC - 1] -= Iceq;
	if (nB > 0) ctx.B[nB - 1] -= Ibeq;
	if (nE > 0) ctx.B[nE - 1] -= Ieeq;
}



void BJT::getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx) const 
{
	int nC = nodes[0], nB = nodes[1], nE = nodes[2];
	std::vector<int> activeNodes;
	if (nC > 0) activeNodes.push_back(nC - 1);
	if (nB > 0) activeNodes.push_back(nB - 1);
	if (nE > 0) activeNodes.push_back(nE - 1);

	for (int r : activeNodes) 
	{
		for (int c : activeNodes) 
		{
			patterns.push_back({r, c});
		}
	}
}

double BJT::calculateCurrent(const std::vector<double>& state, int extraIdx) const 
{
	int nC = nodes[0], nB = nodes[1], nE = nodes[2];
	double vC = (nC > 0) ? state[nC - 1] : 0.0;
	double vB = (nB > 0) ? state[nB - 1] : 0.0;
	double vE = (nE > 0) ? state[nE - 1] : 0.0;
    
	double polarity = isNpn ? 1.0 : -1.0;
	double vbe = polarity * (vB - vE);
	double vbc = polarity * (vB - vC);
    
	double iec = Is * (std::exp(vbc / 0.026) - 1.0);
	double icc = Is * (std::exp(vbe / 0.026) - 1.0);
	double ibr = (Is / Br) * (std::exp(vbc / 0.026) - 1.0);
    
	return polarity * (icc - iec - ibr);
}

double BJT::calculateVoltageDrop(const std::vector<double>& state) const 
{
	int nC = nodes[0], nE = nodes[2];
	double vC = (nC > 0) ? state[nC - 1] : 0.0;
	double vE = (nE > 0) ? state[nE - 1] : 0.0;
	return isNpn ? (vC - vE) : (vE - vC);
}
