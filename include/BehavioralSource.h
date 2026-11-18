// freyalab/BehavioralSource.h
#pragma once

#include "freyalab/CircuitComponents.h"
#include "freyalab/math.h"
#include <complex>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

/**
 * Forward-mode Dual Number for exact automatic Jacobian differentiation:
 *   x = val + deriv * epsilon, where epsilon^2 = 0.
 */
struct Dual {
    double val = 0.0;
    double deriv = 0.0;

    constexpr Dual() = default;
    constexpr Dual(double v, double d = 0.0) : val(v), deriv(d) {}
};

/**
 * Arbitrary Nonlinear Behavioral Source (SPICE B-card):
 *   Bname N+ N- V = <expression>
 *   Bname N+ N- I = <expression>
 *
 * Supports node references V(node), differential voltages V(n+, n-),
 * global .PARAM constants, simulation time TIME, and nonlinear functions:
 *   SIN, COS, TAN, TANH, EXP, LOG, LN, LOG10, SQRT, ABS, POW, MIN, MAX, URAMP, STP
 */
class BehavioralSource : public Component {
public:
    enum class Mode {
        VOLTAGE, // V = f(V_nodes, time) -> allocates 1 branch current row
        CURRENT  // I = f(V_nodes, time) -> stamps directly into KCL node rows
    };

private:
    Mode mode = Mode::CURRENT;
    std::string expression;
    std::map<std::string, double> paramContext;
    std::function<NetNode(const std::string&)> nodeResolver;
    std::vector<NetNode> dependencyNodes; // Unique external NetNodes referenced inside V(...)

    void discoverDependencies();
    Dual evaluateWithSeed(const std::vector<double>* stateV,
                          double timeVal,
                          NetNode seedNode) const;

public:
    BehavioralSource(const std::string& name,
                     const std::vector<NetNode>& terminalNodes,
                     Mode sourceMode,
                     const std::string& expr,
                     const std::map<std::string, double>& params,
                     std::function<NetNode(const std::string&)> resolver);

    int getExtraRowRequirement() const override {
        return (mode == Mode::VOLTAGE) ? 1 : 0;
    }

    void allocateRows(int& globalExtraRowCounter) override;
    void stamp(MnaMatrixContext& ctx) override;
    void getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns,
                            MnaIndex extraIdx = MnaIndex{}) const override;
    void stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                 double omega, const std::vector<double>& dcState) override;

    double calculateVoltageDrop(const std::vector<double>& state) const override;
    double calculateCurrent(const SimulationContext& simCtx,
                            MnaIndex extraIdx = MnaIndex{}) const override;

    double getValue() const override { return 0.0; }
    void setValue(double newValue) override { (void)newValue; }
};

