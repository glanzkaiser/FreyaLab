// BehavioralSource.cpp
#include "freyalab/BehavioralSource.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

// --- Dual-Number Elementary Arithmetic & Nonlinear Chain-Rule Kernels ---
inline Dual operator+(const Dual& a, const Dual& b) { return {a.val + b.val, a.deriv + b.deriv}; }
inline Dual operator-(const Dual& a, const Dual& b) { return {a.val - b.val, a.deriv - b.deriv}; }
inline Dual operator*(const Dual& a, const Dual& b) { return {a.val * b.val, a.deriv * b.val + a.val * b.deriv}; }

inline Dual operator/(const Dual& a, const Dual& b) {
    double safeB = (std::abs(b.val) > 1e-30) ? b.val : ((b.val < 0.0) ? -1e-30 : 1e-30);
    return {a.val / safeB, (a.deriv * safeB - a.val * b.deriv) / (safeB * safeB)};
}

inline Dual dualPow(const Dual& a, const Dual& b) {
    if (std::abs(b.deriv) < 1e-30) {
        double v = std::pow(a.val, b.val);
        double d = (b.val == 0.0 || a.val == 0.0) ? 0.0 : b.val * std::pow(a.val, b.val - 1.0) * a.deriv;
        return {v, d};
    }
    double safeBase = std::max(a.val, 1e-30);
    double v = std::pow(safeBase, b.val);
    double d = v * (b.deriv * std::log(safeBase) + (b.val * a.deriv) / safeBase);
    return {v, d};
}

inline Dual dualSin(const Dual& x)  { return {std::sin(x.val),  std::cos(x.val) * x.deriv}; }
inline Dual dualCos(const Dual& x)  { return {std::cos(x.val), -std::sin(x.val) * x.deriv}; }
inline Dual dualTan(const Dual& x)  {
    double c = std::cos(x.val);
    double sec2 = 1.0 / std::max(c * c, 1e-18);
    return {std::tan(x.val), sec2 * x.deriv};
}
inline Dual dualTanh(const Dual& x) {
    double th = std::tanh(x.val);
    return {th, (1.0 - th * th) * x.deriv};
}
inline Dual dualExp(const Dual& x)  {
    double ex = std::exp(std::clamp(x.val, -80.0, 80.0));
    return {ex, ex * x.deriv};
}
inline Dual dualLog(const Dual& x)  {
    double safeX = std::max(x.val, 1e-30);
    return {std::log(safeX), x.deriv / safeX};
}
inline Dual dualLog10(const Dual& x) {
    double safeX = std::max(x.val, 1e-30);
    constexpr double INV_LN10 = 0.4342944819032518;
    return {std::log10(safeX), (INV_LN10 * x.deriv) / safeX};
}
inline Dual dualSqrt(const Dual& x) {
    double safeX = std::max(x.val, 1e-30);
    double sq = std::sqrt(safeX);
    return {sq, (0.5 * x.deriv) / sq};
}
inline Dual dualAbs(const Dual& x)  {
    double sgn = (x.val > 0.0) ? 1.0 : ((x.val < 0.0) ? -1.0 : 0.0);
    return {std::abs(x.val), sgn * x.deriv};
}
inline Dual dualUramp(const Dual& x) {
    return (x.val > 0.0) ? x : Dual{0.0, 0.0};
}
inline Dual dualStp(const Dual& x) {
    return {(x.val > 0.0) ? 1.0 : 0.0, 0.0};
}

double parseEngineeringLiteral(const std::string& tok) {
    size_t idx = 0;
    double num = std::stod(tok, &idx);
    if (idx >= tok.size()) return num;

    std::string suffix = tok.substr(idx);
    std::transform(suffix.begin(), suffix.end(), suffix.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (suffix == "f")   return num * 1e-15;
    if (suffix == "p")   return num * 1e-12;
    if (suffix == "n")   return num * 1e-9;
    if (suffix == "u")   return num * 1e-6;
    if (suffix == "m")   return num * 1e-3;
    if (suffix == "k")   return num * 1e3;
    if (suffix == "meg") return num * 1e6;
    if (suffix == "g")   return num * 1e9;
    if (suffix == "t")   return num * 1e12;
    return num;
}

class DualBehavioralParser {
public:
    DualBehavioralParser(const std::string& expr,
                         const std::map<std::string, double>& params,
                         const std::function<NetNode(const std::string&)>& resolver,
                         const std::vector<double>* stateV,
                         double timeVal,
                         NetNode seedNode,
                         std::set<NetNode>* collectedDeps)
        : s(expr), pos(0), paramMap(params), resolveNode(resolver),
          state(stateV), t(timeVal), seed(seedNode), deps(collectedDeps) {}

    Dual parse() {
        return parseAddSub();
    }

private:
    std::string s;
    size_t pos;
    const std::map<std::string, double>& paramMap;
    const std::function<NetNode(const std::string&)>& resolveNode;
    const std::vector<double>* state;
    double t;
    NetNode seed;
    std::set<NetNode>* deps;

    void skipWs() {
        while (pos < s.size() && (std::isspace(static_cast<unsigned char>(s[pos])) ||
                                  s[pos] == '{' || s[pos] == '}' || s[pos] == '\'')) {
            ++pos;
        }
    }

    Dual sampleNodeDual(const std::string& rawNodeName) {
        std::string cleanName = rawNodeName;
        cleanName.erase(std::remove_if(cleanName.begin(), cleanName.end(), ::isspace), cleanName.end());
        NetNode node = resolveNode(cleanName);

        if (deps && !node.isGround()) {
            deps->insert(node);
        }

        double v = (state != nullptr) ? getSafeVoltage(state, node) : 0.0;
        double d = (!node.isGround() && node == seed) ? 1.0 : 0.0;
        return Dual{v, d};
    }

    Dual parseAddSub() {
        Dual left = parseMulDiv();
        while (true) {
            skipWs();
            if (pos < s.size() && s[pos] == '+') { ++pos; left = left + parseMulDiv(); }
            else if (pos < s.size() && s[pos] == '-') { ++pos; left = left - parseMulDiv(); }
            else break;
        }
        return left;
    }

    Dual parseMulDiv() {
        Dual left = parsePower();
        while (true) {
            skipWs();
            if (pos < s.size() && s[pos] == '*') { ++pos; left = left * parsePower(); }
            else if (pos < s.size() && s[pos] == '/') { ++pos; left = left / parsePower(); }
            else break;
        }
        return left;
    }

    Dual parsePower() {
        Dual base = parseUnary();
        skipWs();
        if (pos + 1 < s.size() && s[pos] == '*' && s[pos + 1] == '*') {
            pos += 2;
            return dualPow(base, parsePower());
        }
        if (pos < s.size() && s[pos] == '^') {
            ++pos;
            return dualPow(base, parsePower());
        }
        return base;
    }

    Dual parseUnary() {
        skipWs();
        if (pos < s.size() && s[pos] == '+') { ++pos; return parseUnary(); }
        if (pos < s.size() && s[pos] == '-') {
            ++pos;
            Dual u = parseUnary();
            return Dual{-u.val, -u.deriv};
        }
        return parsePrimary();
    }

    Dual parsePrimary() {
        skipWs();
        if (pos >= s.size()) return Dual{0.0, 0.0};

        if (s[pos] == '(') {
            ++pos;
            Dual inner = parseAddSub();
            skipWs();
            if (pos < s.size() && s[pos] == ')') ++pos;
            return inner;
        }

        if (std::isdigit(static_cast<unsigned char>(s[pos])) || s[pos] == '.') {
            size_t start = pos;
            while (pos < s.size() && (std::isalnum(static_cast<unsigned char>(s[pos])) || s[pos] == '.')) {
                if ((s[pos] == 'e' || s[pos] == 'E') && pos + 1 < s.size() &&
                    (s[pos + 1] == '+' || s[pos + 1] == '-')) {
                    pos += 2;
                    continue;
                }
                ++pos;
            }
            return Dual{parseEngineeringLiteral(s.substr(start, pos - start)), 0.0};
        }

        if (std::isalpha(static_cast<unsigned char>(s[pos])) || s[pos] == '_') {
            size_t start = pos;
            while (pos < s.size() && (std::isalnum(static_cast<unsigned char>(s[pos])) || s[pos] == '_')) {
                ++pos;
            }
            std::string ident = s.substr(start, pos - start);
            std::transform(ident.begin(), ident.end(), ident.begin(),
                           [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

            skipWs();
            // Special SPICE nodal voltage reference: V(node) or V(nodePos, nodeNeg)
            if (ident == "V" && pos < s.size() && s[pos] == '(') {
                ++pos;
                size_t closeP = s.find(')', pos);
                std::string inside = (closeP != std::string::npos) ? s.substr(pos, closeP - pos) : "";
                pos = (closeP != std::string::npos) ? (closeP + 1) : s.size();

                size_t comma = inside.find(',');
                if (comma != std::string::npos) {
                    Dual vPos = sampleNodeDual(inside.substr(0, comma));
                    Dual vNeg = sampleNodeDual(inside.substr(comma + 1));
                    return vPos - vNeg;
                }
                return sampleNodeDual(inside);
            }

            // Standard 1-arg or 2-arg mathematical functions
            if (pos < s.size() && s[pos] == '(') {
                ++pos;
                Dual arg1 = parseAddSub();
                Dual arg2{0.0, 0.0};
                bool hasArg2 = false;
                skipWs();
                if (pos < s.size() && s[pos] == ',') {
                    ++pos;
                    arg2 = parseAddSub();
                    hasArg2 = true;
                }
                skipWs();
                if (pos < s.size() && s[pos] == ')') ++pos;

                                if (ident == "SIN")   return dualSin(arg1);
                if (ident == "COS")   return dualCos(arg1);
                if (ident == "TAN")   return dualTan(arg1);
                if (ident == "TANH")  return dualTanh(arg1);
                if (ident == "EXP")   return dualExp(arg1);
                if (ident == "LN" || ident == "LOG") return dualLog(arg1);
                if (ident == "LOG10") return dualLog10(arg1);
                if (ident == "SQRT")  return dualSqrt(arg1);
                if (ident == "ABS")   return dualAbs(arg1);
                if (ident == "URAMP") return dualUramp(arg1);
                if (ident == "STP")   return dualStp(arg1);
                if (ident == "POW" && hasArg2) return dualPow(arg1, arg2);
                if (ident == "MIN" && hasArg2) return (arg1.val <= arg2.val) ? arg1 : arg2;
                if (ident == "MAX" && hasArg2) return (arg1.val >= arg2.val) ? arg1 : arg2;
                return arg1;
            }

            if (ident == "TIME") {
                return Dual{t, 0.0};
            }

            auto it = paramMap.find(ident);
            if (it != paramMap.end()) {
                return Dual{it->second, 0.0};
            }
            return Dual{0.0, 0.0};
        }

        return Dual{0.0, 0.0};
    }
};

} // namespace

// ============================================================================
// BehavioralSource Member Implementations
// ============================================================================

BehavioralSource::BehavioralSource(const std::string& name,
                                   const std::vector<NetNode>& terminalNodes,
                                   Mode sourceMode,
                                   const std::string& expr,
                                   const std::map<std::string, double>& params,
                                   std::function<NetNode(const std::string&)> resolver)
    : Component(name, ComponentType::BEHAVIORAL_SOURCE, terminalNodes),
      mode(sourceMode),
      expression(expr),
      paramContext(params),
      nodeResolver(std::move(resolver)) {
    discoverDependencies();
}

void BehavioralSource::discoverDependencies() {
    std::set<NetNode> uniqueDeps;
    DualBehavioralParser scanner(expression, paramContext, nodeResolver, nullptr, 0.0, NetNode{0}, &uniqueDeps);
    scanner.parse();

    dependencyNodes.clear();
    for (NetNode n : uniqueDeps) {
        if (!n.isGround()) {
            dependencyNodes.push_back(n);
        }
    }
}

Dual BehavioralSource::evaluateWithSeed(const std::vector<double>* stateV,
                                        double timeVal,
                                        NetNode seedNode) const {
    DualBehavioralParser parser(expression, paramContext, nodeResolver, stateV, timeVal, seedNode, nullptr);
    return parser.parse();
}

void BehavioralSource::allocateRows(int& globalExtraRowCounter) {
    if (mode == Mode::VOLTAGE) {
        assignedExtraRowIdx = MnaIndex::fromExtraRow(globalExtraRowCounter++);
    }
}

void BehavioralSource::getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns,
                                          MnaIndex extraIdx) const {
    if (nodes.size() < 2) return;

    MnaIndex mPos = MnaIndex::fromNetNode(nodes[0]);
    MnaIndex mNeg = MnaIndex::fromNetNode(nodes[1]);

    if (mode == Mode::CURRENT) {
        // Current source I = f(V_dep) injects Jacobian entries dI/dV_dep into N+ and N- rows
        for (NetNode dep : dependencyNodes) {
            MnaIndex mDep = MnaIndex::fromNetNode(dep);
            if (!mDep.isValid()) continue;
            if (mPos.isValid()) patterns.push_back({mPos, mDep});
            if (mNeg.isValid()) patterns.push_back({mNeg, mDep});
        }
    } else {
        // Voltage source V = f(V_dep) allocates branch row `idx`
        MnaIndex idx = extraIdx.isValid() ? extraIdx : assignedExtraRowIdx;
        if (!idx.isValid()) return;

        if (mPos.isValid()) {
            patterns.push_back({mPos, idx});
            patterns.push_back({idx, mPos});
        }
        if (mNeg.isValid()) {
            patterns.push_back({mNeg, idx});
            patterns.push_back({idx, mNeg});
        }
        patterns.push_back({idx, idx});

        for (NetNode dep : dependencyNodes) {
            MnaIndex mDep = MnaIndex::fromNetNode(dep);
            if (mDep.isValid()) {
                patterns.push_back({idx, mDep});
            }
        }
    }
}

void BehavioralSource::stamp(MnaMatrixContext& ctx) {
    if (nodes.size() < 2) return;

    MnaIndex mPos = MnaIndex::fromNetNode(nodes[0]);
    MnaIndex mNeg = MnaIndex::fromNetNode(nodes[1]);

    // Evaluate f(V_k, t) and exact partial derivatives df/dV_j via dual-number seeding
    double fVal = evaluateWithSeed(ctx.currentIterV, ctx.time, NetNode{0}).val;
    double linearSum = 0.0;

    std::vector<std::pair<MnaIndex, double>> jacobians;
    jacobians.reserve(dependencyNodes.size());

    for (NetNode dep : dependencyNodes) {
        MnaIndex mDep = MnaIndex::fromNetNode(dep);
        if (!mDep.isValid()) continue;

        double g_j = evaluateWithSeed(ctx.currentIterV, ctx.time, dep).deriv;
        double v_j = ctx.currentIterV ? getSafeVoltage(ctx.currentIterV, dep) : 0.0;

        jacobians.push_back({mDep, g_j});
        linearSum += g_j * v_j;
    }

    double rhsOffset = fVal - linearSum;

    if (mode == Mode::CURRENT) {
        for (const auto& kv : jacobians) {
            if (mPos.isValid()) ctx.G.addValue(mPos.idx, kv.first.idx,  kv.second);
            if (mNeg.isValid()) ctx.G.addValue(mNeg.idx, kv.first.idx, -kv.second);
        }
        addSafeRhs(ctx.B, mPos, -rhsOffset);
        addSafeRhs(ctx.B, mNeg,  rhsOffset);
    } else {
        MnaIndex idx = ctx.assignedExtraIdx.isValid() ? ctx.assignedExtraIdx : assignedExtraRowIdx;
        if (!idx.isValid() || static_cast<size_t>(idx.idx) >= ctx.B.size()) return;

        if (mPos.isValid()) {
            ctx.G.addValue(mPos.idx, idx.idx,  1.0);
            ctx.G.addValue(idx.idx, mPos.idx,  1.0);
        }
        if (mNeg.isValid()) {
            ctx.G.addValue(mNeg.idx, idx.idx, -1.0);
            ctx.G.addValue(idx.idx, mNeg.idx, -1.0);
        }
        ctx.G.addValue(idx.idx, idx.idx, 0.0);

        for (const auto& kv : jacobians) {
            ctx.G.addValue(idx.idx, kv.first.idx, -kv.second);
        }
        setSafeRhs(ctx.B, idx, rhsOffset);
    }
}

void BehavioralSource::stampAC(CsrMatrixComplex& A_c,
                               std::vector<std::complex<double>>& B_c,
                               double omega,
                               const std::vector<double>& dcState) {
    (void)B_c;
    (void)omega;
    if (nodes.size() < 2) return;

    MnaIndex mPos = MnaIndex::fromNetNode(nodes[0]);
    MnaIndex mNeg = MnaIndex::fromNetNode(nodes[1]);

    if (mode == Mode::CURRENT) {
        for (NetNode dep : dependencyNodes) {
            MnaIndex mDep = MnaIndex::fromNetNode(dep);
            if (!mDep.isValid()) continue;
            double g_j = evaluateWithSeed(&dcState, 0.0, dep).deriv;
            if (mPos.isValid()) A_c.addValue(mPos.idx, mDep.idx, std::complex<double>( g_j, 0.0));
            if (mNeg.isValid()) A_c.addValue(mNeg.idx, mDep.idx, std::complex<double>(-g_j, 0.0));
        }
    } else {
        MnaIndex idx = assignedExtraRowIdx;
        if (!idx.isValid()) return;

        if (mPos.isValid()) {
            A_c.addValue(mPos.idx, idx.idx, std::complex<double>(1.0, 0.0));
            A_c.addValue(idx.idx, mPos.idx, std::complex<double>(1.0, 0.0));
        }
        if (mNeg.isValid()) {
            A_c.addValue(mNeg.idx, idx.idx, std::complex<double>(-1.0, 0.0));
            A_c.addValue(idx.idx, mNeg.idx, std::complex<double>(-1.0, 0.0));
        }

        for (const auto& dep : dependencyNodes) {
            MnaIndex mDep = MnaIndex::fromNetNode(dep);
            if (!mDep.isValid()) continue;
            double g_j = evaluateWithSeed(&dcState, 0.0, dep).deriv;
            A_c.addValue(idx.idx, mDep.idx, std::complex<double>(-g_j, 0.0));
        }
    }
}

double BehavioralSource::calculateVoltageDrop(const std::vector<double>& state) const {
    if (nodes.size() < 2) return 0.0;
    return getSafeVoltage(state, nodes[0]) - getSafeVoltage(state, nodes[1]);
}

double BehavioralSource::calculateCurrent(const SimulationContext& simCtx,
                                          MnaIndex extraIdx) const {
    if (mode == Mode::CURRENT) {
        return evaluateWithSeed(&simCtx.state, simCtx.time, NetNode{0}).val;
    }
    MnaIndex idx = extraIdx.isValid() ? extraIdx : assignedExtraRowIdx;
    return getSafeExtraValue(simCtx.state, idx);
}

