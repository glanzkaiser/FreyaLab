#pragma once

#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <utility> 

struct CsrMatrix; 

enum class SourceType { DC, AC, PULSE };
enum class ComponentType { RESISTOR, VOLTAGE_SOURCE, CAPACITOR, INDUCTOR, DIODE, BJT };

struct MnaMatrixContext {
    CsrMatrix& G;
    std::vector<double>& B;
    const std::vector<double>& currentIterV;   // Active Newton-Raphson guess (v_k)
    const std::vector<double>& historicStepV;  // Fully converged solution from previous time-step (v_n)
    double time;
    double dt;
    int extraIdx;
};

class Component {
protected:
    std::string name;
    ComponentType type;
    std::vector<int> nodes;

public:
    Component(const std::string& name, ComponentType type, const std::vector<int>& nodes);
    virtual ~Component() = default;

    std::string getName() const { return name; }
    ComponentType getType() const { return type; }
    const std::vector<int>& getNodes() const { return nodes; }

    virtual int getExtraRowRequirement() const { return 0; }
    virtual void stamp(MnaMatrixContext& ctx) = 0;
    virtual void getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx = -1) const = 0;

    virtual double calculateCurrent(const std::vector<double>& state, int extraIdx = -1) const = 0;
    virtual double calculateVoltageDrop(const std::vector<double>& state) const = 0;
    double calculatePower(const std::vector<double>& state, int extraIdx = -1) const {
        return calculateVoltageDrop(state) * calculateCurrent(state, extraIdx);
    }
};

class Resistor : public Component {
private:
    double resistance;
public:
    Resistor(const std::string& name, int nodeA, int nodeB, double resistance);
    void stamp(MnaMatrixContext& ctx) override;
    void getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx = -1) const override;
    double calculateCurrent(const std::vector<double>& state, int extraIdx = -1) const override;
    double calculateVoltageDrop(const std::vector<double>& state) const override;
};

class VoltageSource : public Component {
private:
    SourceType srcType;
    double value1, value2, delay, rise, fall, width, period;
public:
    VoltageSource(const std::string& name, int nodePos, int nodeNeg, double dcVal);
    VoltageSource(const std::string& name, int nodePos, int nodeNeg, SourceType type, 
                  double v1, double v2, double td, double tr=0, double tf=0, double tw=0, double tp=0);
    int getExtraRowRequirement() const override { return 1; }
    void stamp(MnaMatrixContext& ctx) override;
    void getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx = -1) const override;
    double getCurrentValue(double time) const;
    SourceType getSourceType() const { return srcType; }
    void setDCValue(double newDc) { if (srcType == SourceType::DC) value1 = newDc; }
    double calculateCurrent(const std::vector<double>& state, int extraIdx = -1) const override;
    double calculateVoltageDrop(const std::vector<double>& state) const override;
};

class Capacitor : public Component {
private:
    double capacitance;
    mutable double lastVdrop = 0.0; 
public:
    friend class DCCircuit;
    Capacitor(const std::string& name, int nodeA, int nodeB, double capacitance);
    void stamp(MnaMatrixContext& ctx) override;
    void getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx = -1) const override;
    double calculateCurrent(const std::vector<double>& state, int extraIdx = -1) const override;
    double calculateVoltageDrop(const std::vector<double>& state) const override;
};

class Inductor : public Component {
private:
    double inductance;
public:
    friend class DCCircuit;
    Inductor(const std::string& name, int nodeA, int nodeB, double inductance);
    int getExtraRowRequirement() const override { return 1; }
    void stamp(MnaMatrixContext& ctx) override;
    void getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx = -1) const override;
    double calculateCurrent(const std::vector<double>& state, int extraIdx = -1) const override;
    double calculateVoltageDrop(const std::vector<double>& state) const override;
};

class Diode : public Component {
public:
    double Is, Vt;
    mutable double vd_prev_iter = 0.0; 

    friend class DCCircuit;
    Diode(const std::string& name, int anode, int cathode, double Is = 1e-14, double Vt = 0.026);
    void stamp(MnaMatrixContext& ctx) override;
    void getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx = -1) const override;
    double calculateCurrent(const std::vector<double>& state, int extraIdx = -1) const override;
    double calculateVoltageDrop(const std::vector<double>& state) const override;
    void resetState() { vd_prev_iter = 0.0; }
};

class BJT : public Component {
public:
    bool isNpn;
    double Bf, Br, Is;
    mutable double vbe_prev = 0.0; 
    mutable double vbc_prev = 0.0;

    friend class DCCircuit;
    BJT(const std::string& name, int collector, int base, int emitter, bool isNpn = true, double Bf = 100, double Br = 1, double Is = 1e-14);
    void stamp(MnaMatrixContext& ctx) override;
    void getSparsityPattern(std::vector<std::pair<int, int>>& patterns, int extraIdx = -1) const override;
    double calculateCurrent(const std::vector<double>& state, int extraIdx = -1) const override;
    double calculateVoltageDrop(const std::vector<double>& state) const override;
    void resetState() { vbe_prev = 0.0; vbc_prev = 0.0; }
};

