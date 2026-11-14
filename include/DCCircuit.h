// DCCircuit.h
#pragma once
#include "Netlist.h"
#include <string>
#include <vector>

class DCCircuit {
private:
    Netlist& netlist;

    // Aligns, calculates, and prints component parameters into the console matrix table
    void printOperatingPointTable(const std::vector<double>& state);

public:
    explicit DCCircuit(Netlist& netlistInstance);

    // Resolves static network equations by executing the operating point Newton-Raphson loop
    bool solveMna();

    // Runs a static DC Operating Point analysis and outputs a formatted metric table
    void runOperatingPointAnalysis();

    // Sweeps an independent source value, resolves network matrices, and outputs 3 CSV files
    void runDcSweepAnalysis(const std::string& sourceName, double startVal, double stopVal, double stepVal);
};

