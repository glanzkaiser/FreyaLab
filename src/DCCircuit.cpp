#include "freyalab/DCCircuit.h"
#include "freyalab/CircuitComponents.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

DCCircuit::DCCircuit(Netlist& netlistInstance) : netlist(netlistInstance) {}

bool DCCircuit::solveMna() {
    int systemSize = netlist.getTotalSystemSize();
    if (systemSize == 0) {
        netlist.buildSparsityPattern();
        systemSize = netlist.getTotalSystemSize();
    }

    MnaMatrixContext ctx{
        netlist.getGMatrix(),
        netlist.getBVector(),
        &netlist.getStateV(),
        &netlist.getHistoricV(),
        &netlist.getMidpointV(),
        0.0,         // Time
        0.0,         // dt
        0,           // orderStep: 0 = DC Operating Point
        MnaIndex{}   // assignedExtraIdx: Managed dynamically during Component iteration loop
    };

    return netlist.runNewtonRaphson(ctx, 150, 1e-3, 1e-6);
}

void DCCircuit::printOperatingPointTable(const std::vector<double>& state) {
    std::cout << "\n=======================================================================================\n"
              << "                                DC OPERATING POINT TABLE                               \n"
              << "=======================================================================================\n"
              << std::left << std::setw(12) << "Name"
              << std::setw(15) << "Nodes"
              << std::setw(15) << "V-Drop(V)"
              << std::setw(15) << "Current(A)"
              << std::setw(15) << "Power(W)"
              << std::setw(15) << "Type" << "\n"
              << "---------------------------------------------------------------------------------------\n";

    SimulationContext simCtx{ state, 0.0, 0 };

    netlist.inspectComponents([&state, &simCtx](const Component& comp, MnaIndex extraIdx) {
        double vDrop = comp.calculateVoltageDrop(state);
        double current = comp.calculateCurrent(simCtx, extraIdx);

        // Generator convention sign correction for independent Voltage Sources
        if (comp.getType() == ComponentType::VOLTAGE_SOURCE) {
            current = -current;
        }

        double power = std::abs(vDrop * current);

        std::string typeStr;
        switch (comp.getType()) {
            case ComponentType::RESISTOR:       typeStr = "Resistor"; break;
            case ComponentType::VOLTAGE_SOURCE: typeStr = "V-Source"; break;
            case ComponentType::CAPACITOR:      typeStr = "Capacitor"; break;
            case ComponentType::INDUCTOR:       typeStr = "Inductor"; break;
            case ComponentType::DIODE:          typeStr = "Diode"; break;
            case ComponentType::MOSFET:         typeStr = "MOSFET"; break;
            case ComponentType::FINFET:         typeStr = "FinFET"; break;
            case ComponentType::RIBBONFET:      typeStr = "RibbonFET"; break;
            case ComponentType::GAA_FET:        typeStr = "GAA-FET"; break;
            case ComponentType::BJT:            typeStr = "BJT"; break;
            case ComponentType::DIGITAL_IC:     typeStr = "Logic Gate"; break;
            default:                            typeStr = "Behavioral"; break;
        }

        std::stringstream ss;
        const auto& nodes = comp.getNodes();
        for (size_t i = 0; i < nodes.size(); ++i) {
            ss << nodes[i].id;
            if (i + 1 < nodes.size()) ss << ",";
        }

        std::cout << std::left << std::setw(12) << comp.getName()
                  << std::setw(15) << ss.str()
                  << std::setw(15) << std::fixed << std::setprecision(5) << vDrop
                  << std::setw(15) << current
                  << std::setw(15) << power
                  << std::setw(15) << typeStr << "\n";
    });

    std::cout << "=======================================================================================\n";
    std::cout << "                       DETAILED MULTI-TERMINAL ACTIVE DEVICE TELEMETRY                 \n";
    std::cout << "=======================================================================================\n";
    std::cout << std::left << std::setw(12) << "Device"
              << std::setw(14) << "Model"
              << std::setw(20) << "Term 1 (C/D)"
              << std::setw(20) << "Term 2 (B/G)"
              << std::setw(20) << "Term 3 (E/S)" << "\n";
    std::cout << "---------------------------------------------------------------------------------------\n";

    netlist.inspectComponents([&state, &simCtx](const Component& comp, MnaIndex extraIdx) {
        const auto& nodes = comp.getNodes();
        if (nodes.size() < 3) return;

        ComponentType cType = comp.getType();

        if (cType == ComponentType::BJT) {
            double vC = getSafeVoltage(state, nodes[0]);
            double vB = getSafeVoltage(state, nodes[1]);
            double vE = getSafeVoltage(state, nodes[2]);

            // Use the device's own polarized VBIC current model + base current approximation
            double I_c = comp.calculateCurrent(simCtx, extraIdx);
            double Vt = 0.025851991;
            double vbe_abs = std::abs(vB - vE);
            double vbc_abs = std::abs(vB - vC);
            double If = 1e-16 * (std::exp(std::min(vbe_abs / Vt, 40.0)) - 1.0);
            double Ir = 1e-16 * (std::exp(std::min(vbc_abs / Vt, 40.0)) - 1.0);
            double sign = (I_c >= 0.0) ? 1.0 : -1.0;
            double I_b = sign * ((If / 100.0) + (Ir / 50.0));
            double I_e = -(I_c + I_b);

            std::cout << std::left << std::setw(12) << comp.getName()
                      << std::setw(14) << "VBIC-BJT"
                      << "Ic: " << std::setw(16) << std::scientific << std::setprecision(4) << I_c
                      << "Ib: " << std::setw(16) << I_b
                      << "Ie: " << std::setw(16) << I_e << "\n";
        }
        else if (cType == ComponentType::MOSFET ||
                 cType == ComponentType::FINFET ||
                 cType == ComponentType::RIBBONFET ||
                 cType == ComponentType::GAA_FET) {
            std::string modelLabel = "MOSFET";
            if (cType == ComponentType::FINFET)         modelLabel = "FinFET";
            else if (cType == ComponentType::RIBBONFET) modelLabel = "RibbonFET";
            else if (cType == ComponentType::GAA_FET)   modelLabel = "GAA-FET";

            // Delegates directly to Mosfet::calculateCurrent or BsimCmgFet::calculateCurrent
            double I_d = comp.calculateCurrent(simCtx, extraIdx);
            double I_g = 0.0;
            double I_s = -I_d;

            std::cout << std::left << std::setw(12) << comp.getName()
                      << std::setw(14) << modelLabel
                      << "Id: " << std::setw(16) << std::scientific << std::setprecision(4) << I_d
                      << "Ig: " << std::setw(16) << I_g
                      << "Is: " << std::setw(16) << I_s << "\n";
        }
    });
    std::cout << "=======================================================================================\n\n";
}

void DCCircuit::runOperatingPointAnalysis() {
    std::cout << "Starting DC Operating Point Analysis..." << std::endl;
    if (solveMna()) {
        const auto& state = netlist.getCurrentState();
        printOperatingPointTable(state);
    } else {
        std::cerr << "Error: DC Operating Point failed to find a stable solution." << std::endl;
    }
}

void DCCircuit::runDcSweepAnalysis(const std::string& sourceName, double startVal, double stopVal, double stepVal) {
    std::cout << "Executing DC Parametric Sweep Analysis on source: " << sourceName << std::endl;

    std::ofstream csvV("dc_sweep_voltage.csv");
    std::ofstream csvI("dc_sweep_current.csv");
    std::ofstream csvP("dc_sweep_power.csv");

    if (!csvV.is_open() || !csvI.is_open() || !csvP.is_open()) {
        throw std::runtime_error("Fatal Error: Could not generate DC Sweep CSV file handles.");
    }

    csvV << "SweepValue";
    csvI << "SweepValue";
    csvP << "SweepValue";

    int systemSize = netlist.getTotalSystemSize();
    if (systemSize == 0) {
        netlist.buildSparsityPattern();
        systemSize = netlist.getTotalSystemSize();
    }

    netlist.inspectComponents([&csvV, &csvI, &csvP](const Component& comp, MnaIndex extraIdx) {
        (void)extraIdx;
        csvV << "," << comp.getName() << "_VDrop";
        csvI << "," << comp.getName() << "_Current";
        csvP << "," << comp.getName() << "_Power";
    });
    csvV << "\n";
    csvI << "\n";
    csvP << "\n";

    if (stepVal <= 0.0) stepVal = 0.1;
    double currentSweepVal = startVal;

    std::vector<double> lastValidState = netlist.getCurrentState();

    while (currentSweepVal <= stopVal + (stepVal * 0.01)) {
        // Mutate target DC source (or BSIM-CMG VTH0) sweep value
        netlist.inspectComponents([&sourceName, currentSweepVal](Component& comp, MnaIndex extraIdx) {
            (void)extraIdx;
            if (comp.getName() == sourceName) {
                comp.setValue(currentSweepVal);
            }
        });

        MnaMatrixContext ctx{
            netlist.getGMatrix(),
            netlist.getBVector(),
            &netlist.getStateV(),
            &netlist.getHistoricV(),
            &netlist.getMidpointV(),
            0.0,         // time
            0.0,         // dt (0.0 for static DC)
            0,           // orderStep: 0 = DC Operating Point
            MnaIndex{}   // assignedExtraIdx: managed per-component inside Netlist
        };

        try {
            // Solver continuity: use the last converged solution as initial guess
            if (!lastValidState.empty() && lastValidState.size() == netlist.getStateV().size()) {
                netlist.getStateV() = lastValidState;
            }

            if (!netlist.runNewtonRaphson(ctx, 100, 1e-3, 1e-6)) {
                throw std::runtime_error("Newton-Raphson divergence on sweep element step.");
            }

            const auto& state = netlist.getCurrentState();
            lastValidState = state;

            csvV << currentSweepVal;
            csvI << currentSweepVal;
            csvP << currentSweepVal;

            SimulationContext sweepSimCtx{ state, 0.0, 0 };

            // Collect telemetry across all components (MOSFET, BJT, FinFET, RibbonFET, GAA-FET)
            netlist.inspectComponents([&state, &sweepSimCtx, &csvV, &csvI, &csvP](const Component& comp, MnaIndex extraIdx) {
                double vDrop = comp.calculateVoltageDrop(state);
                double current = comp.calculateCurrent(sweepSimCtx, extraIdx);

                if (comp.getType() == ComponentType::VOLTAGE_SOURCE) {
                    current = -current;
                }

                double power = std::abs(vDrop * current);

                csvV << "," << vDrop;
                csvI << "," << current;
                csvP << "," << power;
            });

            csvV << "\n";
            csvI << "\n";
            csvP << "\n";

        } catch (const std::exception& e) {
            std::cerr << "Warning: Sweep step failure at value " << currentSweepVal
                      << " due to convergence breakdown: " << e.what() << std::endl;
        }

        currentSweepVal += stepVal;
    }

    csvV.close();
    csvI.close();
    csvP.close();

    std::cout << "DC Sweep completed successfully. Matrix datasets committed to CSV." << std::endl;
}
