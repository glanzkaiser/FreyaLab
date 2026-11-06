#include "freyalab/DCCircuit.h"
#include "freyalab/CircuitComponents.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <stdexcept>

DCCircuit::DCCircuit(Netlist& netlistInstance) : netlist(netlistInstance) {}

bool DCCircuit::solveMna() {
    int systemSize = netlist.getTotalSystemSize();
    if (systemSize == 0) {
        netlist.buildSparsityPattern();
        systemSize = netlist.getTotalSystemSize();
    }
    
    // Explicit 1-based mapping topology separation
    int nodeCount = netlist.getNumNodes(); 
    int currentExtraCount = systemSize - nodeCount; 
    
    MnaMatrixContext ctx{ 
        netlist.getGMatrix(), 
        netlist.getBVector(), 
        &netlist.getStateV(), 
        &netlist.getHistoricV(), 
        &netlist.getMidpointV(), 
        0.0,    // Time
        0.0,   // dt
        0,      // orderStep: 0 = DC Operating Point
        -1      // assignedExtraIdx: Managed dynamically during Component iteration loop
    };
    
    return netlist.runNewtonRaphson(ctx, 150, 1e-3, 1e-6);
}

void DCCircuit::printOperatingPointTable(const std::vector<double>& state) {
    std::cout << "\n========================================================================\n"
              << "                      DC OPERATING POINT TABLE                          \n"
              << "========================================================================\n"
              << std::left << std::setw(12) << "Name" 
              << std::setw(15) << "Nodes" 
              << std::setw(15) << "V-Drop(V)" 
              << std::setw(15) << "Current(A)" 
              << std::setw(15) << "Power(W)"
	          << std::setw(15) << "Type" << "\n"
              << "------------------------------------------------------------------------\n";

    SimulationContext simCtx{ state, 0.0, 0 };

    netlist.inspectComponents([&state, &simCtx](const Component& comp, int extraIdx) {
        double vDrop = comp.calculateVoltageDrop(state);
        
        // PERBAIKAN SELESTIAL: Konversikan extraIdx global ke 1-based index 
        // agar selaras dengan fungsi calculateCurrent internal komponen
        //int mappedExtraIdx = (extraIdx != -1) ? (extraIdx + 1) : -1;
        double current = comp.calculateCurrent(simCtx, extraIdx);
        
        // Koreksi arah tanda baca konvensi generator arus khusus Voltage Source
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
            case ComponentType::BJT:            typeStr = "BJT"; break;
            case ComponentType::DIGITAL_IC:     typeStr = "Logic Gate"; break;
            default:                            typeStr = "Behavioral"; break;
        }

        std::stringstream ss;
        const auto& nodes = comp.getNodes();
        for (size_t i = 0; i < nodes.size(); ++i) {
            ss << nodes[i];
            if (i < nodes.size() - 1) ss << ",";
        }

        std::cout << std::left << std::setw(12) << comp.getName()
                  << std::setw(15) << ss.str()
                  << std::setw(15) << std::fixed << std::setprecision(5) << vDrop
                  << std::setw(15) << current
                  << std::setw(15) << power 
                  << std::setw(15) << typeStr << "\n";
    });
    
    std::cout << "======================================================================\n";
    std::cout << "             DETAILED MULTI-TERMINAL ACTIVE DEVICE TELEMETRY          \n";
    std::cout << "======================================================================\n";
    std::cout << std::left << std::setw(12) << "Device" 
              << std::setw(18) << "Term 1 (C/D)" 
              << std::setw(18) << "Term 2 (B/G)" 
              << std::setw(18) << "Term 3 (E/S)" << "\n";
    std::cout << "----------------------------------------------------------------------\n";

    netlist.inspectComponents([&state](const Component& comp, int extraIdx) {
        (void)extraIdx;
        const auto& nodes = comp.getNodes();
        int stateSize = static_cast<int>(state.size());
        
        auto getSafeNodeVoltage = [&](int nodeNum) -> double {
            if (nodeNum <= 0) return 0.0;
            int matrixIdx = nodeNum - 1;
            if (matrixIdx < stateSize) return state[matrixIdx];
            return 0.0;
        };

        if (comp.getType() == ComponentType::BJT && nodes.size() >= 3) {
            double vC = getSafeNodeVoltage(nodes[0]);
            double vB = getSafeNodeVoltage(nodes[1]);
            double vE = getSafeNodeVoltage(nodes[2]);
            
            double Vt = 0.02585;
            double If = 1e-16 * (std::exp((vB - vE) / Vt) - 1.0); 
            double Ir = 1e-16 * (std::exp((vB - vC) / Vt) - 1.0);
            double Ice = (If - Ir) / (1.0 + (If / 0.01)); 
            
            double I_c = Ice - Ir; 
            double I_b = If + Ir; 
            double I_e = -(I_c + I_b); 

            std::cout << std::left << std::setw(12) << comp.getName()
                      << "Ic: " << std::setw(14) << std::scientific << std::setprecision(4) << I_c
                      << "Ib: " << std::setw(14) << I_b
                      << "Ie: " << std::setw(14) << I_e << "\n";
        }
        else if (comp.getType() == ComponentType::MOSFET && nodes.size() >= 3) {
            double vD = getSafeNodeVoltage(nodes[0]);
            double vG = getSafeNodeVoltage(nodes[1]);
            double vS = getSafeNodeVoltage(nodes[2]);

            double vgs = vG - vS; double vds = vD - vS;
            double I_d = 0.0;
            if (vgs > 0.7) {
                if (vds < (vgs - 0.7)) I_d = 1e-3 * ((vgs - 0.7) * vds - 0.5 * vds * vds);
                else I_d = 0.5 * 1e-3 * (vgs - 0.7) * (vgs - 0.7);
            }
            double I_g = 0.0;        
            double I_s = -I_d;

            std::cout << std::left << std::setw(12) << comp.getName()
                      << "Id: " << std::setw(14) << std::scientific << std::setprecision(4) << I_d
                      << "Ig: " << std::setw(14) << I_g
                      << "Is: " << std::setw(14) << I_s << "\n";
        }
    });
    std::cout << "======================================================================\n\n";
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

    // Tulis kolom penanda pertama untuk nilai sapuan parameter
    csvV << "SweepValue"; csvI << "SweepValue"; csvP << "SweepValue";

    int systemSize = netlist.getTotalSystemSize();
    if (systemSize == 0) {
        netlist.buildSparsityPattern();
        systemSize = netlist.getTotalSystemSize();
    }

    // Tulis Header kolom CSV secara horizontal
    netlist.inspectComponents([&csvV, &csvI, &csvP](const Component& comp, int extraIdx) {
        (void)extraIdx;
        csvV << "," << comp.getName() << "_VDrop";
        csvI << "," << comp.getName() << "_Current";
        csvP << "," << comp.getName() << "_Power";
    });
    csvV << "\n"; csvI << "\n"; csvP << "\n";

    if (stepVal <= 0.0) stepVal = 0.1;
    double currentSweepVal = startVal;

    std::vector<double> lastValidState = netlist.getCurrentState();

    // Toleransi komparasi floating point aman mencegah infinite loop
    while (currentSweepVal <= stopVal + (stepVal * 0.01)) {
        
        // Mutasi nilai sumber tegangan/arus DC target sweep menggunakan pemanggilan non-const
        netlist.inspectComponents([&sourceName, currentSweepVal](Component& comp, int extraIdx) {
            (void)extraIdx;
            if (comp.getName() == sourceName) {
                comp.setValue(currentSweepVal);
            }
        });

        // PERBAIKAN: Dapatkan jumlah node murni yang tepat untuk isolasi topologi ekstra
        int nodeCount = netlist.getNumNodes();
        int currentExtraCount = systemSize - nodeCount;
        
        // PERBAIKAN: Hapus operator '&' dan gunakan netlist.getBVector() untuk parameter ke-2!
        MnaMatrixContext ctx{ 
            const_cast<CsrMatrix&>(netlist.getGMatrix()), 
            const_cast<std::vector<double>&>(netlist.getHistoricV()), 
            &netlist.getStateV(), 
            &netlist.getHistoricV(), 
            &netlist.getMidpointV(), 
            0.0, 1e-12, 0, currentExtraCount 
        };
        try {
            // Kontinuitas Solver: Gunakan kondisi konvergen terakhir sebagai tebakan awal langkah baru
            if (!lastValidState.empty() && lastValidState.size() == netlist.getStateV().size()) {
                netlist.getStateV() = lastValidState;
            }

            if (!netlist.runNewtonRaphson(ctx, 100, 1e-3, 1e-6)) {
                throw std::runtime_error("Newton-Raphson divergence on sweep element step.");
            }
            
            const auto& state = netlist.getCurrentState();
            lastValidState = state; // Perbarui kondisi valid untuk iterasi selanjutnya

            csvV << currentSweepVal;
            csvI << currentSweepVal;
            csvP << currentSweepVal;

            // PERBAIKAN: Gunakan inisialisasi agregat standar guna mencegah C++20 standard mismatch error
            SimulationContext sweepSimCtx{ state, 0.0, 0 };

            // Kumpulkan hasil kalkulasi final titik penyapuan ke masing-masing file CSV
            netlist.inspectComponents([&state, &sweepSimCtx, &csvV, &csvI, &csvP](const Component& comp, int extraIdx) {
                double vDrop = comp.calculateVoltageDrop(state);
                double current = comp.calculateCurrent(sweepSimCtx, extraIdx);
                double power = std::abs(vDrop * current);

                csvV << "," << vDrop;
                csvI << "," << current;
                csvP << "," << power;
            });
            
            csvV << "\n"; csvI << "\n"; csvP << "\n";

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
