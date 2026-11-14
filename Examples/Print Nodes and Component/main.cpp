#include "freyalab/Netlist.h"
#include "freyalab/FileParser.h"
#include "freyalab/PulseCircuit.h"
#include "freyalab/CircuitComponents.h"
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

int main() {
    try {
        Netlist myNetlist;
        FileParser myParser;

        // 1. Parse netlist.txt (registers NetNode IDs and string labels automatically)
        myParser.parseNetlistFile("netlist.txt", myNetlist);
        myNetlist.buildSparsityPattern();

        // 2. Print all component NetNode mappings and terminal I/O roles
        myNetlist.printComponentNodes();
        myNetlist.printComponentInputOutputNodes();

        // 3. Auto-discover type-safe differential NetNode pairs for pulse monitoring
        PulseCircuit pulseAnalyzer(myNetlist);
        if (pulseAnalyzer.isPulseCircuit()) {
            std::vector<std::pair<NetNode, NetNode>> monitorPairs;
            pulseAnalyzer.autoGenerateMonitorPairs(monitorPairs);

            const std::vector<double>& currentState = myNetlist.getCurrentState();
            std::cout << "Monitored Differential NetNode Pairs:\n";
            for (const auto& pair : monitorPairs) {
                double vDiff = getSafeVoltage(currentState, pair.first) -
                               getSafeVoltage(currentState, pair.second);
                std::cout << "  V(" << myNetlist.getNodeName(pair.first)
                          << ", " << myNetlist.getNodeName(pair.second) << ") = "
                          << std::fixed << std::setprecision(4) << vDiff << " V\n";
            }
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 1;
    }

    return 0;
}
