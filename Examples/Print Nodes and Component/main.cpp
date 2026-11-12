#include "freyalab/Netlist.h"
#include "freyalab/FileParser.h"
#include "freyalab/PulseCircuit.h"
#include <iostream>
#include <vector>
#include <utility>

int main() {
    Netlist myNetlist;
    FileParser myParser;

    // 1. Read the plain text netlist file 
    // Load the file layout - Node mappings are registered automatically inside the tracker maps
    myParser.parseNetlistFile("netlist.txt", myNetlist);

    // Print all component nodes and string labels securely
    myNetlist.printComponentNodes();
    myNetlist.printComponentInputOutputNodes();

    return 0;
}
