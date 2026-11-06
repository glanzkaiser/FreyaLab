
#include "freyalab/Visualization.h"
#include "freyalab/Netlist.h"
#include <iostream>

int main() {
    std::cout << "Initializing FreyaLab Simulator Engine..." << std::endl;
    
    // Create the mandatory core engine Netlist dependency
    Netlist coreNetlist;
    
    // Initialize the visualizer workspace
    FreyaLabVisualizer visualizer(coreNetlist);
    
    std::cout << "Launching Native OpenGL Engineering Workspace Dashboard..." << std::endl;
    // Launch the GLFW loop context window
    visualizer.runInterfaceWorkspace();
    
    std::cout << "FreyaLab Application Terminated Smoothly." << std::endl;
    return 0;
}
