#include "freyalab/LogicGates.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <chrono>

using namespace std::chrono;
using namespace std;

struct LogicTimeStep 
{
    std::string timeScale;
    uint8_t inputA; 
    uint8_t inputB;
};

int main() {
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	// Declare the input vectors exactly matching your requested timeline
	std::vector<LogicTimeStep> simulationSteps = {
		{ "10ns", 0b0000, 0b0000 },
		{ "20ns", 0b1010, 0b0101 },
		{ "30ns", 0b1100, 0b1010 },
		{ "40ns", 0b1111, 0b1111 }
	};

	// Print out headers exactly matching your target format
	std::cout << std::left << std::setw(9) << "Time" 
	<< std::setw(4) << "A" 
	<< std::setw(6) << " B" << "| "
	<< std::setw(5) << "NOT" 
	<< std::setw(5) << "AND" 
	<< std::setw(5) << "OR" 
	<< std::setw(5) << "XOR" 
	<< std::setw(5) << "NAND" 
	<< std::setw(5) << "NOR" 
	<< "XNOR\n";
	std::cout << "-------------------------------------------------\n";

	// Loop through timescales, dispatching calculations dynamically via the GateType enum
	for (const auto& step : simulationSteps) 
	{
		uint8_t a = step.inputA;
		uint8_t b = step.inputB;

		// FIXED: Mengubah evaluateGate menjadi evaluateBusGate untuk pemrosesan paralel 4-bit
		std::cout << std::left << std::setw(9) << step.timeScale
			<< std::setw(4) << LogicGates::toBinary4Bit(a) << " "
			<< std::setw(6) << LogicGates::toBinary4Bit(b) << "| "
			<< std::setw(5) << LogicGates::toBinary4Bit(LogicGates::evaluateBusGate(GateType::NOT, a))
			<< std::setw(5) << LogicGates::toBinary4Bit(LogicGates::evaluateBusGate(GateType::AND, a, b))
			<< std::setw(5) << LogicGates::toBinary4Bit(LogicGates::evaluateBusGate(GateType::OR, a, b))
			<< std::setw(5) << LogicGates::toBinary4Bit(LogicGates::evaluateBusGate(GateType::XOR, a, b))
			<< std::setw(5) << LogicGates::toBinary4Bit(LogicGates::evaluateBusGate(GateType::NAND, a, b))
			<< std::setw(5) << LogicGates::toBinary4Bit(LogicGates::evaluateBusGate(GateType::NOR, a, b))
			<< LogicGates::toBinary4Bit(LogicGates::evaluateBusGate(GateType::XNOR, a, b)) << "\n";
	}
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	return 0;
}
