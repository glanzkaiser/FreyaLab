
#include "freyalab.h"

Circuit::~Circuit() // Destructor to automatically clean up memory when circuit goes out of scope
{
	/*for (int i = 0; i < nodes.size(); ++i)
		delete nodes[i];

	for (int i = 0; i < components.size(); ++i)
		delete components[i];*/
	resetCircuit();
}

Node* Circuit::addNode(int id)
{
	Node* node = new Node(id);
	nodes.push_back(node);
	return node;
}

Resistor* Circuit::addResistor(double resistance, Node* term1, Node* term2)
{
	Resistor* res = new Resistor(resistance, term1, term2);
	components.push_back(res);
	return res;
}

VoltageSource* Circuit::addVoltageSource(double volt, Node* term1, Node* term2)
{
	VoltageSource* source = new VoltageSource(volt, term1, term2);
	components.push_back(source);
	return source;
}

void Circuit::saveData()
{
	// 1. Get current time
	auto now = std::chrono::system_clock::now();
	auto in_time_t = std::chrono::system_clock::to_time_t(now);

	// 2. Format the time (e.g., save_2026-09-23_05-30-00.txt)
	std::stringstream ss;
	ss << "EC_" << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d_%H-%M-%S") << ".txt";
	std::string filename = ss.str();
	int numNodes = nodes.size();

	Vector2D mat(numNodes, numNodes);
	std::vector<double> b(numNodes);

	// 3. Save the file
	std::ofstream outFile(filename);
	if (outFile.is_open()) 
	{
		int i = 0, j=0;
		for (Component* comp : components)
		{
			if (Resistor* resistor = dynamic_cast<Resistor*>(comp))
			{
				outFile << "R" << i << " "<< resistor->getValue()  << endl;
				i = i+1;
			}
			if (VoltageSource* voltagesource = dynamic_cast<VoltageSource*>(comp))
			{
				outFile << "V" << j << " "<< voltagesource->getValue()  << endl;
				j = j+1;
			}
			
		}

		for (Component* comp : components)
		{
			if (Resistor* resistor = dynamic_cast<Resistor*>(comp))
			{
				outFile << "\n\nThe Voltage Drop on Resistor " << resistor->getValue() << " ohm is: " << resistor->getVoltageDrop() << " V " << endl;
			}
			if (Resistor* resistor = dynamic_cast<Resistor*>(comp))
			{
				outFile << "The current on Resistor " << resistor->getValue() << " ohm is: " << resistor->getCurrent() << " ampere" << endl;
			}

			if (Resistor* resistor = dynamic_cast<Resistor*>(comp))
			{
				Node* term1 = resistor->getTerminal1();
				Node* term2  = resistor->getTerminal2();
				int r = term1->getId();
				int c = term2->getId();

				outFile << "Terminal 1 ID of Resistor " << resistor->getValue() << " = " << r << endl;
				outFile << "Terminal 2 ID of Resistor " << resistor->getValue() << " = " << c << endl;
			}
		}
		
		for (Component* comp : components)
		{
			if (VoltageSource* voltagesource = dynamic_cast<VoltageSource*>(comp))
			{
				outFile << "\n\nThe Voltage Drop on Voltage Source " << voltagesource->getValue() << " V is: " << voltagesource->getVoltageDrop() << " V " << endl;
			}

			
			if (VoltageSource* voltagesource = dynamic_cast<VoltageSource*>(comp))
			{
				Node* term1 = voltagesource->getTerminal1();
				Node* term2  = voltagesource->getTerminal2();
				int r = term1->getId();
				int c = term2->getId();

				outFile << "Terminal 1 ID of Voltage Source " << voltagesource->getValue() << " = " << r << endl;
				outFile << "Terminal 2 ID of Voltage Source " << voltagesource->getValue() << " = " << c << endl;
			}
		}
		int numNodes = nodes.size();
		Vector2D mat(numNodes, numNodes);
		vector<double> b(numNodes);

		for (Component* comp : components)
		{
			comp->fillMNA(mat, b);
		}
		mat.solveEquations(b);

		for (int i = 0; i < numNodes; ++i)	
		{
			nodes[i]->setVolt(b[i]);
			outFile << "Nodes id" << nodes[i]->getId() << " , volt = "<< nodes[i]->getVolt()  << endl;
		}
		outFile << "\nTest Save FreyaLab Circuit.\n";
		outFile.close();
		cout << "Saved successfully as: " << filename << endl;
	}
}

void Circuit::saveasNetlist() // function to save drawn components as netlist
{
	// 1. Get current time
	auto now = std::chrono::system_clock::now();
	auto in_time_t = std::chrono::system_clock::to_time_t(now);

	// 2. Format the time (e.g., save_2026-09-23_05-30-00.txt)
	std::stringstream ss;
	ss << "netlist_" << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d_%H-%M-%S") << ".txt";
	std::string filename = ss.str();
	int numNodes = nodes.size();

	Vector2D mat(numNodes, numNodes);
	std::vector<double> b(numNodes);

	// 3. Save the file
	std::ofstream outFile(filename);
	if (outFile.is_open()) 
	{
		int i = 0, j=0;
		for (Component* comp : components)
		{
			if (Resistor* resistor = dynamic_cast<Resistor*>(comp))
			{
				Node* term1 = resistor->getTerminal1();
				Node* term2  = resistor->getTerminal2();
				int r = term1->getId();
				int c = term2->getId();
				outFile << "R" << i << " " << r  << " " << c << " " << resistor->getValue() << endl;
				i = i+1;
			}
			if (VoltageSource* voltagesource = dynamic_cast<VoltageSource*>(comp))
			{
				Node* term1 = voltagesource->getTerminal1();
				Node* term2  = voltagesource->getTerminal2();
				int r = term1->getId();
				int c = term2->getId();
				outFile << "V" << j << " " << r  << " " << c << " " << voltagesource->getValue()  << endl;
				j = j+1;
			}
			
		}

		outFile.close();
		cout << "Saved successfully as: " << filename << endl;
	}
}

void Circuit::resetCircuit() // function to delete all created nodes and components
{
	// Clean up previous simulation layout
	// 1. Delete all dynamically allocated components
	for (Component* comp : components) 
	{
		if (comp != nullptr) 
		{
			delete comp;
		}
	}
	components.clear(); // Clear the vector elements

	// 2. Delete all dynamically allocated nodes
	for (Node* node : nodes) 
	{
		if (node != nullptr) 
		{
			delete node;
		}
	}
	nodes.clear(); // Clear the vector elements
}


void Circuit::simulateCircuit()
{
	int numNodes = nodes.size();

	Vector2D mat(numNodes, numNodes);
	std::vector<double> b(numNodes);

	for (Component* comp : components)
		comp->fillMNA(mat, b);

	mat.solveEquations(b);

	for (int i = 0; i < numNodes; ++i)
		nodes[i]->setVolt(b[i]);
}

void Circuit::printResult()
{
	for (Component* comp : components)
	{
		if (Resistor* resistor = dynamic_cast<Resistor*>(comp))
		{
			cout << "The Voltage Drop on Resistor " << resistor->getValue() << " ohm is: " << resistor->getVoltageDrop() << " V\n";
		}
	}
}