#pragma once
#include <iostream>
#include "../circuit_modules/include/Resistor.h"
#include "../circuit_modules/include/VoltageSource.h"
class Circuit
{
private:
	std::vector<Node*> nodes;
	std::vector<Component*> components;

public:
	//Circuit(); 		// constructor
	~Circuit(); 	// destructor

	Node* addNode(int id);
	Resistor* addResistor(double resistance, Node* term1, Node* term2);
	VoltageSource* addVoltageSource(double volt, Node* term1, Node* term2);

	void saveData();
	void saveasNetlist();

	void resetCircuit();
	void simulateCircuit();

	void printResult();

	// Wipes the entire circuit grid completely
	void clearAll() 
	{
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
};

