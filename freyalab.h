
// freyalab.h

// include some important headers
#include <SFML/Graphics.hpp>

#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <bits/stdc++.h>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

// phased include headers
// according to class hierarchy
#include "include/Circuit.h"
#include "include/ParseFile.h"

#include "circuit_modules/include/Component.h"
#include "circuit_modules/include/Node.h"
#include "circuit_modules/include/Resistor.h"
#include "circuit_modules/include/VoltageSource.h"

#include "GUI/include/App.h"
#include "GUI/include/AssetManager.h"
#include "GUI/include/DEFINITION.h"
#include "GUI/include/InputManager.h"
#include "GUI/include/MainMenuState.h"
#include "GUI/include/ResistorTexture.h"
#include "GUI/include/SimulationState.h"
#include "GUI/include/SplashState.h"
#include "GUI/include/State.h"
 #include "GUI/include/StateMachine.h"
#include "GUI/include/Wire.h"

#include "LinearAlgebra/include/algebra.h"

#include"circuit_modules/include/Component.h"
using namespace std;