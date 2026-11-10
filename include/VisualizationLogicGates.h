#pragma once

#include "freyalab/Netlist.h"
#include "freyalab/FileParser.h" 
#include <freyalab/glad/glad.h> // Fixed: Direct inclusion of hardware primitive identifiers
#include <string>
#include <vector>
#include <map>

// Defined standard multi-channel 2D vertex footprint matching the OpenGL stream array strides
struct Vertex2D {
    float x, y;
    float r, g, b;
};

struct VisualComponent {
    std::string name;                        // Unique descriptor identifier matching netlist primitives
    std::string typeStr;                     // Device type name descriptor matching class signatures
    float posX = 0.0f;                       // Localized model space horizontal coordinate 
    float posY = 0.0f;                       // Localized model space vertical coordinate
    float rotationAngle = 0.0f;              // Orientation parameter transformation degree tracking offset
    bool isSelected = false;                 // State tracking highlighting context flags
    double scalarValue = 0.0;                // Real electrical rating scale metrics
    std::vector<std::string> connectedNodes; // Sequence tracking terminal connections
};

// ============================================================================
// HIGH-PERFORMANCE MULTI-STIMULUS INDEPENDENT SIGNAL SOURCE TYPES
// ============================================================================
enum class StimulusType {
    DC_SINGLE,      // Pure static constant DC voltage bias emitter
    AC_SINUSOID,    // Small-signal frequency phasor sweep carrier
    PULSE_SQUARE    // High-speed transient switching pulse square wave
};

class FreyaLabVisualizer {
public:
    // Constructor 1: Optimized for pointer-binding frameworks (Visualization.cpp)
    explicit FreyaLabVisualizer(Netlist* netlistInstance);
    
    // Constructor 2: Standard fallback initialization reference wrapper
    explicit FreyaLabVisualizer(Netlist& netlistInstance);

    ~FreyaLabVisualizer() = default;

    // Prevent copy mechanics to isolate hardware-bound OpenGL contexts safely
    FreyaLabVisualizer(const FreyaLabVisualizer&) = delete;
    FreyaLabVisualizer& operator=(const FreyaLabVisualizer&) = delete;

    // --- Public UI Dashboard and Interaction Execution Tracks ---
    void runInterfaceWorkspace();
    void executeActiveSimulation();
    void synchronizeVisualsFromNetlist();
    void regenerateNetlistTextFromVisuals();
    
    // Accessor required by main.cpp file layer boot sequences
    void setRawNetlistBuffer(const std::string& content) { rawNetlistBuffer = content; }
    
    void handleCanvasInteraction(float mouseX, float yMouse, int viewW, int viewH, 
                                 bool isClick, bool isDragging, float dragDx, float dragDy);

public:
    // ========================================================================
    // PUBLIC WORKSPACE STATE STORAGE CELLS & DATA REGISTERS
    // ========================================================================
    Netlist* netlistPtr = nullptr;
    FileParser parser;
    
    std::string rawNetlistBuffer;
    std::map<std::string, VisualComponent> visualComponents;
    std::string selectedCompName;

    int analysisTypeActive = 1; // 1 = Pulse/Digital table, 2 = AC Sweep
    float zoomLevel = 1.0f;
    float scrollOffsetX = 0.0f;
    float scrollOffsetY = 0.0f;

    // --- Interactive Vector Backing Stores (ImPlot / ImGui Grid Backends) ---
    std::vector<double> plotTimeSteps;
    std::vector<double> plotInputWaveform;
    std::vector<double> plotOutputWaveform;
    std::vector<double> plotFrequencies;
    std::vector<double> plotBodeMagnitudeDb;

private:
    // ========================================================================
    // BACKEND HARDWARE INFRASTRUCTURE REGISTERS FOR STIMULUS MULTIPLEXING
    // ========================================================================
    StimulusType activeStimulus = StimulusType::DC_SINGLE; 
    
    // Dynamic tracking cells for user parameter updates matching double precision layout sizes
    double stimDcVolt = 5.0;       
    double stimAcMag = 1.0;        
    double stimAcPhase = 0.0;      
    double stimPulseLow = 0.0;     
    double stimPulseHigh = 5.0;    
    double stimPulseWidth = 10e-6;  
    double stimPulsePeriod = 20e-6; 
};

