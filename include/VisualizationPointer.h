#pragma once

#include "freyalab/Netlist.h"
#include "freyalab/FileParser.h" // Fixed: Secure inclusion path tracking
#include <string>
#include <vector>
#include <map>

// Symmetrical definition guard for OpenGL context bindings cleanly isolated
#ifndef GLUINT_DEFINED
typedef unsigned int GLuint;
#define GLUINT_DEFINED
#endif

// FIX 1: Declared standard multi-channel 2D vertex footprint before usage tracks
struct Vertex2D {
    float x, y;
    float r, g, b;
};

struct VisualComponent {
    std::string name;                     // Unique descriptor identifier matching netlist primitives
    std::string typeStr;                  // Device type name descriptor matching class signatures
    float posX = 0.0f;                    // Localized model space horizontal coordinate 
    float posY = 0.0f;                    // Localized model space vertical coordinate
    float rotationAngle = 0.0f;           // Orientation parameter transformation degree tracking offset
    bool isSelected = false;              // State tracking highlighting context flags
    double scalarValue = 0.0;             // Real electrical rating scale metrics
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
    // Constructor 1: Optimized for pointer-binding frameworks (e.g., Visualization_Pointer.cpp)
    explicit FreyaLabVisualizer(Netlist* netlistInstance);
    
    // Constructor 2: Standard fallback initialization reference wrapper (e.g., Visualization_Reference.cpp)
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
    
    void handleCanvasInteraction(float mouseX, float yMouse, int viewW, int viewH, 
                                 bool isClick, bool isDragging, float dragDx, float dragDy);
    void renderDashboardGUI();

public:
    // ========================================================================
    // PUBLIC WORKSPACE STATE STORAGE CELLS & DATA REGISTERS
    // ========================================================================
    // FIXED LAYER: Backing variables unified into a single pointer format to clear 

    Netlist* netlistPtr = nullptr;
    FileParser parser;
    
    std::string rawNetlistBuffer;
    std::map<std::string, VisualComponent> visualComponents;
    std::string selectedCompName;

    int analysisTypeActive = 1; // 1 = Pulse (TR-BDF2), 2 = AC Sweep (Bode)
    float zoomLevel = 1.0f;
    float scrollOffsetX = 0.0f;
    float scrollOffsetY = 0.0f;

    // --- Interactive Vector Backing Stores (ImPlot Backends) ---
    std::vector<double> plotTimeSteps;
    std::vector<double> plotInputWaveform;
    std::vector<double> plotOutputWaveform;
    std::vector<double> plotFrequencies;
    std::vector<double> plotBodeMagnitudeDb;

private:
    void buildSchematicGeometry(const Netlist& targetNetlist, std::vector<Vertex2D>& geometryVertices);

private:
    // ========================================================================
    // BACKEND HARDWARE INFRASTRUCTURE REGISTERS FOR STIMULUS MULTIPLEXING
    // ========================================================================
    StimulusType activeStimulus = StimulusType::DC_SINGLE; 
    
    // Dynamic tracking cells for user parameter updates inside ImGui controls
    float stimDcVolt = 5.0f;       
    float stimAcMag = 1.0f;        
    float stimAcPhase = 0.0f;      
    float stimPulseLow = 0.0f;     
    float stimPulseHigh = 5.0f;    
    float stimPulseWidth = 10e-6f;  
    float stimPulsePeriod = 20e-6f; 
};


