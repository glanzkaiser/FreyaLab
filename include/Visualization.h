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
     // Konstruktor 1: Untuk Visualization_Pointer.cpp (Menerima Pointer)
    explicit FreyaLabVisualizer(Netlist* netlistInstance);
    
    // Konstruktor 2: Untuk Visualization_Reference.cpp (Menerima Referensi)
    explicit FreyaLabVisualizer(Netlist& netlistInstance);

    ~FreyaLabVisualizer() = default;

    // --- Public UI Dashboard and Interaction Execution Tracks ---
    void runInterfaceWorkspace();
    void executeActiveSimulation();
    void synchronizeVisualsFromNetlist();
    void regenerateNetlistTextFromVisuals();
    
    // FIX 2: Restored missing interactive user input canvas transformation method parameters
    void handleCanvasInteraction(float mouseX, float yMouse, int viewW, int viewH, 
                                 bool isClick, bool isDragging, float dragDx, float dragDy);
    void renderDashboardGUI();

public:
    // ========================================================================
    // PUBLIC WORKSPACE STATE STORAGE CELLS & DATA REGISTERS
    // ========================================================================
    Netlist& netlist;
    FileParser parser;
    
    std::string rawNetlistBuffer;
    std::map<std::string, VisualComponent> visualComponents;
    std::string selectedCompName;

    int analysisTypeActive; // 1 = Pulse, 2 = AC Sweep
    float zoomLevel;
    float scrollOffsetX;
    float scrollOffsetY;

    // --- Interactive Vector Backing Stores (ImPlot Backends) ---
    std::vector<double> plotTimeSteps;
    std::vector<double> plotInputWaveform;
    std::vector<double> plotOutputWaveform;
    std::vector<double> plotFrequencies;
    std::vector<double> plotBodeMagnitudeDb;

private:
    // FIX 3: Restored necessary background geometry builder declaration
    void buildSchematicGeometry(const Netlist& targetNetlist, std::vector<Vertex2D>& geometryVertices);

private:
    // ========================================================================
    // BACKEND HARDWARE INFRASTRUCTURE REGISTERS FOR STIMULUS MULTIPLEXING
    // ========================================================================
    StimulusType activeStimulus = StimulusType::DC_SINGLE; // Default baseline state cell
    
    // Dynamic tracking cells for user parameter updates inside ImGui controls
    float stimDcVolt = 5.0f;       // Volts parameter scale
    float stimAcMag = 1.0f;        // Small-signal Magnitude peak scale
    float stimAcPhase = 0.0f;      // Phasor phase alignment in degrees metric
    float stimPulseLow = 0.0f;     // Low voltage floor level state (V)
    float stimPulseHigh = 5.0f;    // High voltage ceiling level state (V)
    float stimPulseWidth = 10e-6f;  // Pulse length spacing duration tracking (10 us)
    float stimPulsePeriod = 20e-6f; // Signal repeat rate loop metric tracking (20 us)
};

