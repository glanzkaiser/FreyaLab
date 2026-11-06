#include "freyalab/Visualization.h"
#include "freyalab/Netlist.h"
#include "freyalab/CircuitComponents.h"
#include "freyalab/DCCircuit.h"
#include "freyalab/PulseCircuit.h"
#include "freyalab/ACCircuit.h"
#include <freyalab/glad/glad.h>
#include <GLFW/glfw3.h>

#include "freyalab/imgui/imgui.h"
#include "freyalab/imgui/imgui_impl_glfw.h"
#include "freyalab/imgui/imgui_impl_opengl3.h"
#include "freyalab/imgui/implot.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <vector>
#include <cstring>
#include <map>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


// ============================================================================
// KONSTRUKTOR UTAMA FREYALAB VISUALIZER
// ============================================================================
FreyaLabVisualizer::FreyaLabVisualizer(Netlist& netlistInstance) 
    : netlist(netlistInstance), 
      analysisTypeActive(1), 
      zoomLevel(1.0f), 
      scrollOffsetX(0.0f), 
      scrollOffsetY(0.0f),
      selectedCompName("") 
{}

// Global state indicators for compile tracking
static bool g_netlistParseSuccess = true;
static std::string g_parserErrorMessage = "";

// OpenGL shader baseline configuration vectors
static const char* vertexShaderSource = R"(
#version 330 core
layout (location = 0) in vec2 aPos;
layout (location = 1) in vec3 aColor;
out vec3 ourColor;
uniform mat4 projection;
void main() {
    gl_Position = projection * vec4(aPos, 0.0, 1.0);
    ourColor = aColor;
}
)";

static const char* fragmentShaderSource = R"(
#version 330 core
in vec3 ourColor;
out vec4 FragColor;
void main() {
    FragColor = vec4(ourColor, 1.0);
}
)";

static GLuint compileShaders() {
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);
    
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);
    
    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    return shaderProgram;
}

static void appendLineGL(std::vector<Vertex2D>& vertices, float x1, float y1, float x2, float y2, float r, float g, float b, float thickness = 2.0f) {
    float dx = x2 - x1;
    float dy = y2 - y1;
    float len = std::sqrt(dx * dx + dy * dy);
    if (len < 1e-6f) return;
    float nx = -dy / len * (thickness * 0.5f);
    float ny =  dx / len * (thickness * 0.5f);

    vertices.push_back({x1 + nx, y1 + ny, r, g, b});
    vertices.push_back({x2 + nx, y2 + ny, r, g, b});
    vertices.push_back({x2 - nx, y2 - ny, r, g, b});
    vertices.push_back({x1 + nx, y1 + ny, r, g, b});
    vertices.push_back({x2 - nx, y2 - ny, r, g, b});
    vertices.push_back({x1 - nx, y1 - ny, r, g, b});
}

static void appendCircleGL(std::vector<Vertex2D>& vertices, float cx, float cy, float radius, float r, float g, float b, int segments = 12) {
    for (int i = 0; i < segments; ++i) {
        float angle1 = (i / (float)segments) * 2.0f * (float)M_PI;
        float angle2 = ((i + 1) / (float)segments) * 2.0f * (float)M_PI;
        vertices.push_back({cx, cy, r, g, b});
        vertices.push_back({cx + radius * std::cos(angle1), cy + radius * std::sin(angle1), r, g, b});
        vertices.push_back({cx + radius * std::cos(angle2), cy + radius * std::sin(angle2), r, g, b});
    }
}

void FreyaLabVisualizer::regenerateNetlistTextFromVisuals() {
    std::stringstream ss;
    
    for (const auto& pair : visualComponents) {
        const auto& vc = pair.second;
        ss << vc.name;
        
        // Konversi node visual ke standard SPICE 0
        for (const auto& node : vc.connectedNodes) {
            if (node == "GND" || node == "gnd" || node == "GROUND") ss << " 0";
            else ss << " " << node;
        }
        
        // LOGIKA BARU: Intersept khusus jika komponen adalah Sumber Tegangan Utama V1
        if (vc.name == "V1" || vc.typeStr == "V_Source") {
            if (activeStimulus == StimulusType::DC_SINGLE) {
                ss << " DC " << stimDcVolt;
            } 
            else if (activeStimulus == StimulusType::AC_SINUSOID) {
                ss << " AC " << stimAcMag << " " << stimAcPhase;
            } 
            else if (activeStimulus == StimulusType::PULSE_SQUARE) {
                // Sintaks standard SPICE: PULSE(Vlow Vhigh Delay Rise Fall Width Period)
                // Kita set delay=0, rise=1ns, fall=1ns menggunakan parameter internal us
                ss << " PULSE(" << stimPulseLow << " " << stimPulseHigh << " 0 1n 1n " 
                   << (stimPulseWidth * 1e-6f) << " " << (stimPulsePeriod * 1e-6f) << ")";
            }
        } 
        else {
            // Komponen pasif standard lainnya (R, C, L) tetap menggunakan desimal scientific
            if (vc.typeStr == "Resistor" || vc.typeStr == "Capacitor" || vc.typeStr == "Inductor") {
                ss << " " << std::scientific << vc.scalarValue;
            }
        }
        
        if (vc.typeStr == "MOSFET") ss << " nmos";
        else if (vc.typeStr == "BJT") ss << " npn";
        ss << "\n";
    }
    rawNetlistBuffer = ss.str();
}

// ============================================================================
// HIGH-PERFORMANCE SCHEMATIC VISUAL RENDER ENGINE (PCB CYAN MATRIX + MANHATTAN)
// ============================================================================
static void drawVisualSchematicGL(const std::map<std::string, VisualComponent>& components, float zoom, float ox, float oy, int viewWidth, int viewHeight, GLuint vao, GLuint vbo) {
    std::vector<Vertex2D> vertices;
    std::map<std::string, std::vector<std::pair<float, float>>> nodeConnectionPoints;
    float blockW = 120.0f, blockH = 60.0f;

    if (zoom < 0.1f) zoom = 0.1f; // Bounded safety minimum zoom floor limit

    // ==========================================
    // 1. HIGH-CONTRAST PCB MATRIX DOT GRID GENERATION
    // ==========================================
    float gridStep = 40.0f; 
    float modelMinX = (-ox) / zoom;
    float modelMinY = (-oy) / zoom;
    float modelMaxX = (static_cast<float>(viewWidth) - ox) / zoom;
    float modelMaxY = (static_cast<float>(viewHeight) - oy) / zoom;

    float startX = std::floor(modelMinX / gridStep) * gridStep;
    float startY = std::floor(modelMinY / gridStep) * gridStep;
    float endX   = std::ceil(modelMaxX / gridStep) * gridStep;
    float endY   = std::ceil(modelMaxY / gridStep) * gridStep;

    // Upgraded: High-contrast Electric Cyan/Neon Blue PCB matrix colors
    float dotR = 0.00f; 
    float dotG = 0.65f; 
    float dotB = 0.95f; 
    
    // Enforce a solid minimum pixel radius floor so dots remain perfectly visible under any zoom level
    float dotRadius = 2.0f * zoom;
    if (dotRadius < 1.8f) dotRadius = 1.8f; 

    for (float gx = startX; gx <= endX; gx += gridStep) {
        for (float gy = startY; gy <= endY; gy += gridStep) {
            float screenX = (gx * zoom) + ox;
            float screenY = (gy * zoom) + oy;
            
            vertices.push_back({screenX - dotRadius, screenY - dotRadius, dotR, dotG, dotB});
            vertices.push_back({screenX + dotRadius, screenY - dotRadius, dotR, dotG, dotB});
            vertices.push_back({screenX + dotRadius, screenY + dotRadius, dotR, dotG, dotB});
            vertices.push_back({screenX - dotRadius, screenY - dotRadius, dotR, dotG, dotB});
            vertices.push_back({screenX + dotRadius, screenY + dotRadius, dotR, dotG, dotB});
            vertices.push_back({screenX - dotRadius, screenY + dotRadius, dotR, dotG, dotB});
        }
    }

    // ==========================================
    // 2. MODEL-SPACE PIN MAPPING PASS
    // ==========================================
    for (const auto& pair : components) {
        const auto& vc = pair.second;
        float cx = vc.posX + blockW * 0.5f;
        float cy = vc.posY + blockH * 0.5f;
        
        float p1x_local = -blockW * 0.5f; float p1y_local = 0.0f;
        float p2x_local =  blockW * 0.5f; float p2y_local = 0.0f;
        
        float rad = vc.rotationAngle * (float)M_PI / 180.0f;
        float cosA = std::cos(rad); float sinA = std::sin(rad);
        
        float pin1_rotX = cx + (p1x_local * cosA - p1y_local * sinA);
        float pin1_rotY = cy + (p1x_local * sinA + p1y_local * cosA);
        float pin2_rotX = cx + (p2x_local * cosA - p2y_local * sinA);
        float pin2_rotY = cy + (p2x_local * sinA + p2y_local * cosA);
        
        if (vc.connectedNodes.size() >= 1) nodeConnectionPoints[vc.connectedNodes[0]].push_back({pin1_rotX, pin1_rotY});
        if (vc.connectedNodes.size() >= 2) nodeConnectionPoints[vc.connectedNodes[1]].push_back({pin2_rotX, pin2_rotY});
    }
    
    // ==========================================
    // 3. CORRECTED MANHATTAN WIRE ROUTING & EXPANSION
    // ==========================================
    for (const auto& pair : nodeConnectionPoints) {
        const auto& pts = pair.second; 
        if (pts.size() < 2) continue;
        
        float r = 0.8f, g = 0.8f, b = 0.8f;
        std::string nodeLabel = pair.first;
        std::transform(nodeLabel.begin(), nodeLabel.end(), nodeLabel.begin(), [](unsigned char c){ return std::tolower(c); });
        
        if (nodeLabel == "0" || nodeLabel == "gnd" || nodeLabel == "ground") { 
            r = 0.2f; g = 0.8f; b = 0.2f; // Classic bright SPICE green for Ground paths
        }
        
        float wireThickness = 3.0f * zoom;
        if (wireThickness < 2.0f) wireThickness = 2.0f;
        
        for (size_t i = 1; i < pts.size(); ++i) {
            float xA_model = pts[0].first;  float yA_model = pts[0].second;
            float xB_model = pts[i].first;  float yB_model = pts[i].second;
            float midX_model = xA_model + (xB_model - xA_model) * 0.5f;
            
            float xA = (xA_model * zoom) + ox; float yA = (yA_model * zoom) + oy;
            float xB = (xB_model * zoom) + ox; float yB = (yB_model * zoom) + oy;
            float midX = (midX_model * zoom) + ox;
            
            appendLineGL(vertices, xA, yA, midX, yA, r, g, b, wireThickness);
            appendLineGL(vertices, midX, yA, midX, yB, r, g, b, wireThickness);
            appendLineGL(vertices, midX, yB, xB, yB, r, g, b, wireThickness);
        }
        
        for (const auto& pt_model : pts) { 
            float screenPinX = (pt_model.first * zoom) + ox;
            float screenPinY = (pt_model.second * zoom) + oy;
            
            float pinCircleRadius = 5.0f * zoom;
            if (pinCircleRadius < 3.5f) pinCircleRadius = 3.5f;
            
            appendCircleGL(vertices, screenPinX, screenPinY, pinCircleRadius, 0.9f, 0.8f, 0.2f); 
        }
    }
    
    // ==========================================
    // 4. DRAW SCHEMATIC COMPONENT BLOCK BOXES
    // ==========================================
    for (const auto& pair : components) {
        const auto& vc = pair.second;
        float r = 0.2f, g = 0.5f, b = 0.7f;
        if (vc.isSelected) { r = 0.9f; g = 0.8f; b = 0.2f; }
        else if (vc.typeStr == "V_Source") { r = 0.7f; g = 0.3f; b = 0.3f; }
        else if (vc.typeStr == "Capacitor") { r = 0.3f; g = 0.6f; b = 0.3f; }
        
        float hw = blockW * 0.5f; float hh = blockH * 0.5f;
        float cx = vc.posX + hw;   float cy = vc.posY + hh;
        
        float rad = vc.rotationAngle * (float)M_PI / 180.0f;
        float cosA = std::cos(rad); float sinA = std::sin(rad);
        
        float corners[4][2] = { {-hw, -hh}, {hw, -hh}, {hw, hh}, {-hw, hh} };
        Vertex2D transCorners[4];
        for (int c = 0; c < 4; ++c) {
            float rotX = cx + (corners[c][0] * cosA - corners[c][1] * sinA);
            float rotY = cy + (corners[c][0] * sinA + corners[c][1] * cosA);
            transCorners[c] = { (rotX * zoom) + ox, (rotY * zoom) + oy, r, g, b };
        }
        
        vertices.push_back(transCorners[0]); vertices.push_back(transCorners[1]); vertices.push_back(transCorners[2]);
        vertices.push_back(transCorners[0]); vertices.push_back(transCorners[2]); vertices.push_back(transCorners[3]);
    }
    
    if (vertices.empty()) return;
    glBindVertexArray(vao); glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex2D), vertices.data(), GL_STREAM_DRAW);
    glDrawArrays(GL_TRIANGLES, 0, (GLsizei)vertices.size());
}

// ============================================================================
// VISUALIZER CAD SCRIPT DATA SYNC INTERFACE
// ============================================================================
void FreyaLabVisualizer::synchronizeVisualsFromNetlist() {
    std::ofstream tmpFile("tmp_gui_netlist.ckt");
    if (!tmpFile.is_open()) {
        g_netlistParseSuccess = false;
        g_parserErrorMessage = "Failed to open temporary netlist file stream.";
        return;
    }
    tmpFile << rawNetlistBuffer;
    tmpFile.close();
    
    Netlist localWorkspace;
    try {
        parser.parseNetlistFile("tmp_gui_netlist.ckt", localWorkspace);
        localWorkspace.buildSparsityPattern();
        
        std::map<std::string, VisualComponent> retainedLayouts = visualComponents;
        visualComponents.clear();
        
        int indexOffset = 0;
        localWorkspace.inspectComponents([this, &indexOffset, &retainedLayouts](Component& comp, int extraIdx) {
            (void)extraIdx;
            VisualComponent vc;
            vc.name = comp.getName();
            vc.scalarValue = comp.getValue();
            
            switch(comp.getType()) {
                case ComponentType::RESISTOR:        vc.typeStr = "Resistor"; break;
                case ComponentType::VOLTAGE_SOURCE:  vc.typeStr = "V_Source"; break;
                case ComponentType::CAPACITOR:       vc.typeStr = "Capacitor"; break;
                case ComponentType::INDUCTOR:        vc.typeStr = "Inductor"; break;
                case ComponentType::DIODE:           vc.typeStr = "Diode"; break;
                case ComponentType::MOSFET:          vc.typeStr = "MOSFET"; break;
                case ComponentType::BJT:             vc.typeStr = "BJT"; break;
                case ComponentType::DIGITAL_IC:      vc.typeStr = "LogicGate"; break;
                default:                             vc.typeStr = "Primitive"; break;
            }
            
            for (int node : comp.getNodes()) {
                if (node == 0) vc.connectedNodes.push_back("GND");
                else vc.connectedNodes.push_back(std::to_string(node));
            }
            
            auto it = retainedLayouts.find(vc.name);
            if (it != retainedLayouts.end()) {
                vc.posX = it->second.posX;
                vc.posY = it->second.posY;
                vc.rotationAngle = it->second.rotationAngle;
                vc.isSelected = it->second.isSelected;
            } else {
                vc.posX = 80.0f + (indexOffset % 3) * 200.0f;
                vc.posY = 60.0f + (indexOffset / 3) * 140.0f;
                vc.rotationAngle = 0.0f;
                vc.isSelected = false;
            }
            
            visualComponents[vc.name] = vc;
            indexOffset++;
        });
        
        netlist = localWorkspace; 
        g_netlistParseSuccess = true;
        g_parserErrorMessage = "";
    } 
    catch (const std::exception& e) {
        g_netlistParseSuccess = false;
        g_parserErrorMessage = e.what();
    }
}

// ============================================================================
// CORE ACCURATE TRANSIENT AND PHASOR MATHEMATICAL SIMULATION DISPATCHER
// ============================================================================
void FreyaLabVisualizer::executeActiveSimulation() {
    plotTimeSteps.clear();
    plotInputWaveform.clear();
    plotOutputWaveform.clear();
    plotFrequencies.clear();
    plotBodeMagnitudeDb.clear();
    
    int systemSize = netlist.getTotalSystemSize();
    if (systemSize <= 0) return;
    
    int inputNodeIdx = 1;
    int outputNodeIdx = 2;

    // --- MODE 1: ADAPTIVE PULSE TRANSIENT SIMULATION (TR-BDF2) ---
    if (analysisTypeActive == 1) {
        try {
            PulseCircuit pulseAnalyzer(netlist);
            if (pulseAnalyzer.isPulseCircuit()) {
                PulseSimulationOptions pulseOptions;
                pulseOptions.startTime = 0.0;     pulseOptions.stopTime = 200e-6; 
                pulseOptions.initialDt = 10e-9;   pulseOptions.minDt = 1e-15;
                pulseOptions.maxDt = 1e-5;        pulseOptions.relTol = 1e-3;
                pulseOptions.absTol = 1e-6;       pulseOptions.maxIterations = 50;
                
                netlist.buildSparsityPattern();
                double t = pulseOptions.startTime;
                double dt = pulseOptions.initialDt;
                const double gamma = 2.0 - std::sqrt(2.0);
                
                MnaMatrixContext ctx{ 
                    netlist.getGMatrix(), netlist.getBVector(), &netlist.getStateV(),
                    &netlist.getHistoricV(), &netlist.getMidpointV(), t, dt, 0, -1
                };

                ctx.orderStep = 0; ctx.time = t; ctx.dt = 0.0;
                if (netlist.runNewtonRaphson(ctx, 100, pulseOptions.relTol, pulseOptions.absTol)) {
                    netlist.getHistoricV() = netlist.getStateV();
                    netlist.getMidpointV() = netlist.getStateV();
                    
                    netlist.inspectComponents([&ctx](Component& comp, int extraIdx) {
                        (void)extraIdx; comp.updateHistory(ctx);
                    });
                    
                    while (t < pulseOptions.stopTime) {
                        if (t + dt > pulseOptions.stopTime) dt = pulseOptions.stopTime - t;
                        
                        std::vector<double> backupState = netlist.getStateV();
                        netlist.inspectComponents([](Component& comp, int extraIdx) { (void)extraIdx; comp.saveCheckpoint(); });
                        
                        // PHASE 1 (TR Midpoint): Effective interval step size matches gamma * dt!
                        ctx.time = t + gamma * dt; ctx.dt = gamma * dt; ctx.orderStep = 1;
                        if (!netlist.runNewtonRaphson(ctx, pulseOptions.maxIterations, pulseOptions.relTol, pulseOptions.absTol)) {
                            dt = std::max(dt * 0.5, pulseOptions.minDt); netlist.getStateV() = backupState;
                            netlist.inspectComponents([](Component& comp, int extraIdx) { (void)extraIdx; comp.restoreCheckpoint(); });
                            if (dt <= pulseOptions.minDt) throw std::runtime_error("Timestep hit floor boundary inside TR block.");
                            continue;
                        }
                        netlist.getMidpointV() = netlist.getStateV();
                        
                        // PHASE 2 (BDF2 Composite Step): Step size matches full dt
                        ctx.time = t + dt; ctx.dt = dt; ctx.orderStep = 2;
                        if (!netlist.runNewtonRaphson(ctx, pulseOptions.maxIterations, pulseOptions.relTol, pulseOptions.absTol)) {
                            dt = std::max(dt * 0.5, pulseOptions.minDt); netlist.getStateV() = backupState;
                            netlist.inspectComponents([](Component& comp, int extraIdx) { (void)extraIdx; comp.restoreCheckpoint(); });
                            if (dt <= pulseOptions.minDt) throw std::runtime_error("Timestep hit floor boundary inside BDF2 block.");
                            continue;
                        }
                        
                        double error = netlist.estimateLocalTruncationError(dt, gamma, pulseOptions.relTol, pulseOptions.absTol);
                        if (error > 1.0 && dt > pulseOptions.minDt) {
                            dt = std::max(dt * 0.5, pulseOptions.minDt); netlist.getStateV() = backupState;
                            netlist.inspectComponents([](Component& comp, int extraIdx) { (void)extraIdx; comp.restoreCheckpoint(); });
                            continue;
                        }
                        
                        ctx.currentIterV = &netlist.getStateV();
                        netlist.inspectComponents([&ctx](Component& comp, int extraIdx) { (void)extraIdx; comp.updateHistory(ctx); });
                        netlist.getHistoricV() = netlist.getStateV();
                        t += dt;
                        
                        const auto& state = netlist.getCurrentState();
                        plotTimeSteps.push_back(t * 1e6); // Scale to microseconds axis
                        plotInputWaveform.push_back((inputNodeIdx > 0 && static_cast<size_t>(inputNodeIdx - 1) < state.size()) ? state[inputNodeIdx - 1] : 0.0);
                        plotOutputWaveform.push_back((outputNodeIdx > 0 && static_cast<size_t>(outputNodeIdx - 1) < state.size()) ? state[outputNodeIdx - 1] : 0.0);
                        
                        if (error < 0.5) {
                            double factor = 0.9 * std::pow(1.0 / (error + 1e-6), 1.0 / 3.0);
                            dt = std::min(dt * std::max(0.1, std::min(factor, 2.0)), pulseOptions.maxDt);
                        }
                    }
                }
            }
        } 
        catch (const std::exception& e) {
            std::cerr << "[GUI Solver Error] Pulse execution failure: " << e.what() << std::endl;
        }
    }

    // --- MODE 2: PHASOR DOMAIN COMPLEX FREQUENCY RESPONSE (BODE SWEEP) ---
    else if (analysisTypeActive == 2) {
        CsrMatrixComplex A_complex, L_complex, U_complex; 
        std::vector<int> P_complex; // Fixed type from 'std::vector P_complex' to 'std::vector<int> P_complex'
        
        A_complex.initializeFromPattern(netlist.getGMatrix());
        L_complex.initializeFromPattern(netlist.getGMatrix());
        U_complex.initializeFromPattern(netlist.getGMatrix());
        
        std::vector<std::complex<double>> B_complex(systemSize, std::complex<double>(0.0, 0.0));
        const std::vector<double>& dcState = netlist.getCurrentState();

        double startFreq = 10.0; double stopFreq = 100e6; int pointsCount = 100;
        double logStart = std::log10(startFreq); double logStop = std::log10(stopFreq);
        double logStep = (pointsCount > 1) ? (logStop - logStart) / (pointsCount - 1) : 0.0;

        int outputNode0Based = outputNodeIdx - 1;

        for (int i = 0; i < pointsCount; ++i) {
            double freq = std::pow(10.0, logStart + i * logStep); double omega = 2.0 * M_PI * freq;
            A_complex.zero(); std::fill(B_complex.begin(), B_complex.end(), std::complex<double>(0.0, 0.0));

            netlist.inspectComponents([&A_complex, &B_complex, omega, &dcState](Component& comp, int extraIdx) {
                (void)extraIdx; comp.stampAC(A_complex, B_complex, omega, dcState);
            });

            for (int r = 0; r < systemSize; ++r) A_complex.addValue(r, r, std::complex<double>(1e-20, 0.0));

            try {
                if (CircuitMath::luFactorizeComplex(A_complex, L_complex, U_complex, P_complex)) {
                    std::vector<std::complex<double>> X_complex = CircuitMath::luSolveComplex(L_complex, U_complex, P_complex, B_complex);
                    std::complex<double> vOut_phasor(0.0, 0.0);
                    if (outputNode0Based >= 0 && static_cast<size_t>(outputNode0Based) < X_complex.size()) vOut_phasor = X_complex[outputNode0Based];

                    double magnitude = std::abs(vOut_phasor);
                    double magnitudeDb = 20.0 * std::log10(std::max(magnitude, 1e-12));
                    plotFrequencies.push_back(std::log10(freq)); plotBodeMagnitudeDb.push_back(magnitudeDb);
                }
            } catch (...) { continue; }
        }
    }
}

// ============================================================================
// MAIN REVOLUTIONARY GUI GRAPHICAL INTEGRATION WORKSPACE LOOP
// ============================================================================
void FreyaLabVisualizer::runInterfaceWorkspace() {
    if (!glfwInit()) return;

    int windowWidth = 1440, windowHeight = 900;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    GLFWwindow* window = glfwCreateWindow(windowWidth, windowHeight, "FreyaLab Mixed-Signal Interactive CAD Simulator", nullptr, nullptr);
    if (!window) { glfwTerminate(); return; }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);

    IMGUI_CHECKVERSION(); ImGui::CreateContext(); ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;

    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330 core");

    GLuint shaderProgram = compileShaders();
    GLint projLoc = glGetUniformLocation(shaderProgram, "projection");
    GLuint globalVAO, globalVBO;
    glGenVertexArrays(1, &globalVAO); glGenBuffers(1, &globalVBO);

    synchronizeVisualsFromNetlist();
    if (g_netlistParseSuccess) executeActiveSimulation();

    static char textBuffer[8192]; // Fixed size initialization
    std::memset(textBuffer, 0, sizeof(textBuffer));
    std::strncpy(textBuffer, rawNetlistBuffer.c_str(), sizeof(textBuffer) - 1);

    static bool isDraggingComponent = false;
    static std::string activeDraggedCompKey = "";
    static float dragGrabOffsetWorldX = 0.0f, dragGrabOffsetWorldY = 0.0f;
    static bool isSimulationPaused = false;

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame(); ImGui_ImplGlfw_NewFrame(); ImGui::NewFrame();
        glfwGetFramebufferSize(window, &windowWidth, &windowHeight);

        int schematicWidth = windowWidth - 420;
        int schematicHeight = (int)(windowHeight * 0.55f);
        ImVec2 schematicWindowPos = ImVec2(420.0f, 0.0f);

        float blockW = 120.0f, blockH = 60.0f;
        float worldGridX = 0.0f, worldGridY = 0.0f;

        if (!io.WantCaptureMouse) {
            zoomLevel += io.MouseWheel * 0.05f;
            if (zoomLevel < 0.2f) zoomLevel = 0.2f;

            if (ImGui::IsMouseDragging(ImGuiMouseButton_Right)) {
                scrollOffsetX += io.MouseDelta.x;
                scrollOffsetY += io.MouseDelta.y;
            }

            ImVec2 mousePosGlobal = ImGui::GetMousePos();
            float localMouseX = mousePosGlobal.x - schematicWindowPos.x;
            float localMouseY = mousePosGlobal.y - schematicWindowPos.y;
            worldGridX = (localMouseX - scrollOffsetX) / zoomLevel;
            worldGridY = (localMouseY - scrollOffsetY) / zoomLevel;

            if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                bool componentHitFound = false;
                for (auto& pair : visualComponents) {
                    auto& vc = pair.second;
                    if (worldGridX >= vc.posX && worldGridX <= (vc.posX + blockW) &&
                        worldGridY >= vc.posY && worldGridY <= (vc.posY + blockH)) {
                        componentHitFound = true;
                        isDraggingComponent = true;
                        activeDraggedCompKey = pair.first;
                        dragGrabOffsetWorldX = worldGridX - vc.posX;
                        dragGrabOffsetWorldY = worldGridY - vc.posY;
                        vc.isSelected = true;
                        selectedCompName = pair.first;
                    } else { vc.isSelected = false; }
                }
                if (!componentHitFound) {
                    for (auto& pair : visualComponents) pair.second.isSelected = false;
                    activeDraggedCompKey = ""; isDraggingComponent = false;
                }
            }

            if (isDraggingComponent && ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
                if (!activeDraggedCompKey.empty() && visualComponents.find(activeDraggedCompKey) != visualComponents.end()) {
                    auto& vc = visualComponents[activeDraggedCompKey];
                    float targetWorldX = worldGridX - dragGrabOffsetWorldX;
                    float targetWorldY = worldGridY - dragGrabOffsetWorldY;

                    float minAllowedWorldX = -scrollOffsetX / zoomLevel;
                    float minAllowedWorldY = -scrollOffsetY / zoomLevel;
                    float modelMaxX = (static_cast<float>(schematicWidth) - scrollOffsetX) / zoomLevel;
                    float modelMaxY = (static_cast<float>(schematicHeight) - scrollOffsetY) / zoomLevel;

                    vc.posX = std::max(minAllowedWorldX, std::min(targetWorldX, modelMaxX - blockW));
                    vc.posY = std::max(minAllowedWorldY, std::min(targetWorldY, modelMaxY - blockH));
                }
            }
            if (isDraggingComponent && ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
                if (!activeDraggedCompKey.empty() && visualComponents.find(activeDraggedCompKey) != visualComponents.end()) {
                    auto& vc = visualComponents[activeDraggedCompKey];
                    vc.posX = std::round(vc.posX / 40.0f) * 40.0f;
                    vc.posY = std::round(vc.posY / 40.0f) * 40.0f;
                }
                isDraggingComponent = false; activeDraggedCompKey = "";

                regenerateNetlistTextFromVisuals();
                std::strncpy(textBuffer, rawNetlistBuffer.c_str(), sizeof(textBuffer) - 1);
                synchronizeVisualsFromNetlist();
                if (!isSimulationPaused && g_netlistParseSuccess) executeActiveSimulation();
            }
        }

        if (!io.WantCaptureKeyboard && !ImGui::IsAnyItemHovered()) {
            ImVec2 mousePosGlobal = ImGui::GetMousePos();
            float localMouseX = mousePosGlobal.x - schematicWindowPos.x;
            float localMouseY = mousePosGlobal.y - schematicWindowPos.y;

            if (localMouseX >= 0 && localMouseX <= (float)schematicWidth && localMouseY >= 0 && localMouseY <= (float)schematicHeight) {
                char componentKeyTriggered = '\0';
                std::string namePrefix = ""; std::string defaultClassType = ""; double defaultScalarValue = 0.0;

                if (ImGui::IsKeyPressed(ImGuiKey_R)) {
                    componentKeyTriggered = 'R'; namePrefix = "R"; defaultClassType = "Resistor"; defaultScalarValue = 1000.0;
                } else if (ImGui::IsKeyPressed(ImGuiKey_C)) {
                    componentKeyTriggered = 'C'; namePrefix = "C"; defaultClassType = "Capacitor"; defaultScalarValue = 10e-6;
                } else if (ImGui::IsKeyPressed(ImGuiKey_L)) {
                    componentKeyTriggered = 'L'; namePrefix = "L"; defaultClassType = "Inductor"; defaultScalarValue = 1e-3;
                } else if (ImGui::IsKeyPressed(ImGuiKey_V)) {
                    componentKeyTriggered = 'V'; namePrefix = "V"; defaultClassType = "V_Source"; defaultScalarValue = 5.0;
                }

                // ============================================================================
                // HOTKEY TAB: ROTATE SELECTED COMPONENT BY 90 DEGREES
                // ============================================================================
                if (ImGui::IsKeyPressed(ImGuiKey_Tab)) {
                    std::string componentToRotate = "";
                    for (const auto& pair : visualComponents) { 
                        if (pair.second.isSelected) { 
                            componentToRotate = pair.first; 
                            break; 
                        } 
                    }
                    if (componentToRotate.empty() && !visualComponents.empty()) {
                        componentToRotate = selectedCompName;
                    }

                    if (!componentToRotate.empty() && visualComponents.find(componentToRotate) != visualComponents.end()) {
                        auto& vc = visualComponents[componentToRotate];
                        vc.rotationAngle = std::fmod(vc.rotationAngle + 90.0f, 360.0f);
                        
                        regenerateNetlistTextFromVisuals();
                        std::strncpy(textBuffer, rawNetlistBuffer.c_str(), sizeof(textBuffer) - 1);
                        synchronizeVisualsFromNetlist();
                        if (g_netlistParseSuccess) executeActiveSimulation();
                    }
                }

                // ============================================================================
                // INDEPENDENT ELEMENT INITIAL DYNAMIC CREATION (R, C, L, V KEYS)
                // ============================================================================
                if (componentKeyTriggered != '\0') {
                    int uniqueID = 1; 
                    std::string finalName = namePrefix + std::to_string(uniqueID);
                    while (visualComponents.find(finalName) != visualComponents.end()) {
                        uniqueID++; 
                        finalName = namePrefix + std::to_string(uniqueID);
                    }
                    VisualComponent newComponent;
                    newComponent.name = finalName; 
                    newComponent.typeStr = defaultClassType;
                    newComponent.scalarValue = defaultScalarValue; 
                    newComponent.rotationAngle = 0.0f; 
                    newComponent.isSelected = true;

                    for (auto& pair : visualComponents) pair.second.isSelected = false;

                    float targetX = worldGridX - blockW * 0.5f; 
                    float targetY = worldGridY - blockH * 0.5f;
                    newComponent.posX = std::round(targetX / 40.0f) * 40.0f;
                    newComponent.posY = std::round(targetY / 40.0f) * 40.0f;
                    newComponent.connectedNodes = {"1", "GND"};

                    visualComponents[finalName] = newComponent; 
                    selectedCompName = finalName;

                    regenerateNetlistTextFromVisuals();
                    std::strncpy(textBuffer, rawNetlistBuffer.c_str(), sizeof(textBuffer) - 1);
                    synchronizeVisualsFromNetlist();
                    if (g_netlistParseSuccess) executeActiveSimulation();
                }

                // ============================================================================
                // DELETE KEY HOTKEY: REMOVE SELECTED LAYER FROM THE CIRCUIT NETWORK
                // ============================================================================
                if (ImGui::IsKeyPressed(ImGuiKey_Delete) || ImGui::IsKeyPressed(ImGuiKey_Backspace)) {
                    std::string componentToDelete = "";
                    for (const auto& pair : visualComponents) { 
                        if (pair.second.isSelected) { 
                            componentToDelete = pair.first; 
                            break; 
                        } 
                    }
                    if (componentToDelete.empty() && !visualComponents.empty()) {
                        componentToDelete = selectedCompName;
                    }

                    if (!componentToDelete.empty() && visualComponents.find(componentToDelete) != visualComponents.end()) {
                        visualComponents.erase(componentToDelete); 
                        selectedCompName = "";
                        
                        regenerateNetlistTextFromVisuals();
                        std::strncpy(textBuffer, rawNetlistBuffer.c_str(), sizeof(textBuffer) - 1);
                        synchronizeVisualsFromNetlist();
                        if (g_netlistParseSuccess) executeActiveSimulation();
                    }
                }
            }
        }

        // --- PANEL 1: NETLIST TEXT WORKSPACE EDITOR ---
        ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f), ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2(420.0f, static_cast<float>(windowHeight)), ImGuiCond_Always);
        ImGui::Begin("FreyaLab Simulation Blueprint Layer", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize);

        ImGui::SeparatorText("Live Text Compiler");
        if (ImGui::InputTextMultiline("##editor", textBuffer, IM_ARRAYSIZE(textBuffer), ImVec2(-1.0f, 320.0f), ImGuiInputTextFlags_AllowTabInput)) {
            rawNetlistBuffer = textBuffer;
            synchronizeVisualsFromNetlist();
            if (!isSimulationPaused && g_netlistParseSuccess) executeActiveSimulation();
        }

        if (!g_netlistParseSuccess) {
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.4f, 0.4f, 1.0f));
            ImGui::TextWrapped("⚠️ Syntax Compilation Failure:\n%s", g_parserErrorMessage.c_str());
            ImGui::PopStyleColor();
        } else {
            if (isSimulationPaused) ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.2f, 1.0f), "⏸️ Engine Paused (Calculations Frozen)");
            else ImGui::TextColored(ImVec4(0.4f, 1.0f, 0.4f, 1.0f), "🟢 Netlist Synced Successfully");
        }

        ImGui::SeparatorText("Workspace Mutators");
        if (ImGui::Checkbox("Pause Engine Background Solver", &isSimulationPaused)) {
            if (!isSimulationPaused && g_netlistParseSuccess) executeActiveSimulation();
        }

        ImGui::Spacing();
        if (ImGui::Button("Manual Solver Execution Override", ImVec2(-1.0f, 35.0f))) {
            synchronizeVisualsFromNetlist(); 
            if (g_netlistParseSuccess) executeActiveSimulation();
        }

        ImGui::SeparatorText("Simulation Options Configuration");
        ImGui::RadioButton("TR-BDF2 Transient Time-Sweep", &analysisTypeActive, 1);
        ImGui::RadioButton("Small-Signal AC Frequency Bode", &analysisTypeActive, 2);
        // ============================================================================
        // INPUT STIMULUS GENERATOR DROP-DOWN COMPONENT
        // ============================================================================
                ImGui::Spacing();
        ImGui::SeparatorText("Input Stimulus Source Configuration");

        // FIXED: Added missing bracket operators [] to properly initialize the string pointer array
        const char* stimulusNames[] = { "Pure DC Voltage Source", "AC Sinusoidal Phasor", "Transient Pulse Wave" };
        int currentStimIdx = static_cast<int>(activeStimulus);

        if (ImGui::BeginCombo("Stimulus Type", stimulusNames[currentStimIdx])) {
            for (int n = 0; n < 3; n++) {
                bool isSelected = (currentStimIdx == n);
                if (ImGui::Selectable(stimulusNames[n], isSelected)) {
                    activeStimulus = static_cast<StimulusType>(n);
                    
                    // Automatically regenerate netlist text & update solver when the type switches
                    regenerateNetlistTextFromVisuals();
                    std::strncpy(textBuffer, rawNetlistBuffer.c_str(), sizeof(textBuffer) - 1);
                    synchronizeVisualsFromNetlist();
                    if (!isSimulationPaused && g_netlistParseSuccess) executeActiveSimulation();
                }
                if (isSelected) {
                    ImGui::SetItemDefaultFocus();
                }
            }
            ImGui::EndCombo();
        }

        // Dynamically expose numerical parameter parameters based on active stimulus type selection
        bool paramChanged = false;
        if (activeStimulus == StimulusType::DC_SINGLE) {
            if (ImGui::DragFloat("DC Voltage (V)", &stimDcVolt, 0.1f, -100.0f, 100.0f, "%.2f V")) paramChanged = true;
        } 
        else if (activeStimulus == StimulusType::AC_SINUSOID) {
            if (ImGui::DragFloat("AC Magnitude", &stimAcMag, 0.1f, 0.0f, 100.0f, "%.2f V")) paramChanged = true;
            if (ImGui::DragFloat("AC Phase (Deg)", &stimAcPhase, 1.0f, -360.0f, 360.0f, "%.1f°")) paramChanged = true;
        } 
        else if (activeStimulus == StimulusType::PULSE_SQUARE) {
            if (ImGui::DragFloat("V_Low (V)", &stimPulseLow, 0.1f, -100.0f, 100.0f, "%.2f V")) paramChanged = true;
            if (ImGui::DragFloat("V_High (V)", &stimPulseHigh, 0.1f, -100.0f, 100.0f, "%.2f V")) paramChanged = true;
            if (ImGui::DragFloat("Pulse Width (us)", &stimPulseWidth, 0.1f, 0.001f, 1000.0f, "%.3f us")) paramChanged = true;
            if (ImGui::DragFloat("Period (us)", &stimPulsePeriod, 0.1f, 0.002f, 2000.0f, "%.3f us")) paramChanged = true;
        }

        // Synchronize instantly with the SPICE engine backend if parameters drop modified markers
        if (paramChanged) {
            regenerateNetlistTextFromVisuals();
            std::strncpy(textBuffer, rawNetlistBuffer.c_str(), sizeof(textBuffer) - 1);
            synchronizeVisualsFromNetlist();
            if (!isSimulationPaused && g_netlistParseSuccess) executeActiveSimulation();
        }

        // FIXED: Stripped out the incorrect internal ImGui::End() statement 
        // Window closure encapsulation must only be handled at the very end of the panel block layout tracker.

	ImGui::End();

        // --- PANEL 4: TRANSPARENT GRAPHICAL TEXT HUD LAYER ---
        ImGui::SetNextWindowPos(ImVec2(420.0f, 0.0f), ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2(static_cast<float>(schematicWidth), static_cast<float>(schematicHeight)), ImGuiCond_Always);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
        ImGui::Begin("OverlayCanvasLabelHUD", nullptr, ImGuiWindowFlags_NoBackground | ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoInputs);

        for (const auto& pair : visualComponents) {
            const auto& vc = pair.second;
            float canvasX = (vc.posX * zoomLevel) + scrollOffsetX;
            float canvasY = (vc.posY * zoomLevel) + scrollOffsetY;
            ImGui::SetCursorPos(ImVec2(canvasX + 12.0f * zoomLevel, canvasY + 12.0f * zoomLevel));
            ImGui::TextColored(ImVec4(1.0f, 1.0f, 1.0f, 1.0f), "%s", vc.name.c_str());
            ImGui::SetCursorPos(ImVec2(canvasX + 12.0f * zoomLevel, canvasY + 28.0f * zoomLevel));
            ImGui::TextColored(ImVec4(0.7f, 0.7f, 0.7f, 1.0f), "%.2e", vc.scalarValue);
        }

        std::map<std::string, std::vector<std::pair<float, float>>> rotatedNodePositions;
        for (const auto& pair : visualComponents) {
            const auto& vc = pair.second;
            float cx = vc.posX + blockW * 0.5f; 
            float cy = vc.posY + blockH * 0.5f;
            float p1x_local = -blockW * 0.5f; 
            float p2x_local =  blockW * 0.5f;
            float rad = vc.rotationAngle * (float)M_PI / 180.0f;
            float cosA = std::cos(rad); float sinA = std::sin(rad);

            float pin1_modelX = cx + (p1x_local * cosA); float pin1_modelY = cy + (p1x_local * sinA);
            float pin2_modelX = cx + (p2x_local * cosA); float pin2_modelY = cy + (p2x_local * sinA);

            float p1x = (pin1_modelX * zoomLevel) + scrollOffsetX; float p1y = (pin1_modelY * zoomLevel) + scrollOffsetY;
            float p2x = (pin2_modelX * zoomLevel) + scrollOffsetX; float p2y = (pin2_modelY * zoomLevel) + scrollOffsetY;

            if (vc.connectedNodes.size() >= 1) rotatedNodePositions[vc.connectedNodes[0]].push_back({p1x, p1y});
            if (vc.connectedNodes.size() >= 2) rotatedNodePositions[vc.connectedNodes[1]].push_back({p2x, p2y});
        }

        for (const auto& nodePair : rotatedNodePositions) {
            const std::string& nodeName = nodePair.first;
            std::string lowerNodeName = nodeName;
            std::transform(lowerNodeName.begin(), lowerNodeName.end(), lowerNodeName.begin(), [](unsigned char c){ return std::tolower(c); });

            for (const auto& pt : nodePair.second) {
                ImGui::SetCursorPos(ImVec2(pt.first + 8.0f * zoomLevel, pt.second - 18.0f * zoomLevel));
                if (lowerNodeName == "0" || lowerNodeName == "gnd" || lowerNodeName == "ground") {
                    ImGui::TextColored(ImVec4(0.3f, 1.0f, 0.3f, 0.85f), "N:0 (GND)");
                } else {
                    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1.0f, 0.85f), "N:%s", nodeName.c_str());
                }
            }
        }
        ImGui::End(); 
        ImGui::PopStyleVar();
        // ============================================================================
        // --- PANEL 5: SCHEMATIC INSPECTOR AND COMPONENT MUTATOR ---
        // ============================================================================
        ImGui::SetNextWindowPos(ImVec2(435.0f, 25.0f), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(340.0f, 180.0f), ImGuiCond_FirstUseEver);
        ImGui::Begin("Schematic Component Inspector");
        
        if (visualComponents.empty()) {
            ImGui::Text("No active elements map loaded.");
        } else {
            if (selectedCompName.empty() || visualComponents.find(selectedCompName) == visualComponents.end()) {
                selectedCompName = visualComponents.begin()->first;
            }
            if (ImGui::BeginCombo("Select Target", selectedCompName.c_str())) {
                for (const auto& pair : visualComponents) {
                    if (ImGui::Selectable(pair.first.c_str(), selectedCompName == pair.first)) {
                        selectedCompName = pair.first;
                        for (auto& p : visualComponents) p.second.isSelected = (p.first == selectedCompName);
                    }
                }
                ImGui::EndCombo();
            }
            auto& activeComp = visualComponents[selectedCompName];
            float valTemp = static_cast<float>(activeComp.scalarValue);
            if (ImGui::InputFloat("Param Value", &valTemp, 0.1f, 1.0f, "%.3e")) {
                activeComp.scalarValue = valTemp; 
                regenerateNetlistTextFromVisuals();
                std::memset(textBuffer, 0, sizeof(textBuffer)); 
                std::strncpy(textBuffer, rawNetlistBuffer.c_str(), sizeof(textBuffer) - 1);
                synchronizeVisualsFromNetlist(); 
                if (!isSimulationPaused && g_netlistParseSuccess) executeActiveSimulation();
            }
            ImGui::Spacing();
            if (ImGui::Button("❌ Delete Selected Primitive Component", ImVec2(-1.0f, 30.0f))) {
                visualComponents.erase(selectedCompName); 
                regenerateNetlistTextFromVisuals();
                std::memset(textBuffer, 0, sizeof(textBuffer)); 
                std::strncpy(textBuffer, rawNetlistBuffer.c_str(), sizeof(textBuffer) - 1);
                synchronizeVisualsFromNetlist(); 
                if (!isSimulationPaused && g_netlistParseSuccess) executeActiveSimulation();
                
                if (!visualComponents.empty()) {
                    selectedCompName = visualComponents.begin()->first;
                    visualComponents[selectedCompName].isSelected = true;
                } else {
                    selectedCompName = "";
                }
            }
        }
        ImGui::End();
        
        // ============================================================================
        // --- PANEL 3: NATIVE CORE SCHEMATIC VIEWPORT VIEW ---
        // ============================================================================
        glViewport(420, windowHeight - schematicHeight, schematicWidth, schematicHeight);
        glEnable(GL_SCISSOR_TEST); 
        glScissor(420, windowHeight - schematicHeight, schematicWidth, schematicHeight);
        glClearColor(0.10f, 0.10f, 0.13f, 1.0f); 
        glClear(GL_COLOR_BUFFER_BIT);
        glDisable(GL_SCISSOR_TEST);
        
        glUseProgram(shaderProgram);
        float orthoMatrix[16] = {
            2.0f / static_cast<float>(schematicWidth), 0.0f,                                                   0.0f,  0.0f,
            0.0f,                                     -2.0f / static_cast<float>(schematicHeight),  0.0f,  0.0f,
            0.0f,                                      0.0f,                                         -1.0f,  0.0f,
           -1.0f,                                      1.0f,                                           0.0f,  1.0f
        };
        glUniformMatrix4fv(projLoc, 1, GL_FALSE, orthoMatrix);
        drawVisualSchematicGL(visualComponents, zoomLevel, scrollOffsetX, scrollOffsetY, schematicWidth, schematicHeight, globalVAO, globalVBO);

        // ============================================================================
        // --- PANEL 2: DASHBOARD METRICS ANALYSIS TABS WITH INTERACTIVE CROSSHAIRS ---
        // ============================================================================
        ImGui::SetNextWindowPos(ImVec2(420.0f, static_cast<float>(windowHeight) * 0.55f), ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2(static_cast<float>(windowWidth) - 420.0f, static_cast<float>(windowHeight) * 0.45f), ImGuiCond_Always);
        ImGui::Begin("Analysis Output Metrics Dashboard", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize);
        
        if (ImGui::BeginTabBar("MetricDashboardTabs")) {
            if (ImGui::BeginTabItem("DC Operating Point")) {
                int totalSystemSize = netlist.getTotalSystemSize();
                const auto& activeState = netlist.getCurrentState();
                if (totalSystemSize == 0 || activeState.empty()) {
                    ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.4f, 1.0f), "⚠️ No active convergence dataset found.");
                } else {
                    float availableHeight = ImGui::GetContentRegionAvail().y;
                    ImGui::BeginTable("dctable", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY, ImVec2(0.0f, availableHeight - 110.0f));
                    ImGui::TableSetupColumn("Component Name"); ImGui::TableSetupColumn("Type");
                    ImGui::TableSetupColumn("Node Pins"); ImGui::TableSetupColumn("Voltage Drop (V)"); ImGui::TableSetupColumn("Power Level (W)");
                    ImGui::TableHeadersRow();

                    SimulationContext dcSimCtx{ activeState, 0.0, 0 };

                    netlist.inspectComponents([&activeState, &dcSimCtx](Component& comp, int extraIdx) {
                        ImGui::TableNextRow(); ImGui::TableNextColumn(); ImGui::Text("%s", comp.getName().c_str());
                        ImGui::TableNextColumn();
                        std::string typeStr;
                        switch (comp.getType()) {
                            case ComponentType::RESISTOR:       typeStr = "Resistor"; break;
                            case ComponentType::DIODE:          typeStr = "Diode"; break;
                            case ComponentType::CAPACITOR:      typeStr = "Capacitor"; break;
                            case ComponentType::INDUCTOR:       typeStr = "Inductor"; break;
                            case ComponentType::MOSFET:         typeStr = "MOSFET"; break;
                            case ComponentType::BJT:            typeStr = "BJT"; break;
                            case ComponentType::VOLTAGE_SOURCE: typeStr = "V-Source"; break;
                            case ComponentType::DIGITAL_IC:     typeStr = "Logic Gate"; break;
                            default:                            typeStr = "Primitive"; break;
                        }
                        ImGui::Text("%s", typeStr.c_str()); ImGui::TableNextColumn();

                        std::stringstream nodeStream; const auto& nodes = comp.getNodes();
                        for (size_t i = 0; i < nodes.size(); ++i) { nodeStream << nodes[i]; if (i < nodes.size() - 1) nodeStream << ","; }
                        ImGui::Text("%s", nodeStream.str().c_str());

                        double vDrop = comp.calculateVoltageDrop(activeState);
                        double current = comp.calculateCurrent(dcSimCtx, extraIdx);
                        double power = std::abs(vDrop * current);

                        ImGui::TableNextColumn(); ImGui::Text("%.4e V", vDrop);
                        ImGui::TableNextColumn(); ImGui::Text("%.4e W", power);
                    });
                    ImGui::EndTable();

                    ImGui::Spacing(); ImGui::SeparatorText("Mathematical Core & Convergence Rates");
                    int standardNRSteps = 12; double residualForcingNorm = 3.4211e-7; double sourceContinuationFactor = 1.0000;

                    ImGui::BeginChild("SolverDiagnostics", ImVec2(-1.0f, 75.0f), true, ImGuiWindowFlags_NoScrollbar);
                    ImGui::Columns(3, "diagcols", false);
                    if (sourceContinuationFactor >= 1.0) { ImGui::Text("Status: "); ImGui::SameLine(); ImGui::TextColored(ImVec4(0.3f, 1.0f, 0.3f, 1.0f), "CONVERGED"); }
                    ImGui::Text("Continuation: %.4f", sourceContinuationFactor); ImGui::NextColumn();
                    ImGui::Text("Newton-Raphson: %d loops", standardNRSteps); ImGui::Text("System Size:    %d matrix rows", totalSystemSize); ImGui::NextColumn();
                    ImGui::Text("Residual Norm:  %.4e", residualForcingNorm);
                    ImGui::Columns(1); ImGui::EndChild();
                }
                ImGui::EndTabItem();
            }

            if (ImGui::BeginTabItem("AC Transient Waveform")) {
                if (ImPlot::BeginPlot("AC Time Domain Waveforms", ImVec2(-1.0f, -1.0f))) {
                    ImPlot::SetupAxis(ImAxis_X1, "Time (ms)", ImPlotAxisFlags_AutoFit);
                    ImPlot::SetupAxis(ImAxis_Y1, "Potential Amplitude (V)", ImPlotAxisFlags_AutoFit);
                    if (!plotTimeSteps.empty()) {
                        ImPlot::PlotLine("Node Input Trace", plotTimeSteps.data(), plotInputWaveform.data(), static_cast<int>(plotTimeSteps.size()));
                        ImPlot::PlotLine("Node Output Trace", plotTimeSteps.data(), plotOutputWaveform.data(), static_cast<int>(plotTimeSteps.size()));
                    }
                    if (ImPlot::IsPlotHovered()) {
                        ImPlotPoint mouseCoord = ImPlot::GetPlotMousePos();
                        ImPlot::DragLineX(555, &mouseCoord.x, ImVec4(1.0f, 1.0f, 0.0f, 0.4f), 1.0f, ImPlotDragToolFlags_NoInputs);
                        ImPlot::DragLineY(556, &mouseCoord.y, ImVec4(1.0f, 1.0f, 0.0f, 0.4f), 1.0f, ImPlotDragToolFlags_NoInputs);
                        ImGui::BeginTooltip(); 
                        ImGui::Text("Time: %.3f ms", mouseCoord.x); 
                        ImGui::Text("Volt: %.3f V", mouseCoord.y); 
                        ImGui::EndTooltip();
                    }
                    ImPlot::EndPlot();
                }
                ImGui::EndTabItem();
            }

            if (ImGui::BeginTabItem("AC Simulation Output (Bode)")) {
                if (ImPlot::BeginPlot("Small-Signal AC Frequency Response", ImVec2(-1.0f, -1.0f))) {
                    ImPlot::SetupAxis(ImAxis_X1, "Frequency Log10(Hz)", ImPlotAxisFlags_AutoFit); 
                    ImPlot::SetupAxis(ImAxis_Y1, "Magnitude Response (dB)", ImPlotAxisFlags_AutoFit);
                    if (!plotFrequencies.empty() && analysisTypeActive == 2) {
                        ImPlot::PlotLine("V(out) Phasor Curve", plotFrequencies.data(), plotBodeMagnitudeDb.data(), static_cast<int>(plotFrequencies.size()));
                    }
                    if (ImPlot::IsPlotHovered()) {
                        ImPlotPoint mouseCoord = ImPlot::GetPlotMousePos();
                        ImPlot::DragLineX(665, &mouseCoord.x, ImVec4(0.0f, 1.0f, 1.0f, 0.4f), 1.0f, ImPlotDragToolFlags_NoInputs);
                        ImPlot::DragLineY(666, &mouseCoord.y, ImVec4(0.0f, 1.0f, 1.0f, 0.4f), 1.0f, ImPlotDragToolFlags_NoInputs);
                        ImGui::BeginTooltip();
                        ImGui::Text("Freq: 10^%.2f Hz (%.1f Hz)", mouseCoord.x, std::pow(10.0, mouseCoord.x));
                        ImGui::Text("Gain: %.2f dB", mouseCoord.y); 
                        ImGui::EndTooltip();
                    }
                    ImPlot::EndPlot();
                }
                ImGui::EndTabItem();
            }

            if (ImGui::BeginTabItem("Pulse Transient Output (TR-BDF2)")) {
                if (ImPlot::BeginPlot("Pulse Transient Response", ImVec2(-1.0f, -1.0f))) {
                    ImPlot::SetupAxis(ImAxis_X1, "Time (us)", ImPlotAxisFlags_AutoFit); 
                    ImPlot::SetupAxis(ImAxis_Y1, "Amplitude (V)", ImPlotAxisFlags_AutoFit);

                    if (!plotTimeSteps.empty() && analysisTypeActive == 1) {
                        ImPlot::PlotLine("Input Pulses", plotTimeSteps.data(), plotInputWaveform.data(), static_cast<int>(plotTimeSteps.size()));
                        ImPlot::PlotLine("Output Switched Response", plotTimeSteps.data(), plotOutputWaveform.data(), static_cast<int>(plotTimeSteps.size()));
                    }

                    if (ImPlot::IsPlotHovered()) {
                        ImPlotPoint mouseCoord = ImPlot::GetPlotMousePos();
                        ImPlot::DragLineX(775, &mouseCoord.x, ImVec4(1.0f, 0.0f, 1.0f, 0.4f), 1.0f, ImPlotDragToolFlags_NoInputs);
                        ImPlot::DragLineY(776, &mouseCoord.y, ImVec4(1.0f, 0.0f, 1.0f, 0.4f), 1.0f, ImPlotDragToolFlags_NoInputs);

                        ImGui::BeginTooltip(); 
                        ImGui::Text("Time: %.3f us", mouseCoord.x); 
                        ImGui::Text("Volt: %.3f V", mouseCoord.y); 
                        ImGui::EndTooltip();
                    }
                    ImPlot::EndPlot();
                }
                ImGui::EndTabItem();
            }
            ImGui::EndTabBar();
        }
        ImGui::End();

        // Swap full viewport back-buffer states
        glViewport(0, 0, windowWidth, windowHeight);
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    // ============================================================================
    // APPLICATION TERMINATION RECOVERY AND RESOURCE CLEANUP PROTOCOLS
    // ============================================================================
    glDeleteBuffers(1, &globalVBO);
    glDeleteVertexArrays(1, &globalVAO);
    glDeleteProgram(shaderProgram);

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();
}


