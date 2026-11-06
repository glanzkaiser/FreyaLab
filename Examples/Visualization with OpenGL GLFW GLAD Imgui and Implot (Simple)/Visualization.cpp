#include "freyalab/Visualization.h"
#include "freyalab/Netlist.h"
#include "freyalab/CircuitComponents.h"
#include "freyalab/DCCircuit.h"
#include "freyalab/PulseCircuit.h"
#include "freyalab/ACCircuit.h"
#include <freyalab/glad/glad.h>
#include <GLFW/glfw3.h>

// ImGui v1.93 & ImPlot v1.1 Integration Headers
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

// ============================================================================
// FREYALABVISUALIZER CLASS MEMBER IMPLEMENTATIONS
// ============================================================================

FreyaLabVisualizer::FreyaLabVisualizer(Netlist& netlistInstance) : netlist(netlistInstance) {
    rawNetlistBuffer = "* New FreyaLab Interactive Project\nV1 1 0 5.0\nR1 1 2 1k\nR2 1 2 5k\nC1 2 0 10u\nR3 3 0 10k\nV2 2 3 2.5u\n";
    analysisTypeActive = 1; // Default to Transient Sweep
}

void FreyaLabVisualizer::synchronizeVisualsFromNetlist() {
    visualComponents.clear();
    std::ofstream tmpFile("tmp_gui_netlist.ckt");
    if (!tmpFile.is_open()) {
        std::cerr << "Error: Unable to create temporary netlist cache configuration." << std::endl;
        return;
    }
    tmpFile << rawNetlistBuffer;
    tmpFile.close();
    
    netlist = Netlist();
    try {
        parser.parseNetlistFile("tmp_gui_netlist.ckt", netlist);
        netlist.buildSparsityPattern();
        int indexOffset = 0;
        
        netlist.inspectComponents([this, &indexOffset](const Component& comp, int extraIdx) {
            (void)extraIdx;
            VisualComponent vc;
            vc.name = comp.getName();
            vc.typeStr = "Unknown";
            vc.scalarValue = comp.getValue();
            
            switch(comp.getType()) {
                case ComponentType::RESISTOR:        vc.typeStr = "Resistor"; break;
                case ComponentType::VOLTAGE_SOURCE:  vc.typeStr = "V_Source"; break;
                case ComponentType::CAPACITOR:       vc.typeStr = "Capacitor"; break;
                case ComponentType::INDUCTOR:        vc.typeStr = "Inductor"; break;
                case ComponentType::DIODE:           vc.typeStr = "Diode"; break;
                case ComponentType::DIGITAL_IC:      vc.typeStr = "LogicGate"; break;
                default: break;
            }
            
            for (int node : comp.getNodes()) {
                vc.connectedNodes.push_back(std::to_string(node));
            }
            
            vc.posX = 50.0f + (indexOffset % 2) * 200.0f;
            vc.posY = 50.0f + (indexOffset / 2) * 120.0f;
            visualComponents[vc.name] = vc;
            indexOffset++;
        });
    } 
    catch (const std::exception& e) {
        std::cerr << "Visual Synchronizer Matrix Drop: " << e.what() << std::endl;
    }
}

void FreyaLabVisualizer::executeActiveSimulation() {
    plotTimeSteps.clear();
    plotInputWaveform.clear();
    plotOutputWaveform.clear();
    plotFrequencies.clear();
    plotBodeMagnitudeDb.clear();
    
    int systemSize = netlist.getTotalSystemSize();
    if (systemSize <= 0) return;
    
    if (analysisTypeActive == 1) {
        double t = 0.0, stopT = 1e-3;
        int step = 0;
        while (t < stopT && step < 200) {
            plotTimeSteps.push_back(t * 1e6);
            plotInputWaveform.push_back(5.0 * std::sin(2.0 * M_PI * 1000.0 * t));
            plotOutputWaveform.push_back(4.2 * (1.0 - std::exp(-t / 1e-4)));
            t += 5e-5;
            step++;
        }
    } 
    else if (analysisTypeActive == 2) {
        CsrMatrixComplex A_complex, L_complex, U_complex;
        std::vector<int> P_complex;
        A_complex.initializeFromPattern(netlist.getGMatrix());
        std::vector<std::complex<double>> B_complex(systemSize, std::complex<double>(0.0, 0.0));
        const std::vector<double>& dcState = netlist.getCurrentState();

        double startFreq = 10.0;
        double stopFreq = 100e6;
        int pointsCount = 100;
        double logStart = std::log10(startFreq);
        double logStop = std::log10(stopFreq);
        double logStep = (pointsCount > 1) ? (logStop - logStart) / (pointsCount - 1) : 0.0;
        
        for (int i = 0; i < pointsCount; ++i) {
            double freq = std::pow(10.0, logStart + i * logStep);
            double omega = 2.0 * M_PI * freq;
            A_complex.zero();
            std::fill(B_complex.begin(), B_complex.end(), std::complex<double>(0.0, 0.0));
            
            netlist.inspectComponents([&A_complex, &B_complex, omega, &dcState](const Component& comp, int extraIdx) {
                (void)extraIdx;
                const_cast<Component&>(comp).stampAC(A_complex, B_complex, omega, dcState);
            });
            
            try {
                if (CircuitMath::luFactorizeComplex(A_complex, L_complex, U_complex, P_complex)) {
                    std::vector<std::complex<double>> X_complex = CircuitMath::luSolveComplex(L_complex, U_complex, P_complex, B_complex);
                    std::complex<double> vOut_phasor(0.0, 0.0);
                    if (X_complex.size() >= 2) {
                        vOut_phasor = X_complex[1];
                    } else if (!X_complex.empty()) {
                        vOut_phasor = X_complex[0];
                    }
                    double magnitude = std::abs(vOut_phasor);
                    double magnitudeDb = 20.0 * std::log10(std::max(magnitude, 1e-12));
                    plotFrequencies.push_back(std::log10(freq));
                    plotBodeMagnitudeDb.push_back(magnitudeDb);
                }
            } 
            catch (...) {
                continue;
            }
        }
    }
}

void FreyaLabVisualizer::runInterfaceWorkspace() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return;
    }
    
    int windowWidth = 1440;
    int windowHeight = 900;
    
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    GLFWwindow* window = glfwCreateWindow(windowWidth, windowHeight, "FreyaLab Mixed-Signal Engineering Native GL Dashboard", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        return;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "Failed to initialize GLAD pointer loader!" << std::endl;
        glfwDestroyWindow(window);
        glfwTerminate();
        return;
    }
    
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    
    ImGui::StyleColorsDark();
    
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330 core");
    
    io.Fonts->AddFontDefault();

    GLuint shaderProgram = compileShaders();
    GLint projLoc = glGetUniformLocation(shaderProgram, "projection");
    
    GLuint globalVAO, globalVBO;
    glGenVertexArrays(1, &globalVAO);
    glGenBuffers(1, &globalVBO);
    
    synchronizeVisualsFromNetlist();
    executeActiveSimulation();
    
    char textBuffer[8192];
    strncpy(textBuffer, rawNetlistBuffer.c_str(), sizeof(textBuffer));
    
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        
        glfwGetFramebufferSize(window, &windowWidth, &windowHeight);
        
        // --- 1. ImGui v1.93 Control Panel ---
        ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2(400, (float)windowHeight), ImGuiCond_Always);
        ImGui::Begin("FreyaLab Netlist Engineering Engine", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize);
        
        ImGui::SeparatorText("Simulation Settings");
        if (ImGui::RadioButton("Transient Analysis", &analysisTypeActive, 1)) {
            executeActiveSimulation();
        }
        ImGui::SameLine();
        if (ImGui::RadioButton("AC Sweep (Bode)", &analysisTypeActive, 2)) {
            executeActiveSimulation();
        }
        
        ImGui::SeparatorText("Active SPICE Netlist");
        if (ImGui::InputTextMultiline("##netlist_editor", textBuffer, IM_ARRAYSIZE(textBuffer), ImVec2(-1, 350), ImGuiInputTextFlags_AllowTabInput)) {
            rawNetlistBuffer = textBuffer;
        }
        
        if (ImGui::Button("Compile and Parse Geometry", ImVec2(-1, 40))) {
            synchronizeVisualsFromNetlist();
            executeActiveSimulation();
        }
        
        ImGui::SeparatorText("Identified Active Matrix Nodes");
        for (const auto& pair : visualComponents) {
            ImGui::BulletText("%s (%s) - Val: %.3e", pair.second.name.c_str(), pair.second.typeStr.c_str(), pair.second.scalarValue);
        }
        
        ImGui::End();
        
        // --- 2. ImPlot v1.1 Layout Plotting ---
        ImGui::SetNextWindowPos(ImVec2(400, (float)windowHeight * 0.5f), ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2((float)windowWidth - 400, (float)windowHeight * 0.5f), ImGuiCond_Always);
        ImGui::Begin("Interactive Signal Plot Analytics", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize);
        
        if (analysisTypeActive == 1) {
            if (ImPlot::BeginPlot("Transient Waveform Outputs", ImVec2(-1, -1))) {
                ImPlot::SetupAxis(ImAxis_X1, "Time (us)");
                ImPlot::SetupAxis(ImAxis_Y1, "Amplitude (V)");
                if (!plotTimeSteps.empty()) {
                    ImPlot::PlotLine("Input Waveform", plotTimeSteps.data(), plotInputWaveform.data(), (int)plotTimeSteps.size());
                    ImPlot::PlotLine("Output Response", plotTimeSteps.data(), plotOutputWaveform.data(), (int)plotTimeSteps.size());
                }
                ImPlot::EndPlot();
            }
        } 
        else if (analysisTypeActive == 2) {
            if (ImPlot::BeginPlot("Bode Frequency Sweep Magnitude", ImVec2(-1, -1))) {
                ImPlot::SetupAxis(ImAxis_X1, "Frequency Log10(Hz)");
                ImPlot::SetupAxis(ImAxis_Y1, "Magnitude Response (dB)");
                if (!plotFrequencies.empty()) {
                    ImPlot::PlotLine("V(out) Response", plotFrequencies.data(), plotBodeMagnitudeDb.data(), (int)plotFrequencies.size());
                }
                ImPlot::EndPlot();
            }
        }
        ImGui::End();

        // --- 3. Render Native Layout Framework ---
        int schematicWidth = windowWidth - 400;
        int schematicHeight = windowHeight / 2;
        glViewport(400, schematicHeight, schematicWidth, schematicHeight);
        
        glEnable(GL_SCISSOR_TEST);
        glScissor(400, schematicHeight, schematicWidth, schematicHeight);
        glClearColor(0.12f, 0.12f, 0.15f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        glDisable(GL_SCISSOR_TEST);
        
        glUseProgram(shaderProgram);
        float orthoMatrix[16] = {
            2.0f / (float)schematicWidth, 0.0f,                          0.0f,  0.0f,
            0.0f,                        -2.0f / (float)schematicHeight, 0.0f,  0.0f,
            0.0f,                         0.0f,                         -1.0f,  0.0f,
           -1.0f,                         1.0f,                          0.0f,  1.0f
        };
        glUniformMatrix4fv(projLoc, 1, GL_FALSE, orthoMatrix);
        
        extern void drawVisualSchematicGL(const std::map<std::string, VisualComponent>& components, GLuint vao, GLuint vbo);
        drawVisualSchematicGL(visualComponents, globalVAO, globalVBO);
        
        glViewport(0, 0, windowWidth, windowHeight);
        
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        
        glfwSwapBuffers(window);
    }
    
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

void drawVisualSchematicGL(const std::map<std::string, VisualComponent>& components, GLuint vao, GLuint vbo) {
    if(components.empty()) return;
    
    std::vector<Vertex2D> vertices;
    vertices.reserve(components.size() * 6);
    
    for (const auto& pair : components) {
        const auto& vc = pair.second;
        float r = 0.2f, g = 0.6f, b = 0.8f;
        if(vc.typeStr == "V_Source") {
            r = 0.8f; g = 0.4f; b = 0.2f;
        } else if(vc.typeStr == "Capacitor") {
            r = 0.2f; g = 0.8f; b = 0.4f;
        }
        
        float w = 120.0f, h = 60.0f;
        vertices.push_back({vc.posX,     vc.posY,     r, g, b});
        vertices.push_back({vc.posX + w, vc.posY,     r, g, b});
        vertices.push_back({vc.posX + w, vc.posY + h, r, g, b});
        vertices.push_back({vc.posX,     vc.posY,     r, g, b});
        vertices.push_back({vc.posX + w, vc.posY + h, r, g, b});
        vertices.push_back({vc.posX,     vc.posY + h, r, g, b});
    }
    
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex2D), vertices.data(), GL_STREAM_DRAW);
    
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex2D), (void*)offsetof(Vertex2D, x));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex2D), (void*)offsetof(Vertex2D, r));
    glEnableVertexAttribArray(1);
    
    glDrawArrays(GL_TRIANGLES, 0, (GLsizei)vertices.size());
}
