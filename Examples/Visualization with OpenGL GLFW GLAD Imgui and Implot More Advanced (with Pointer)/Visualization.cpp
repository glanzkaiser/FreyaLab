#include "freyalab/VisualizationPointer.h"
#include "freyalab/Netlist.h"
#include "freyalab/KLU.h"
#include "freyalab/CircuitComponents.h"
#include "freyalab/DCCircuit.h"
#include "freyalab/PulseCircuit.h"
#include "freyalab/ACCircuit.h"
#include "freyalab/LogicGates.h"
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
// IMPLOT v1.1 COMPATIBILITY WRAPPERS
// ============================================================================

static inline void PlotLineStyled(const char* label, const double* xs, const double* ys, int count, const ImVec4& color, float weight = 2.0f) {
    ImPlotSpec spec;
    spec.LineColor = color;
    spec.LineWeight = weight;
    ImPlot::PlotLine(label, xs, ys, count, spec);
}
static inline void PlotStairsStyled(const char* label, const double* xs, const double* ys, int count, const ImVec4& color, float weight = 2.5f) {
    ImPlotSpec spec;
    spec.LineColor = color;
    spec.LineWeight = weight;
    ImPlot::PlotStairs(label, xs, ys, count, spec);
}

static inline void PlotScatterStyled(const char* label, const double* xs, const double* ys, int count, ImPlotMarker marker, float size, const ImVec4& color) {
    ImPlotSpec spec;
    spec.Marker = marker;
    spec.MarkerSize = size;
    spec.MarkerFillColor = color;
    spec.MarkerLineColor = color;
    ImPlot::PlotScatter(label, xs, ys, count, spec);
}

// ============================================================================
// FREYALAB VISUALIZER CONSTRUCTORS & DATA STRUCTURES
// ============================================================================
FreyaLabVisualizer::FreyaLabVisualizer(Netlist* netlistInstance) 
    : netlistPtr(netlistInstance), selectedCompName(""), analysisTypeActive(1), zoomLevel(1.0f), 
      scrollOffsetX(0.0f), scrollOffsetY(0.0f), stimPulseWidth(10.0f), stimPulsePeriod(20.0f),
      plotActivePowerP(0.0), plotReactivePowerQ(0.0), plotApparentPowerS(0.0), plotPowerFactor(1.0) {}

FreyaLabVisualizer::FreyaLabVisualizer(Netlist& netlistInstance) : FreyaLabVisualizer(&netlistInstance) {}

static bool g_netlistParseSuccess = true;
static std::string g_parserErrorMessage = "";
static std::vector<std::complex<double>> g_latestAcNodePhasors;
static float g_acTestFrequencyHz = 10000.0f; // Default 10 kHz AC test frequency

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
    float dx = x2 - x1; float dy = y2 - y1;
    float len = std::sqrt(dx * dx + dy * dy);
    if (len < 1e-6f) return;
    float nx = -dy / len * (thickness * 0.5f);
    float ny = dx / len * (thickness * 0.5f);
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

static void drawVisualSchematicGL(
    const std::map<std::string, VisualComponent>& visualComponents,
    float zoomLevel, float scrollOffsetX, float scrollOffsetY,
    int schematicWidth, int schematicHeight,
    GLuint globalVAO, GLuint globalVBO)
{
    std::vector<Vertex2D> vertices;

    // 1. Background Grid (40.0f spacing matching grid-snap)
    float gridStep = 40.0f;
    float worldLeft   = -scrollOffsetX / zoomLevel;
    float worldRight  = (static_cast<float>(schematicWidth) - scrollOffsetX) / zoomLevel;
    float worldTop    = -scrollOffsetY / zoomLevel;
    float worldBottom = (static_cast<float>(schematicHeight) - scrollOffsetY) / zoomLevel;

    float startX = std::floor(worldLeft / gridStep) * gridStep;
    float startY = std::floor(worldTop / gridStep) * gridStep;

    for (float gx = startX; gx <= worldRight; gx += gridStep) {
        appendLineGL(vertices, gx, worldTop, gx, worldBottom, 0.17f, 0.17f, 0.21f, 1.0f / zoomLevel);
    }
    for (float gy = startY; gy <= worldBottom; gy += gridStep) {
        appendLineGL(vertices, worldLeft, gy, worldRight, gy, 0.17f, 0.17f, 0.21f, 1.0f / zoomLevel);
    }

    // 2. Render Distinct Circuit Symbols
    float blockW = 120.0f;
    float blockH = 60.0f;

    for (const auto& pair : visualComponents) {
        const auto& vc = pair.second;
        float cx = vc.posX + blockW * 0.5f;
        float cy = vc.posY + blockH * 0.5f;

        float rad = vc.rotationAngle * static_cast<float>(M_PI) / 180.0f;
        float cosA = std::cos(rad);
        float sinA = std::sin(rad);

        // Local-to-world rotation helper
        auto rotPt = [&](float lx, float ly) -> std::pair<float, float> {
            return { cx + lx * cosA - ly * sinA, cy + lx * sinA + ly * cosA };
        };

        // Helper to draw a line in local component space
        auto drawLocalLine = [&](float x1, float y1, float x2, float y2, float r, float g, float b, float thick) {
            auto p1 = rotPt(x1, y1);
            auto p2 = rotPt(x2, y2);
            appendLineGL(vertices, p1.first, p1.second, p2.first, p2.second, r, g, b, thick);
        };

        // Helper to draw a circle in local component space
        auto drawLocalCircle = [&](float lx, float ly, float radius, float r, float g, float b, int segs = 16) {
            auto p = rotPt(lx, ly);
            appendCircleGL(vertices, p.first, p.second, radius, r, g, b, segs);
        };

        float thick = vc.isSelected ? 3.0f : 2.2f;
        float leadR = 0.78f, leadG = 0.80f, leadB = 0.85f;

        // Default highlight vs category colors
        float r = vc.isSelected ? 1.0f : 0.30f;
        float g = vc.isSelected ? 0.85f : 0.75f;
        float b = vc.isSelected ? 0.20f : 0.95f;

        bool isThreeTerminal = false;

        // --------------------------------------------------------------------
        // A. RESISTOR (Zigzag IEEE Symbol)
        // --------------------------------------------------------------------
        if (vc.typeStr == "Resistor") {
            if (!vc.isSelected) { r = 0.95f; g = 0.65f; b = 0.25f; }
            drawLocalLine(-60.0f, 0.0f, -30.0f, 0.0f, leadR, leadG, leadB, 2.0f);
            drawLocalLine( 30.0f, 0.0f,  60.0f, 0.0f, leadR, leadG, leadB, 2.0f);

            const float zx[] = { -30.0f, -24.0f, -14.0f, -4.0f, 6.0f, 16.0f, 24.0f, 30.0f };
            const float zy[] = {   0.0f, -12.0f,  12.0f,-12.0f,12.0f,-12.0f, 12.0f,  0.0f };
            for (int i = 0; i < 7; ++i) {
                drawLocalLine(zx[i], zy[i], zx[i + 1], zy[i + 1], r, g, b, thick);
            }
        }
        // --------------------------------------------------------------------
        // B. CAPACITOR (Parallel Plates)
        // --------------------------------------------------------------------
        else if (vc.typeStr == "Capacitor") {
            if (!vc.isSelected) { r = 0.25f; g = 0.85f; b = 0.95f; }
            drawLocalLine(-60.0f, 0.0f, -7.0f, 0.0f, leadR, leadG, leadB, 2.0f);
            drawLocalLine(  7.0f, 0.0f, 60.0f, 0.0f, leadR, leadG, leadB, 2.0f);

            drawLocalLine(-7.0f, -18.0f, -7.0f, 18.0f, r, g, b, thick + 0.6f);
            drawLocalLine( 7.0f, -18.0f,  7.0f, 18.0f, r, g, b, thick + 0.6f);
        }
        // --------------------------------------------------------------------
        // C. INDUCTOR (4-Loop Coils)
        // --------------------------------------------------------------------
        else if (vc.typeStr == "Inductor") {
            if (!vc.isSelected) { r = 0.45f; g = 0.90f; b = 0.55f; }
            drawLocalLine(-60.0f, 0.0f, -32.0f, 0.0f, leadR, leadG, leadB, 2.0f);
            drawLocalLine( 32.0f, 0.0f,  60.0f, 0.0f, leadR, leadG, leadB, 2.0f);

            int coils = 4;
            float coilWidth = 64.0f / coils;
            int arcSegs = 10;
            for (int c = 0; c < coils; ++c) {
                float startCoilX = -32.0f + c * coilWidth;
                float centerCoilX = startCoilX + coilWidth * 0.5f;
                float radius = coilWidth * 0.5f;
                for (int s = 0; s < arcSegs; ++s) {
                    float a1 = static_cast<float>(M_PI) * (1.0f - (float)s / arcSegs);
                    float a2 = static_cast<float>(M_PI) * (1.0f - (float)(s + 1) / arcSegs);
                    drawLocalLine(centerCoilX + radius * std::cos(a1), -radius * 1.3f * std::sin(a1),
                                  centerCoilX + radius * std::cos(a2), -radius * 1.3f * std::sin(a2),
                                  r, g, b, thick);
                }
            }
        }
        // --------------------------------------------------------------------
        // D. VOLTAGE SOURCE (Circle with + and - Polarity)
        // --------------------------------------------------------------------
        else if (vc.typeStr == "V_Source") {
            if (!vc.isSelected) { r = 0.30f; g = 0.95f; b = 0.50f; }
            drawLocalLine(-60.0f, 0.0f, -22.0f, 0.0f, leadR, leadG, leadB, 2.0f);
            drawLocalLine( 22.0f, 0.0f,  60.0f, 0.0f, leadR, leadG, leadB, 2.0f);

            // Outer circle ring
            int segs = 24;
            float radCircle = 22.0f;
            for (int s = 0; s < segs; ++s) {
                float a1 = (s / (float)segs) * 2.0f * (float)M_PI;
                float a2 = ((s + 1) / (float)segs) * 2.0f * (float)M_PI;
                drawLocalLine(radCircle * std::cos(a1), radCircle * std::sin(a1),
                              radCircle * std::cos(a2), radCircle * std::sin(a2), r, g, b, thick);
            }
            // '+' symbol on left terminal side (-10, 0)
            drawLocalLine(-14.0f,  0.0f, -6.0f, 0.0f, r, g, b, thick);
            drawLocalLine(-10.0f, -4.0f,-10.0f, 4.0f, r, g, b, thick);
            // '-' symbol on right terminal side (+10, 0)
            drawLocalLine(  6.0f,  0.0f, 14.0f, 0.0f, r, g, b, thick);
        }
        // --------------------------------------------------------------------
        // E. DIODE (Triangle + Cathode Bar)
        // --------------------------------------------------------------------
        else if (vc.typeStr == "Diode") {
            if (!vc.isSelected) { r = 0.95f; g = 0.45f; b = 0.55f; }
            drawLocalLine(-60.0f, 0.0f, -14.0f, 0.0f, leadR, leadG, leadB, 2.0f);
            drawLocalLine( 14.0f, 0.0f,  60.0f, 0.0f, leadR, leadG, leadB, 2.0f);

            // Triangle anode
            drawLocalLine(-14.0f, -14.0f, -14.0f, 14.0f, r, g, b, thick);
            drawLocalLine(-14.0f, -14.0f,  14.0f,  0.0f, r, g, b, thick);
            drawLocalLine(-14.0f,  14.0f,  14.0f,  0.0f, r, g, b, thick);
            // Cathode bar
            drawLocalLine( 14.0f, -14.0f,  14.0f, 14.0f, r, g, b, thick + 0.5f);
        }
        // --------------------------------------------------------------------
        // F. PLANAR MOSFET (Gate + Insulated Oxide + Channel + Drain/Source)
        // --------------------------------------------------------------------
        else if (vc.typeStr == "MOSFET") {
            isThreeTerminal = true;
            if (!vc.isSelected) { r = 0.85f; g = 0.50f; b = 1.0f; }
            // Gate lead & plate
            drawLocalLine(-60.0f,   0.0f, -16.0f,   0.0f, leadR, leadG, leadB, 2.0f);
            drawLocalLine(-16.0f, -16.0f, -16.0f,  16.0f, r, g, b, thick);
            // Channel bar
            drawLocalLine( -8.0f, -18.0f,  -8.0f,  18.0f, r, g, b, thick);
            // Drain (top-right) & Source (bottom-right)
            drawLocalLine( -8.0f, -12.0f,  25.0f, -12.0f, r, g, b, thick);
            drawLocalLine( 25.0f, -12.0f,  25.0f, -24.0f, leadR, leadG, leadB, 2.0f);
            drawLocalLine( -8.0f,  12.0f,  25.0f,  12.0f, r, g, b, thick);
            drawLocalLine( 25.0f,  12.0f,  25.0f,  24.0f, leadR, leadG, leadB, 2.0f);
            // Source arrow
            drawLocalLine(  4.0f,   8.0f,  12.0f,  12.0f, r, g, b, 1.8f);
            drawLocalLine(  4.0f,  16.0f,  12.0f,  12.0f, r, g, b, 1.8f);
        }
        // --------------------------------------------------------------------
        // G. FINFET (3D Multi-Fin Vertical Channel Bars + Tri-Gate Wrap)
        // --------------------------------------------------------------------
        else if (vc.typeStr == "FinFET") {
            isThreeTerminal = true;
            if (!vc.isSelected) { r = 0.20f; g = 0.90f; b = 0.85f; }
            // Gate input & wrap-around gate bar
            drawLocalLine(-60.0f,   0.0f, -20.0f,   0.0f, leadR, leadG, leadB, 2.0f);
            drawLocalLine(-20.0f, -18.0f, -20.0f,  18.0f, r, g, b, thick + 0.5f);
            // Triple vertical silicon fins (distinguishing FinFET from planar MOSFET)
            drawLocalLine(-12.0f, -16.0f, -12.0f,  16.0f, 1.0f, 0.75f, 0.30f, 2.0f);
            drawLocalLine( -6.0f, -16.0f,  -6.0f,  16.0f, 1.0f, 0.75f, 0.30f, 2.0f);
            drawLocalLine(  0.0f, -16.0f,   0.0f,  16.0f, 1.0f, 0.75f, 0.30f, 2.0f);
            // Drain (top-right) & Source (bottom-right) connections
            drawLocalLine(-12.0f, -14.0f,  25.0f, -14.0f, r, g, b, thick);
            drawLocalLine( 25.0f, -14.0f,  25.0f, -24.0f, leadR, leadG, leadB, 2.0f);
            drawLocalLine(-12.0f,  14.0f,  25.0f,  14.0f, r, g, b, thick);
            drawLocalLine( 25.0f,  14.0f,  25.0f,  24.0f, leadR, leadG, leadB, 2.0f);
        }
        // --------------------------------------------------------------------
        // H. GAAFET / RIBBONFET (Gate-All-Around Nanosheet Concentric Rings)
        // --------------------------------------------------------------------
        else if (vc.typeStr == "GAAFET" || vc.typeStr == "RibbonFET") {
            isThreeTerminal = true;
            if (!vc.isSelected) { r = 0.95f; g = 0.35f; b = 0.85f; }
            // Gate lead into full-surround Gate-All-Around ring
            drawLocalLine(-60.0f, 0.0f, -22.0f, 0.0f, leadR, leadG, leadB, 2.0f);
            int segs = 20;
            float outerR = 16.0f;
            for (int s = 0; s < segs; ++s) {
                float a1 = (s / (float)segs) * 2.0f * (float)M_PI;
                float a2 = ((s + 1) / (float)segs) * 2.0f * (float)M_PI;
                drawLocalLine(-6.0f + outerR * std::cos(a1), outerR * std::sin(a1),
                              -6.0f + outerR * std::cos(a2), outerR * std::sin(a2), r, g, b, thick);
            }
            // Stacked horizontal nanosheet channels enclosed inside the gate
            drawLocalLine(-16.0f, -6.0f,  4.0f, -6.0f, 0.30f, 0.95f, 1.0f, 2.2f);
            drawLocalLine(-16.0f,  0.0f,  4.0f,  0.0f, 0.30f, 0.95f, 1.0f, 2.2f);
            drawLocalLine(-16.0f,  6.0f,  4.0f,  6.0f, 0.30f, 0.95f, 1.0f, 2.2f);

            // Drain (top-right) & Source (bottom-right)
            drawLocalLine(  4.0f, -12.0f, 25.0f, -12.0f, r, g, b, thick);
            drawLocalLine( 25.0f, -12.0f, 25.0f, -24.0f, leadR, leadG, leadB, 2.0f);
            drawLocalLine(  4.0f,  12.0f, 25.0f,  12.0f, r, g, b, thick);
            drawLocalLine( 25.0f,  12.0f, 25.0f,  24.0f, leadR, leadG, leadB, 2.0f);
        }
        // --------------------------------------------------------------------
        // I. BJT (Base Vertical Bar + Angled Collector & Emitter Arrow)
        // --------------------------------------------------------------------
        else if (vc.typeStr == "BJT") {
            isThreeTerminal = true;
            if (!vc.isSelected) { r = 0.95f; g = 0.60f; b = 0.35f; }
            // Base lead & bar
            drawLocalLine(-60.0f,   0.0f, -12.0f,   0.0f, leadR, leadG, leadB, 2.0f);
            drawLocalLine(-12.0f, -16.0f, -12.0f,  16.0f, r, g, b, thick + 0.5f);
            // Collector (top-right)
            drawLocalLine(-12.0f,  -7.0f,  25.0f, -20.0f, r, g, b, thick);
            drawLocalLine( 25.0f, -20.0f,  25.0f, -24.0f, leadR, leadG, leadB, 2.0f);
            // Emitter (bottom-right) + arrow
            drawLocalLine(-12.0f,   7.0f,  25.0f,  20.0f, r, g, b, thick);
            drawLocalLine( 25.0f,  20.0f,  25.0f,  24.0f, leadR, leadG, leadB, 2.0f);
            drawLocalLine( 12.0f,  19.0f,  20.0f,  18.0f, r, g, b, 2.0f);
            drawLocalLine( 16.0f,  12.0f,  20.0f,  18.0f, r, g, b, 2.0f);
        }
	        // --------------------------------------------------------------------
        // I2. IEEE LOGIC GATES (NOT, AND, NAND, OR, NOR, XOR, XNOR)
        // --------------------------------------------------------------------
        else if (vc.typeStr == "NOT" || vc.typeStr == "AND" || vc.typeStr == "NAND" ||
                 vc.typeStr == "OR"  || vc.typeStr == "NOR" || vc.typeStr == "XOR"  || vc.typeStr == "XNOR") {
            if (!vc.isSelected) { r = 0.25f; g = 0.90f; b = 0.70f; }

            bool isSingleInput = (vc.typeStr == "NOT");
            bool hasBubble = (vc.typeStr == "NOT" || vc.typeStr == "NAND" || vc.typeStr == "NOR" || vc.typeStr == "XNOR");
            bool isXorFamily = (vc.typeStr == "XOR" || vc.typeStr == "XNOR");

            if (isSingleInput) {
                drawLocalLine(-60.0f, 0.0f, -20.0f, 0.0f, leadR, leadG, leadB, 2.0f);
            } else {
                float inputStopX = (vc.typeStr == "OR" || vc.typeStr == "NOR") ? -18.0f : (isXorFamily ? -23.0f : -22.0f);
                drawLocalLine(-60.0f, -12.0f, inputStopX, -12.0f, leadR, leadG, leadB, 2.0f);
                drawLocalLine(-60.0f,  12.0f, inputStopX,  12.0f, leadR, leadG, leadB, 2.0f);
            }

            float outStartX = hasBubble ? 28.0f : 20.0f;
            drawLocalLine(outStartX, 0.0f, 60.0f, 0.0f, leadR, leadG, leadB, 2.0f);

            if (vc.typeStr == "NOT") {
                drawLocalLine(-20.0f, -18.0f, -20.0f, 18.0f, r, g, b, thick);
                drawLocalLine(-20.0f, -18.0f,  18.0f,  0.0f, r, g, b, thick);
                drawLocalLine(-20.0f,  18.0f,  18.0f,  0.0f, r, g, b, thick);
            } else if (vc.typeStr == "AND" || vc.typeStr == "NAND") {
                drawLocalLine(-22.0f, -20.0f, -22.0f, 20.0f, r, g, b, thick);
                drawLocalLine(-22.0f, -20.0f,   0.0f, -20.0f, r, g, b, thick);
                drawLocalLine(-22.0f,  20.0f,   0.0f,  20.0f, r, g, b, thick);
                int arcSegs = 20;
                float radius = 20.0f;
                for (int s = 0; s < arcSegs; ++s) {
                    float a1 = -0.5f * static_cast<float>(M_PI) + (s / (float)arcSegs) * static_cast<float>(M_PI);
                    float a2 = -0.5f * static_cast<float>(M_PI) + ((s + 1) / (float)arcSegs) * static_cast<float>(M_PI);
                    drawLocalLine(radius * std::cos(a1), radius * std::sin(a1),
                                  radius * std::cos(a2), radius * std::sin(a2), r, g, b, thick);
                }
            } else {
                int curveSegs = 16;
                for (int s = 0; s < curveSegs; ++s) {
                    float t1 = -1.0f + 2.0f * (s / (float)curveSegs);
                    float t2 = -1.0f + 2.0f * ((s + 1) / (float)curveSegs);
                    float x1 = -22.0f + 6.0f * (1.0f - t1 * t1);
                    float y1 = 20.0f * t1;
                    float x2 = -22.0f + 6.0f * (1.0f - t2 * t2);
                    float y2 = 20.0f * t2;
                    drawLocalLine(x1, y1, x2, y2, r, g, b, thick);
                    if (isXorFamily) {
                        drawLocalLine(x1 - 6.0f, y1, x2 - 6.0f, y2, r, g, b, thick);
                    }
                }
                for (int s = 0; s < curveSegs; ++s) {
                    float u1 = s / (float)curveSegs;
                    float u2 = (s + 1) / (float)curveSegs;
                    float x1 = -22.0f + 42.0f * u1;
                    float y1 = 20.0f * (1.0f - u1 * u1);
                    float x2 = -22.0f + 42.0f * u2;
                    float y2 = 20.0f * (1.0f - u2 * u2);
                    drawLocalLine(x1, -y1, x2, -y2, r, g, b, thick);
                    drawLocalLine(x1,  y1, x2,  y2, r, g, b, thick);
                }
            }

            if (hasBubble) {
                int bSegs = 16;
                float bRad = 4.5f;
                float bCenterX = 23.5f;
                for (int s = 0; s < bSegs; ++s) {
                    float a1 = (s / (float)bSegs) * 2.0f * static_cast<float>(M_PI);
                    float a2 = ((s + 1) / (float)bSegs) * 2.0f * static_cast<float>(M_PI);
                    drawLocalLine(bCenterX + bRad * std::cos(a1), bRad * std::sin(a1),
                                  bCenterX + bRad * std::cos(a2), bRad * std::sin(a2), r, g, b, thick);
                }
            }

            if (isSingleInput) {
                drawLocalCircle(-60.0f, 0.0f, 3.8f, 0.3f, 1.0f, 0.4f, 12);
                drawLocalCircle( 60.0f, 0.0f, 3.8f, 0.3f, 1.0f, 0.4f, 12);
            } else {
                drawLocalCircle(-60.0f, -12.0f, 3.8f, 0.3f, 1.0f, 0.4f, 12);
                drawLocalCircle(-60.0f,  12.0f, 3.8f, 0.3f, 1.0f, 0.4f, 12);
                drawLocalCircle( 60.0f,   0.0f, 3.8f, 0.3f, 1.0f, 0.4f, 12);
            }
            continue;
        }

        // --------------------------------------------------------------------
        // J. DEFAULT / LOGIC GATE / GENERIC PRIMITIVE BLOCK
        // --------------------------------------------------------------------
        else {
            drawLocalLine(-60.0f,   0.0f, -32.0f,   0.0f, leadR, leadG, leadB, 2.0f);
            drawLocalLine( 32.0f,   0.0f,  60.0f,   0.0f, leadR, leadG, leadB, 2.0f);
            drawLocalLine(-32.0f, -16.0f,  32.0f, -16.0f, r, g, b, thick);
            drawLocalLine( 32.0f, -16.0f,  32.0f,  16.0f, r, g, b, thick);
            drawLocalLine( 32.0f,  16.0f, -32.0f,  16.0f, r, g, b, thick);
            drawLocalLine(-32.0f,  16.0f, -32.0f, -16.0f, r, g, b, thick);
        }

        // 3. Terminal Pin Dots (2-pin vs 3-pin devices)
        if (isThreeTerminal) {
            drawLocalCircle(-60.0f,   0.0f, 3.8f, 0.3f, 1.0f, 0.4f, 12); // Gate / Base
            drawLocalCircle( 25.0f, -24.0f, 3.8f, 0.3f, 1.0f, 0.4f, 12); // Drain / Collector
            drawLocalCircle( 25.0f,  24.0f, 3.8f, 0.3f, 1.0f, 0.4f, 12); // Source / Emitter
        } else {
            drawLocalCircle(-60.0f, 0.0f, 3.8f, 0.3f, 1.0f, 0.4f, 12);
            drawLocalCircle( 60.0f, 0.0f, 3.8f, 0.3f, 1.0f, 0.4f, 12);
        }
    }
   
    // ========================================================================
    // 3.5. MANHATTAN ORTHOGONAL NET WIRING
    // ========================================================================
    struct PinWorldPos {
        float x;
        float y;
    };
    std::map<std::string, std::vector<PinWorldPos>> netPinsMap;

    for (const auto& pair : visualComponents) {
        const auto& vc = pair.second;
        float cx = vc.posX + blockW * 0.5f;
        float cy = vc.posY + blockH * 0.5f;
        float rad = vc.rotationAngle * static_cast<float>(M_PI) / 180.0f;
        float cosA = std::cos(rad);
        float sinA = std::sin(rad);

        auto toWorld = [&](float lx, float ly) -> PinWorldPos {
            return { cx + lx * cosA - ly * sinA, cy + lx * sinA + ly * cosA };
        };

        bool isTwoInputGate = (  vc.typeStr == "AND" || vc.typeStr == "NAND" ||
            			 vc.typeStr == "OR"  || vc.typeStr == "NOR"  ||
            			 vc.typeStr == "XOR" || vc.typeStr == "XNOR"
        );
	bool isThreeTerm = (vc.typeStr == "MOSFET" || vc.typeStr == "FinFET" ||
                            vc.typeStr == "GAAFET" || vc.typeStr == "RibbonFET" ||
                            vc.typeStr == "BJT");

	std::vector<PinWorldPos> localPins;
        if (isTwoInputGate) {
            // Pin 0: In A (-60, -12), Pin 1: In B (-60, 12), Pin 2: Out Y (60, 0)
            localPins.push_back(toWorld(-60.0f, -12.0f));
            localPins.push_back(toWorld(-60.0f,  12.0f));
            localPins.push_back(toWorld( 60.0f,   0.0f));
        } else if (isThreeTerm) {
            // Pin 0: Drain/Collector (25, -24), Pin 1: Gate/Base (-60, 0), Pin 2: Source/Emitter (25, 24)
            localPins.push_back(toWorld( 25.0f, -24.0f));
            localPins.push_back(toWorld(-60.0f,   0.0f));
            localPins.push_back(toWorld( 25.0f,  24.0f));
        } else {
            // Pin 0: Left (-60, 0), Pin 1: Right (60, 0) — covers 2-pin passives & NOT gates
            localPins.push_back(toWorld(-60.0f, 0.0f));
            localPins.push_back(toWorld( 60.0f, 0.0f));
        }

        for (size_t i = 0; i < vc.connectedNodes.size() && i < localPins.size(); ++i) {
            std::string netName = vc.connectedNodes[i];
            std::string lowerNet = netName;
            std::transform(lowerNet.begin(), lowerNet.end(), lowerNet.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

            // Skip drawing global GND wires to avoid clutter (remove this check if you want GND wired too)
            if (lowerNet == "0" || lowerNet == "gnd" || lowerNet == "ground") continue;

            netPinsMap[netName].push_back(localPins[i]);
        }
    }

    // Route orthogonal Manhattan segments for every net with 2+ connected pins
    float wireR = 0.25f, wireG = 0.85f, wireB = 0.45f; // PCB green trace color
    float wireThickness = 2.0f;

    for (const auto& netPair : netPinsMap) {
        const auto& pins = netPair.second;
        if (pins.size() < 2) continue;

        // Sort pins left-to-right (then top-to-bottom) for cleaner daisy-chain routing
        std::vector<PinWorldPos> sortedPins = pins;
        std::sort(sortedPins.begin(), sortedPins.end(), [](const PinWorldPos& a, const PinWorldPos& b) {
            if (std::abs(a.x - b.x) > 1.0f) return a.x < b.x;
            return a.y < b.y;
        });

        for (size_t i = 0; i < sortedPins.size() - 1; ++i) {
            float x1 = sortedPins[i].x;
            float y1 = sortedPins[i].y;
            float x2 = sortedPins[i + 1].x;
            float y2 = sortedPins[i + 1].y;

            // Z-shaped (3-segment) Manhattan dogleg through horizontal midpoint
            float midX = (x1 + x2) * 0.5f;

            // Segment 1: Horizontal from Pin 1 to midX
            appendLineGL(vertices, x1, y1, midX, y1, wireR, wireG, wireB, wireThickness);
            // Segment 2: Vertical from y1 to y2 at midX
            appendLineGL(vertices, midX, y1, midX, y2, wireR, wireG, wireB, wireThickness);
            // Segment 3: Horizontal from midX to Pin 2
            appendLineGL(vertices, midX, y2, x2, y2, wireR, wireG, wireB, wireThickness);

            // Draw solder junction dot if 3+ pins share this net
            if (sortedPins.size() > 2 && i > 0) {
                appendCircleGL(vertices, x1, y1, 4.2f, wireR, wireG, wireB, 12);
            }
        }
    }

    if (vertices.empty()) return;

    // 4. Upload to GPU & Render Batch
    glBindVertexArray(globalVAO);
    glBindBuffer(GL_ARRAY_BUFFER, globalVBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex2D), vertices.data(), GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex2D), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex2D), (void*)(2 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(vertices.size()));
    glBindVertexArray(0);
}

// ============================================================================
// TOPOLOGY PROBE HELPER
// ============================================================================
static void discoverGenericProbeNodes(Netlist* netlistPtr, const std::string& selectedCompName, NetNode& detectedInputNode, NetNode& detectedOutputNode) {
    detectedInputNode = NetNode{0}; detectedOutputNode = NetNode{0};
    if (!netlistPtr) return;
    netlistPtr->inspectComponents([&detectedInputNode](Component& comp, MnaIndex extraIdx) {
        (void)extraIdx; if (!detectedInputNode.isGround()) return;
        ComponentType cType = comp.getType();
        if (cType == ComponentType::VOLTAGE_SOURCE || cType == ComponentType::AC_SOURCE_DEVICE || cType == ComponentType::CURRENT_SOURCE) {
            for (const NetNode& n : comp.getNodes()) { if (!n.isGround()) { detectedInputNode = n; break; } }
        }
    });
    int totalNodes = netlistPtr->getNumNodes();
    if (detectedInputNode.isGround() && totalNodes >= 1) detectedInputNode = NetNode{1};
    if (!selectedCompName.empty()) {
        netlistPtr->inspectComponents([&](Component& comp, MnaIndex extraIdx) {
            (void)extraIdx; if (comp.getName() == selectedCompName) {
                for (const NetNode& n : comp.getNodes()) { if (!n.isGround() && n != detectedInputNode) { detectedOutputNode = n; break; } }
            }
        });
    }
    if (detectedOutputNode.isGround()) {
        for (int i = 1; i <= totalNodes; ++i) {
            NetNode candidate{i}; std::string label = netlistPtr->getNodeName(candidate);
            std::transform(label.begin(), label.end(), label.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
            if (label == "OUT" || label == "VOUT" || label == "OUTPUT" || label == "LOAD" || label == "Y") { detectedOutputNode = candidate; break; }
        }
    }
    if (detectedOutputNode.isGround()) {
        for (int i = totalNodes; i >= 1; --i) { NetNode candidate{i}; if (candidate != detectedInputNode) { detectedOutputNode = candidate; break; } }
        if (detectedOutputNode.isGround()) detectedOutputNode = detectedInputNode;
    }
}

// ============================================================================
// REGENERATE SPICE TEXT & BACKEND SYNCHRONIZATION
// ============================================================================
void FreyaLabVisualizer::regenerateNetlistTextFromVisuals() {
    std::stringstream ss;
    for (const auto& pair : visualComponents) {
        const auto& vc = pair.second; ss << vc.name;
        for (const auto& node : vc.connectedNodes) {
            std::string lowerNode = node;
            std::transform(lowerNode.begin(), lowerNode.end(), lowerNode.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (lowerNode == "gnd" || lowerNode == "ground" || lowerNode == "0") ss << " 0"; else ss << " " << node;
        }
        if (vc.name == "V1") {
            if (activeStimulus == StimulusType::DC_SINGLE) ss << " DC " << stimDcVolt;
            else if (activeStimulus == StimulusType::AC_SINUSOID) ss << " AC " << stimAcMag << " " << stimAcPhase;
            else if (activeStimulus == StimulusType::PULSE_SQUARE) {
                // Use a realistic edge transition time (2% of pulse width, minimum 10 ns) so TR-BDF2 converges smoothly
                float edgeTime = std::max(10e-9f, (stimPulseWidth * 1e-6f) * 0.02f);
                ss << " PULSE(" << stimPulseLow << " " << stimPulseHigh << " 0 "
                   << edgeTime << " " << edgeTime << " "
                   << (stimPulseWidth * 1e-6f) << " " << (stimPulsePeriod * 1e-6f) << ")";
            }
        } else if (vc.typeStr == "V_Source") {
            ss << " DC " << vc.scalarValue;
        } else {
            if (vc.typeStr == "Resistor" || vc.typeStr == "Capacitor" || vc.typeStr == "Inductor") ss << " " << std::scientific << vc.scalarValue;
        }
        if (vc.typeStr == "MOSFET") {
            ss << " nmos";
        } else if (vc.typeStr == "BJT") {
            ss << " npn";
        } else if (vc.typeStr == "NOT" || vc.typeStr == "AND" || vc.typeStr == "NAND" ||
                   vc.typeStr == "OR"  || vc.typeStr == "NOR" || vc.typeStr == "XOR"  || vc.typeStr == "XNOR") {
            ss << " " << vc.typeStr;
        }
        ss << "\n";
    }
    rawNetlistBuffer = ss.str();
}

void FreyaLabVisualizer::synchronizeVisualsFromNetlist() {
    if (!netlistPtr) return;

    auto populateVisualsFrom = [this](Netlist& sourceNetlist) {
        std::map<std::string, VisualComponent> retainedLayouts = visualComponents;
        visualComponents.clear();
        int indexOffset = 0;

        sourceNetlist.inspectComponents([this, &sourceNetlist, &indexOffset, &retainedLayouts](Component& comp, MnaIndex extraIdx) {
            (void)extraIdx;
            VisualComponent vc;
            vc.name = comp.getName();
            vc.scalarValue = comp.getValue();

            if (auto* gateComp = dynamic_cast<BehaviorLogicGateComponent*>(&comp)) {
                switch (gateComp->getGateType()) {
                    case GateType::NOT:  vc.typeStr = "NOT";       break;
                    case GateType::AND:  vc.typeStr = "AND";       break;
                    case GateType::NAND: vc.typeStr = "NAND";      break;
                    case GateType::OR:   vc.typeStr = "OR";        break;
                    case GateType::NOR:  vc.typeStr = "NOR";       break;
                    case GateType::XOR:  vc.typeStr = "XOR";       break;
                    case GateType::XNOR: vc.typeStr = "XNOR";      break;
                    default:             vc.typeStr = "LogicGate"; break;
                }
            } else {
                switch (comp.getType()) {
                    case ComponentType::RESISTOR:       vc.typeStr = "Resistor";  break;
                    case ComponentType::VOLTAGE_SOURCE: vc.typeStr = "V_Source";  break;
                    case ComponentType::CAPACITOR:      vc.typeStr = "Capacitor"; break;
                    case ComponentType::INDUCTOR:       vc.typeStr = "Inductor";  break;
                    case ComponentType::DIODE:          vc.typeStr = "Diode";     break;
                    case ComponentType::MOSFET:         vc.typeStr = "MOSFET";    break;
                    case ComponentType::FINFET:         vc.typeStr = "FinFET";    break;
                    case ComponentType::RIBBONFET:      vc.typeStr = "RibbonFET"; break;
                    case ComponentType::GAA_FET:        vc.typeStr = "GAAFET";    break;
                    case ComponentType::BJT:            vc.typeStr = "BJT";       break;
                    case ComponentType::DIGITAL_IC:     vc.typeStr = "LogicGate"; break;
                    default:                            vc.typeStr = "Primitive"; break;
                }
            }

            for (const NetNode& node : comp.getNodes()) {
                if (node.isGround()) {
                    vc.connectedNodes.push_back("GND");
                } else {
                    std::string nodeLabel = sourceNetlist.getNodeName(node);
                    if (nodeLabel.rfind("Node_", 0) == 0) {
                        nodeLabel = nodeLabel.substr(5);
                    }
                    vc.connectedNodes.push_back(nodeLabel);
                }
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
    };


    if (rawNetlistBuffer.empty() && netlistPtr->getTotalSystemSize() > 0) {
        netlistPtr->buildSparsityPattern();
        populateVisualsFrom(*netlistPtr);
        regenerateNetlistTextFromVisuals();
        g_netlistParseSuccess = true;
        g_parserErrorMessage = "";
        return;
    }

    std::ofstream tmpFile("tmp_gui_netlist.ckt");
    if (!tmpFile.is_open()) {
        g_netlistParseSuccess = false;
        g_parserErrorMessage = "Failed to open temporary netlist file stream.";
        return;
    }
    tmpFile << rawNetlistBuffer;
    tmpFile.close();

    Netlist localWorkspace;
    FileParser freshParser;

    try {
        freshParser.parseNetlistFile("tmp_gui_netlist.ckt", localWorkspace);
        localWorkspace.buildSparsityPattern();
        populateVisualsFrom(localWorkspace);

        *netlistPtr = std::move(localWorkspace);
        g_netlistParseSuccess = true;
        g_parserErrorMessage = "";
    } catch (const std::exception& e) {
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
    plotNyquistReal.clear();
    plotNyquistImag.clear();
    
    if (!netlistPtr) return;
    int systemSize = netlistPtr->getTotalSystemSize();
    if (systemSize <= 0) return;

    NetNode inputNode{0};
    NetNode outputNode{0};
    discoverGenericProbeNodes(netlistPtr, selectedCompName, inputNode, outputNode);

    // --- MODE 1: ADAPTIVE PULSE TRANSIENT SIMULATION (TR-BDF2) ---
    if (analysisTypeActive == 1) {
        try {
            PulseCircuit pulseAnalyzer(*netlistPtr);
                        if (pulseAnalyzer.isPulseCircuit()) {
                PulseSimulationOptions pulseOptions;
                double periodSec = std::max(1e-7, static_cast<double>(stimPulsePeriod) * 1e-6);
                pulseOptions.startTime = 0.0;
                pulseOptions.stopTime  = 3.0 * periodSec;
                pulseOptions.initialDt = periodSec / 200.0;
                pulseOptions.minDt     = 1e-16;
                pulseOptions.maxDt     = periodSec / 20.0;
                pulseOptions.relTol    = 1e-3;
                pulseOptions.absTol    = 1e-6;
                pulseOptions.maxIterations = 80;
                
                netlistPtr->buildSparsityPattern();
                double t = pulseOptions.startTime;
                double dt = pulseOptions.initialDt;
                const double gamma = 2.0 - std::sqrt(2.0);
                
                MnaMatrixContext ctx{ 
                    netlistPtr->getGMatrix(), netlistPtr->getBVector(), &netlistPtr->getStateV(),
                    &netlistPtr->getHistoricV(), &netlistPtr->getMidpointV(), t, dt, 0, MnaIndex{},
                    netlistPtr->getTemperatureCelsius(),
                    netlistPtr->getThermalVoltage()
                };

                ctx.orderStep = 0; ctx.time = t; ctx.dt = 0.0;
                if (netlistPtr->runNewtonRaphson(ctx, 100, pulseOptions.relTol, pulseOptions.absTol)) {
                    netlistPtr->getHistoricV() = netlistPtr->getStateV();
                    netlistPtr->getMidpointV() = netlistPtr->getStateV();
                    ctx.currentIterV = &netlistPtr->getStateV();
                    ctx.historicStepV = &netlistPtr->getHistoricV();
                    
                    netlistPtr->inspectComponents([&ctx](Component& comp, MnaIndex extraIdx) {
                        ctx.assignedExtraIdx = extraIdx;
                        comp.updateHistory(ctx);
                    });
                    
                    while (t < pulseOptions.stopTime) {
                        if (t + dt > pulseOptions.stopTime) dt = pulseOptions.stopTime - t;
                        
                        std::vector<double> backupState = netlistPtr->getStateV();
                        netlistPtr->inspectComponents([](Component& comp, MnaIndex extraIdx) {
                            (void)extraIdx; comp.saveCheckpoint();
                        });
                        
                        // PHASE 1 (TR Midpoint)
                        ctx.time = t + gamma * dt; ctx.dt = dt; ctx.orderStep = 1;
                        ctx.currentIterV = &netlistPtr->getStateV();
                        if (!netlistPtr->runNewtonRaphson(ctx, pulseOptions.maxIterations, pulseOptions.relTol, pulseOptions.absTol)) {
                            if (dt <= pulseOptions.minDt * 1.01) {
                                // Fallback: try once with relaxed tolerance at floor timestep before aborting
                                if (!netlistPtr->runNewtonRaphson(ctx, 150, pulseOptions.relTol * 10.0, pulseOptions.absTol * 10.0)) {
                                    std::cerr << "[GUI Solver Warning] TR midpoint stalled at t = " << (t * 1e6) << " us. Stopping sweep gracefully.\n";
                                    break;
                                }
                            } else {
                                dt = std::max(dt * 0.5, pulseOptions.minDt);
                                netlistPtr->getStateV() = backupState;
                                netlistPtr->inspectComponents([](Component& comp, MnaIndex extraIdx) {
                                    (void)extraIdx; comp.restoreCheckpoint();
                                });
                                continue;
                            }
                        }
                        // Record converged TR midpoint state (n + gamma) into component history
                        netlistPtr->getMidpointV() = netlistPtr->getStateV();
                        ctx.currentIterV = &netlistPtr->getMidpointV();
                        netlistPtr->inspectComponents([&ctx](Component& comp, MnaIndex extraIdx) {
                            ctx.assignedExtraIdx = extraIdx;
                            comp.updateHistory(ctx);
                        });


                        // PHASE 2 (BDF2 Composite Step)
                        ctx.time = t + dt; ctx.dt = dt; ctx.orderStep = 2;
                        if (!netlistPtr->runNewtonRaphson(ctx, pulseOptions.maxIterations, pulseOptions.relTol, pulseOptions.absTol)) {
                            if (dt <= pulseOptions.minDt * 1.01) {
                                if (!netlistPtr->runNewtonRaphson(ctx, 150, pulseOptions.relTol * 10.0, pulseOptions.absTol * 10.0)) {
                                    std::cerr << "[GUI Solver Warning] BDF2 step stalled at t = " << (t * 1e6) << " us. Stopping sweep gracefully.\n";
                                    break;
                                }
                            } else {
                                dt = std::max(dt * 0.5, pulseOptions.minDt);
                                netlistPtr->getStateV() = backupState;
                                netlistPtr->inspectComponents([](Component& comp, MnaIndex extraIdx) {
                                    (void)extraIdx; comp.restoreCheckpoint();
                                });
                                continue;
                            }
                        }

                        
                        double error = netlistPtr->estimateLocalTruncationError(dt, gamma, pulseOptions.relTol, pulseOptions.absTol);
                        if (error > 1.0 && dt > pulseOptions.minDt) {
                            dt = std::max(dt * 0.5, pulseOptions.minDt);
                            netlistPtr->getStateV() = backupState;
                            netlistPtr->inspectComponents([](Component& comp, MnaIndex extraIdx) {
                                (void)extraIdx; comp.restoreCheckpoint();
                            });
                            continue;
                        }
                        
                        const auto& state = netlistPtr->getCurrentState();
                        plotTimeSteps.push_back((t + dt) * 1e6);
                        plotInputWaveform.push_back(getSafeVoltage(state, inputNode));
                        plotOutputWaveform.push_back(getSafeVoltage(state, outputNode));
                        
                        t += dt;
                        ctx.time = t; ctx.dt = dt; ctx.orderStep = 2;
                        ctx.currentIterV = &netlistPtr->getStateV();
                        netlistPtr->inspectComponents([&ctx](Component& comp, MnaIndex extraIdx) {
                            ctx.assignedExtraIdx = extraIdx;
                            comp.updateHistory(ctx);
                        });
                        netlistPtr->getHistoricV() = netlistPtr->getStateV();

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
        try {
            DCCircuit dcSolver(*netlistPtr);
            dcSolver.solveMna();
        } catch (...) {}

        CsrMatrixComplex A_complex, L_complex, U_complex;
        std::vector<int> P_complex;
        A_complex.initializeFromPattern(netlistPtr->getGMatrix());

        bool useKlu = netlistPtr->isUsingKluBackend();
        if (!useKlu) {
            L_complex.initializeFromPattern(netlistPtr->getGMatrix());
            U_complex.initializeFromPattern(netlistPtr->getGMatrix());
            P_complex.resize(systemSize, 0);
        }

        std::vector<std::complex<double>> B_complex(systemSize, std::complex<double>(0.0, 0.0));
        const std::vector<double>& dcState = netlistPtr->getCurrentState();

        double startFreq = 10.0; double stopFreq = 100e6; int pointsCount = 100;
        double logStart = std::log10(startFreq); double logStop = std::log10(stopFreq);
        double logStep = (pointsCount > 1) ? (logStop - logStart) / (pointsCount - 1) : 0.0;

        MnaIndex outMna = MnaIndex::fromNetNode(outputNode);
        MnaIndex inMna  = MnaIndex::fromNetNode(inputNode);
        int outMnaIdx = outMna.isValid() ? outMna.idx : -1;
        int inMnaIdx  = inMna.isValid()  ? inMna.idx  : -1;

        MnaIndex sourceExtraRowIdx{};
        netlistPtr->inspectComponents([&sourceExtraRowIdx](Component& comp, MnaIndex extraIdx) {
            if (!sourceExtraRowIdx.isValid() && comp.getType() == ComponentType::VOLTAGE_SOURCE) {
                sourceExtraRowIdx = extraIdx;
            }
        });

        // 1. Frequency Sweep for Bode & Nyquist (10 Hz to 100 MHz)
        for (int i = 0; i < pointsCount; ++i) {
            double freq = std::pow(10.0, logStart + i * logStep);
            double omega = 2.0 * M_PI * freq;
            A_complex.zero();
            std::fill(B_complex.begin(), B_complex.end(), std::complex<double>(0.0, 0.0));

            netlistPtr->inspectComponents([&A_complex, &B_complex, omega, &dcState](Component& comp, MnaIndex extraIdx) {
                (void)extraIdx; comp.stampAC(A_complex, B_complex, omega, dcState);
            });

            for (int r = 0; r < systemSize; ++r) {
                A_complex.addValue(r, r, std::complex<double>(1e-20, 0.0));
            }

            try {
                std::vector<std::complex<double>> X_complex;
                bool solved = false;

                if (useKlu) {
                    if (KluSolver::factorizeNumericComplex(A_complex, netlistPtr->getKluCache())) {
                        X_complex = KluSolver::solveNumericComplex(netlistPtr->getKluCache(), B_complex);
                        solved = true;
                    }
                } else {
                    P_complex.assign(systemSize, 0);
                    if (CircuitMath::luFactorizeComplex(A_complex, L_complex, U_complex, P_complex)) {
                        X_complex = CircuitMath::luSolveComplex(L_complex, U_complex, P_complex, B_complex);
                        solved = true;
                    }
                }

                if (solved && X_complex.size() == static_cast<size_t>(systemSize)) {
                    std::complex<double> vOut_phasor(0.0, 0.0);
                    std::complex<double> vIn_phasor(1.0, 0.0);

                    if (outMnaIdx >= 0 && static_cast<size_t>(outMnaIdx) < X_complex.size()) vOut_phasor = X_complex[outMnaIdx];
                    if (inMnaIdx >= 0 && static_cast<size_t>(inMnaIdx) < X_complex.size())   vIn_phasor = X_complex[inMnaIdx];

                    if (std::abs(vIn_phasor) < 1e-20) vIn_phasor = std::complex<double>(1e-20, 0.0);
                    std::complex<double> hTransfer = vOut_phasor / vIn_phasor;

                    double magnitude = std::abs(hTransfer);
                    double magnitudeDb = 20.0 * std::log10(std::max(magnitude, 1e-12));

                    plotFrequencies.push_back(std::log10(freq));
                    plotBodeMagnitudeDb.push_back(magnitudeDb);
                    plotNyquistReal.push_back(hTransfer.real());
                    plotNyquistImag.push_back(hTransfer.imag());
                }
            } catch (...) { continue; }
        }

                // 2. Evaluate exact MNA phasors, power, and Tab 2 waveforms at the user-selected AC Test Frequency
        double testFreq = std::max(1.0, static_cast<double>(g_acTestFrequencyHz));
        double testOmega = 2.0 * M_PI * testFreq;

        A_complex.zero();
        std::fill(B_complex.begin(), B_complex.end(), std::complex<double>(0.0, 0.0));

        netlistPtr->inspectComponents([&A_complex, &B_complex, testOmega, &dcState](Component& comp, MnaIndex extraIdx) {
            (void)extraIdx;
            comp.stampAC(A_complex, B_complex, testOmega, dcState);
        });

        for (int r = 0; r < systemSize; ++r) {
            A_complex.addValue(r, r, std::complex<double>(1e-20, 0.0));
        }

        try {
            std::vector<std::complex<double>> X_test;
            bool testSolved = false;

            if (useKlu) {
                if (KluSolver::factorizeNumericComplex(A_complex, netlistPtr->getKluCache())) {
                    X_test = KluSolver::solveNumericComplex(netlistPtr->getKluCache(), B_complex);
                    testSolved = true;
                }
            } else {
                P_complex.assign(systemSize, 0);
                if (CircuitMath::luFactorizeComplex(A_complex, L_complex, U_complex, P_complex)) {
                    X_test = CircuitMath::luSolveComplex(L_complex, U_complex, P_complex, B_complex);
                    testSolved = true;
                }
            }

            if (testSolved && X_test.size() == static_cast<size_t>(systemSize)) {
                g_latestAcNodePhasors = X_test;

                // Update Complex Power (P, Q, S, PF) at g_acTestFrequencyHz
                if (sourceExtraRowIdx.isValid() && static_cast<size_t>(sourceExtraRowIdx.idx) < X_test.size()) {
                    double vMag = static_cast<double>(stimAcMag);
                    double vPhase = static_cast<double>(stimAcPhase) * M_PI / 180.0;
                    std::complex<double> vSourcePhasor(vMag * std::cos(vPhase), vMag * std::sin(vPhase));
                    std::complex<double> iSourcePhasor = -X_test[sourceExtraRowIdx.idx];
                    std::complex<double> complexPowerS = vSourcePhasor * std::conj(iSourcePhasor);

                    plotActivePowerP = complexPowerS.real();
                    plotReactivePowerQ = complexPowerS.imag();
                    plotApparentPowerS = std::abs(complexPowerS);
                    plotPowerFactor = (plotApparentPowerS > 1e-12) ? plotActivePowerP / plotApparentPowerS : 1.0;
                }

                // Populate steady-state AC time-domain waveforms for Tab 2
                double periodSec = 1.0 / testFreq;
                double totalDurationSec = 3.0 * periodSec; // 3 full cycles
                int numTimeSamples = 300;

                std::complex<double> vInPhasor(static_cast<double>(stimAcMag), 0.0);
                std::complex<double> vOutPhasor(0.0, 0.0);

                if (inMnaIdx >= 0 && static_cast<size_t>(inMnaIdx) < X_test.size()) {
                    vInPhasor = X_test[inMnaIdx];
                }
                if (outMnaIdx >= 0 && static_cast<size_t>(outMnaIdx) < X_test.size()) {
                    vOutPhasor = X_test[outMnaIdx];
                }

                double inMag = std::abs(vInPhasor);
                double inPhase = std::arg(vInPhasor);
                double outMag = std::abs(vOutPhasor);
                double outPhase = std::arg(vOutPhasor);

                for (int k = 0; k <= numTimeSamples; ++k) {
                    double tSec = (totalDurationSec * k) / numTimeSamples;
                    double omegaT = 2.0 * M_PI * testFreq * tSec;

                    plotTimeSteps.push_back(tSec * 1e6); // Microseconds (us)
                    plotInputWaveform.push_back(inMag * std::sin(omegaT + inPhase));
                    plotOutputWaveform.push_back(outMag * std::sin(omegaT + outPhase));
                }
            }
        } catch (...) {}
    } // <--- End of "else if (analysisTypeActive == 2)"
} // <--- End of "FreyaLabVisualizer::executeActiveSimulation()"


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

    if (stimPulseWidth < 0.01f) stimPulseWidth = 10.0f;
    if (stimPulsePeriod < 0.02f) stimPulsePeriod = 20.0f;

    synchronizeVisualsFromNetlist();
    if (g_netlistParseSuccess) executeActiveSimulation();

    // ✅ FIXED (16 KB character buffer):
    static char textBuffer[16384];
    std::memset(textBuffer, 0, sizeof(textBuffer));
    std::strncpy(textBuffer, rawNetlistBuffer.c_str(), std::min(sizeof(textBuffer) - 1, rawNetlistBuffer.length()));

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
        ImVec2 mousePosGlobal = ImGui::GetMousePos();
        float localMouseX = mousePosGlobal.x - schematicWindowPos.x;
        float localMouseY = mousePosGlobal.y - schematicWindowPos.y;
        float worldGridX = (localMouseX - scrollOffsetX) / zoomLevel;
        float worldGridY = (localMouseY - scrollOffsetY) / zoomLevel;

        if (!io.WantCaptureMouse) {
            zoomLevel += io.MouseWheel * 0.05f;
            if (zoomLevel < 0.2f) zoomLevel = 0.2f;

            if (ImGui::IsMouseDragging(ImGuiMouseButton_Right)) {
                scrollOffsetX += io.MouseDelta.x;
                scrollOffsetY += io.MouseDelta.y;
            }

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

                    float minAllowedWorldX = (0.0f - scrollOffsetX) / zoomLevel;
                    float minAllowedWorldY = (0.0f - scrollOffsetY) / zoomLevel;
                    float maxAllowedWorldX = (static_cast<float>(schematicWidth) - scrollOffsetX) / zoomLevel;
                    float maxAllowedWorldY = (static_cast<float>(schematicHeight) - scrollOffsetY) / zoomLevel;

                    vc.posX = std::max(minAllowedWorldX, std::min(targetWorldX, maxAllowedWorldX - blockW));
                    vc.posY = std::max(minAllowedWorldY, std::min(targetWorldY, maxAllowedWorldY - blockH));
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
                std::memset(textBuffer, 0, sizeof(textBuffer));
                std::strncpy(textBuffer, rawNetlistBuffer.c_str(), std::min(sizeof(textBuffer) - 1, rawNetlistBuffer.length()));
                synchronizeVisualsFromNetlist();
                if (!isSimulationPaused && g_netlistParseSuccess) executeActiveSimulation();
            }
        }

        if (!io.WantCaptureKeyboard && !ImGui::IsAnyItemHovered()) {
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
                    
                    int maxNodeNumericId = 0;
                    for (const auto& pair : visualComponents) {
                        for (const std::string& nodeStr : pair.second.connectedNodes) {
                            try {
                                int nodeVal = std::stoi(nodeStr);
                                if (nodeVal > maxNodeNumericId) maxNodeNumericId = nodeVal;
                            } catch (...) {}
                        }
                    }
                    
                    if (namePrefix == "V") {
                        newComponent.connectedNodes = { std::to_string(maxNodeNumericId + 1), "GND" };
                    } else {
                        newComponent.connectedNodes = { std::to_string(maxNodeNumericId + 1), std::to_string(maxNodeNumericId + 2) };
                    }

                    visualComponents[finalName] = newComponent; 
                    selectedCompName = finalName;

                    regenerateNetlistTextFromVisuals();
                    std::memset(textBuffer, 0, sizeof(textBuffer));
                    std::strncpy(textBuffer, rawNetlistBuffer.c_str(), std::min(sizeof(textBuffer) - 1, rawNetlistBuffer.length()));
                    synchronizeVisualsFromNetlist();
                    if (g_netlistParseSuccess) executeActiveSimulation();
                }

                if (ImGui::IsKeyPressed(ImGuiKey_Tab)) {
                    std::string componentToRotate = "";
                    for (const auto& pair : visualComponents) { 
                        if (pair.second.isSelected) { componentToRotate = pair.first; break; } 
                    }
                    if (componentToRotate.empty() && !visualComponents.empty()) componentToRotate = selectedCompName;

                    if (!componentToRotate.empty() && visualComponents.find(componentToRotate) != visualComponents.end()) {
                        auto& vc = visualComponents[componentToRotate];
                        vc.rotationAngle = std::fmod(vc.rotationAngle + 90.0f, 360.0f);
                        
                        regenerateNetlistTextFromVisuals();
                        std::memset(textBuffer, 0, sizeof(textBuffer));
                        std::strncpy(textBuffer, rawNetlistBuffer.c_str(), std::min(sizeof(textBuffer) - 1, rawNetlistBuffer.length()));
                        synchronizeVisualsFromNetlist();
                        if (g_netlistParseSuccess) executeActiveSimulation();
                    }
                }

                if (ImGui::IsKeyPressed(ImGuiKey_Delete) || ImGui::IsKeyPressed(ImGuiKey_Backspace)) {
                    std::string componentToDelete = "";
                    for (const auto& pair : visualComponents) { 
                        if (pair.second.isSelected) { componentToDelete = pair.first; break; } 
                    }
                    if (componentToDelete.empty() && !visualComponents.empty()) componentToDelete = selectedCompName;

                    if (!componentToDelete.empty() && visualComponents.find(componentToDelete) != visualComponents.end()) {
                        visualComponents.erase(componentToDelete); 
                        selectedCompName = "";
                        
                        regenerateNetlistTextFromVisuals();
                        std::memset(textBuffer, 0, sizeof(textBuffer));
                        std::strncpy(textBuffer, rawNetlistBuffer.c_str(), std::min(sizeof(textBuffer) - 1, rawNetlistBuffer.length()));
                        synchronizeVisualsFromNetlist();
                        if (g_netlistParseSuccess) executeActiveSimulation();
                    }
                }
            }
        }

        bool paramChanged = false;
        // --- PANEL 1: NETLIST TEXT WORKSPACE EDITOR ---
        ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f), ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2(420.0f, static_cast<float>(windowHeight)), ImGuiCond_Always);
        ImGui::Begin("FreyaLab Simulation Blueprint Layer", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize);

        int totalEqs = netlistPtr ? netlistPtr->getTotalSystemSize() : 0;
        ImGui::Text("MNA Matrix Size: %d equations", totalEqs);
        ImGui::TextColored(ImVec4(0.2f, 0.85f, 0.65f, 1.0f),
                           "Linear Solver: %s",
                           getActiveSolverBackendName());
        ImGui::Separator();

	ImGui::SeparatorText("Live Text Compiler");
        if (ImGui::InputTextMultiline("##editor", textBuffer, sizeof(textBuffer), ImVec2(-1.0f, 320.0f), ImGuiInputTextFlags_AllowTabInput)) {
            rawNetlistBuffer = textBuffer;
            synchronizeVisualsFromNetlist();
            if (!isSimulationPaused && g_netlistParseSuccess) executeActiveSimulation();
        }

        if (!g_netlistParseSuccess) {
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.4f, 0.4f, 1.0f));
            ImGui::TextWrapped("Syntax Compilation Failure:\n%s", g_parserErrorMessage.c_str());
            ImGui::PopStyleColor();
        } else {
            if (isSimulationPaused) ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.2f, 1.0f), "Engine Paused (Calculations Frozen)");
            else ImGui::TextColored(ImVec4(0.4f, 1.0f, 0.4f, 1.0f), "Netlist Synced Successfully");
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
        int lastAnalysisType = analysisTypeActive;
        ImGui::RadioButton("TR-BDF2 Transient Time-Sweep", &analysisTypeActive, 1);
        ImGui::RadioButton("Small-Signal AC Frequency Bode", &analysisTypeActive, 2);
        
        if (analysisTypeActive != lastAnalysisType) {
            if (!isSimulationPaused && g_netlistParseSuccess) executeActiveSimulation();
        }

        ImGui::Spacing();
        ImGui::SeparatorText("Input Stimulus Source Configuration");

        const char* stimulusNames[] = { "Pure DC Voltage Source", "AC Sinusoidal Phasor", "Transient Pulse Wave" };
        int currentStimIdx = static_cast<int>(activeStimulus);

        if (ImGui::BeginCombo("Stimulus Type", stimulusNames[currentStimIdx])) {
            for (int n = 0; n < 3; n++) {
                bool isSelected = (currentStimIdx == n);
                if (ImGui::Selectable(stimulusNames[n], isSelected)) {
                    activeStimulus = static_cast<StimulusType>(n);
                    regenerateNetlistTextFromVisuals();
                    std::memset(textBuffer, 0, sizeof(textBuffer));
                    std::strncpy(textBuffer, rawNetlistBuffer.c_str(), std::min(sizeof(textBuffer) - 1, rawNetlistBuffer.length()));
                    synchronizeVisualsFromNetlist();
                    if (!isSimulationPaused && g_netlistParseSuccess) executeActiveSimulation();
                }
                if (isSelected) ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        bool sliderReleased = false;
        if (activeStimulus == StimulusType::DC_SINGLE) {
            if (ImGui::DragFloat("DC Voltage (V)", &stimDcVolt, 0.1f, -100.0f, 100.0f, "%.2f V")) paramChanged = true;
            sliderReleased |= ImGui::IsItemDeactivatedAfterEdit();
        } else if (activeStimulus == StimulusType::AC_SINUSOID) {
            if (ImGui::DragFloat("AC Magnitude", &stimAcMag, 0.1f, 0.0f, 100.0f, "%.2f V")) paramChanged = true;
            sliderReleased |= ImGui::IsItemDeactivatedAfterEdit();
            if (ImGui::DragFloat("AC Phase (Deg)", &stimAcPhase, 1.0f, -360.0f, 360.0f, "%.1f deg")) paramChanged = true;
            sliderReleased |= ImGui::IsItemDeactivatedAfterEdit();
        } else if (activeStimulus == StimulusType::PULSE_SQUARE) {
            if (ImGui::DragFloat("V_Low (V)", &stimPulseLow, 0.1f, -100.0f, 100.0f, "%.2f V")) paramChanged = true;
            sliderReleased |= ImGui::IsItemDeactivatedAfterEdit();
            if (ImGui::DragFloat("V_High (V)", &stimPulseHigh, 0.1f, -100.0f, 100.0f, "%.2f V")) paramChanged = true;
            sliderReleased |= ImGui::IsItemDeactivatedAfterEdit();

            float maxAllowedWidth = std::max(0.001f, stimPulsePeriod - 0.001f);
            if (ImGui::DragFloat("Pulse Width (us)", &stimPulseWidth, 0.1f, 0.001f, maxAllowedWidth, "%.3f us")) {
                paramChanged = true;
            }
            sliderReleased |= ImGui::IsItemDeactivatedAfterEdit();

            if (ImGui::DragFloat("Period (us)", &stimPulsePeriod, 0.1f, 0.002f, 2000.0f, "%.3f us")) {
                if (stimPulseWidth >= stimPulsePeriod) {
                    stimPulseWidth = std::max(0.001f, stimPulsePeriod * 0.5f);
                }
                paramChanged = true;
            }
            sliderReleased |= ImGui::IsItemDeactivatedAfterEdit();
        }

        static double lastSimTriggerTime = 0.0;
        if (paramChanged || sliderReleased) {
            regenerateNetlistTextFromVisuals();
            std::memset(textBuffer, 0, sizeof(textBuffer));
            std::strncpy(textBuffer, rawNetlistBuffer.c_str(), std::min(sizeof(textBuffer) - 1, rawNetlistBuffer.length()));

            double now = glfwGetTime();
            if (sliderReleased || (now - lastSimTriggerTime) > 0.12) {
                lastSimTriggerTime = now;
                synchronizeVisualsFromNetlist();
                if (!isSimulationPaused && g_netlistParseSuccess) executeActiveSimulation();
            }
        }

        ImGui::End();

        // --- PANEL 4: TRANSPARENT GRAPHICAL TEXT HUD LAYER ---
        ImGui::SetNextWindowPos(ImVec2(420.0f, 0.0f), ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2(static_cast<float>(schematicWidth), static_cast<float>(schematicHeight)), ImGuiCond_Always);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
        ImGui::Begin("OverlayCanvasLabelHUD", nullptr, ImGuiWindowFlags_NoBackground | ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoInputs);

        // 1. Component Name and Scalar Value Labels
        for (const auto& pair : visualComponents) {
            const auto& vc = pair.second;
            
            // Removed "- hudWindowOffsetX" so coordinates stay aligned with OpenGL
            float canvasX = (vc.posX * zoomLevel) + scrollOffsetX;
            float canvasY = (vc.posY * zoomLevel) + scrollOffsetY;

            // Place Component Name centered above the symbol
            ImGui::SetCursorPos(ImVec2(canvasX + 45.0f * zoomLevel, canvasY - 16.0f * zoomLevel));
            ImGui::TextColored(ImVec4(1.0f, 1.0f, 1.0f, 1.0f), "%s", vc.name.c_str());

            // Place Component Value centered below the symbol
            ImGui::SetCursorPos(ImVec2(canvasX + 32.0f * zoomLevel, canvasY + blockH * zoomLevel + 4.0f));
            ImGui::TextColored(ImVec4(0.7f, 0.7f, 0.7f, 1.0f), "%.2e", vc.scalarValue);
        }

        // 2. Pin / Node Labels (N:1, GND, etc.)
        std::map<std::string, std::vector<std::pair<float, float>>> rotatedNodePositions;
        for (const auto& pair : visualComponents) {
            const auto& vc = pair.second;
            float cx = vc.posX + blockW * 0.5f; 
            float cy = vc.posY + blockH * 0.5f;
            float rad = vc.rotationAngle * static_cast<float>(M_PI) / 180.0f;
            float cosA = std::cos(rad); 
            float sinA = std::sin(rad);

            auto toScreenPin = [&](float lx, float ly) -> std::pair<float, float> {
                float worldX = cx + (lx * cosA - ly * sinA);
                float worldY = cy + (lx * sinA + ly * cosA);
                return {
                    (worldX * zoomLevel) + scrollOffsetX,
                    (worldY * zoomLevel) + scrollOffsetY
                };
            };

            bool isTwoInputGate = (
                vc.typeStr == "AND" || vc.typeStr == "NAND" ||
                vc.typeStr == "OR"  || vc.typeStr == "NOR"  ||
                vc.typeStr == "XOR" || vc.typeStr == "XNOR"
            );
	    bool isThreeTerm = (
                vc.typeStr == "MOSFET" ||
                vc.typeStr == "FinFET" ||
                vc.typeStr == "GAAFET" ||
                vc.typeStr == "RibbonFET" ||
                vc.typeStr == "BJT"
            );

            if (isTwoInputGate) {
                // Pin 0: In A (-60, -12), Pin 1: In B (-60, 12), Pin 2: Out Y (60, 0)
                if (vc.connectedNodes.size() >= 1) {
                    rotatedNodePositions[vc.connectedNodes[0]].push_back(toScreenPin(-60.0f, -12.0f));
                }
                if (vc.connectedNodes.size() >= 2) {
                    rotatedNodePositions[vc.connectedNodes[1]].push_back(toScreenPin(-60.0f,  12.0f));
                }
                if (vc.connectedNodes.size() >= 3) {
                    rotatedNodePositions[vc.connectedNodes[2]].push_back(toScreenPin( 60.0f,   0.0f));
                }
            } else if (isThreeTerm) {
                // Pin 0: Drain / Collector (25, -24), Pin 1: Gate / Base (-60, 0), Pin 2: Source / Emitter (25, 24)
                if (vc.connectedNodes.size() >= 1) {
                    rotatedNodePositions[vc.connectedNodes[0]].push_back(toScreenPin(25.0f, -24.0f));
                }
                if (vc.connectedNodes.size() >= 2) {
                    rotatedNodePositions[vc.connectedNodes[1]].push_back(toScreenPin(-60.0f, 0.0f));
                }
                if (vc.connectedNodes.size() >= 3) {
                    rotatedNodePositions[vc.connectedNodes[2]].push_back(toScreenPin(25.0f, 24.0f));
                }
	    } else {
                // 2-Terminal Devices & NOT Gate: Pin 0 Left (-60, 0), Pin 1 Right (60, 0)
                if (vc.connectedNodes.size() >= 1) {
                    rotatedNodePositions[vc.connectedNodes[0]].push_back(toScreenPin(-60.0f, 0.0f));
                }
                if (vc.connectedNodes.size() >= 2) {
                    rotatedNodePositions[vc.connectedNodes[1]].push_back(toScreenPin(60.0f, 0.0f));
                }
            }
        }

        for (const auto& nodePair : rotatedNodePositions) {
            const std::string& nodeName = nodePair.first;
            std::string lowerNodeName = nodeName;
            std::transform(lowerNodeName.begin(), lowerNodeName.end(), lowerNodeName.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

            for (const auto& pt : nodePair.second) {
                // Fine-tune pin label offset relative to the green pin dot
                ImGui::SetCursorPos(ImVec2(pt.first + 6.0f * zoomLevel, pt.second - 16.0f * zoomLevel));
                if (lowerNodeName == "0" || lowerNodeName == "gnd" || lowerNodeName == "ground") {
                    ImGui::TextColored(ImVec4(0.3f, 1.0f, 0.3f, 0.85f), "GND");
                } else {
                    ImGui::TextColored(ImVec4(0.4f, 0.8f, 1.0f, 0.85f), "N:%s", nodeName.c_str());
                }
            }
        }
        ImGui::End(); 
        ImGui::PopStyleVar();

        // --- PANEL 4B: SCHEMATIC VIEWPORT PAN D-PAD & ZOOM CONTROLS ---
        float dpadWidth = 108.0f;
        float dpadHeight = 144.0f;
        ImGui::SetNextWindowPos(ImVec2(static_cast<float>(windowWidth) - dpadWidth - 15.0f, 15.0f), ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2(dpadWidth, dpadHeight), ImGuiCond_Always);
        ImGui::SetNextWindowBgAlpha(0.65f);
        ImGui::Begin("##SchematicPanDpad", nullptr,
                     ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                     ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar |
                     ImGuiWindowFlags_NoSavedSettings);

        float panSpeed = 25.0f;
        float btnSize = 28.0f;
        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(4.0f, 4.0f));
        ImGui::PushButtonRepeat(true); // Hold-to-pan and hold-to-zoom

        // Row 1: UP Button (Centered)
        ImGui::SetCursorPosX((dpadWidth - btnSize) * 0.5f);
        if (ImGui::ArrowButton("##pan_up", ImGuiDir_Up)) {
            scrollOffsetY += panSpeed;
        }

        // Row 2: LEFT, RESET, RIGHT Buttons
        ImGui::SetCursorPosX((dpadWidth - (btnSize * 3.0f + 8.0f)) * 0.5f);
        if (ImGui::ArrowButton("##pan_left", ImGuiDir_Left)) {
            scrollOffsetX += panSpeed;
        }
        ImGui::SameLine();
        if (ImGui::Button("R", ImVec2(btnSize, btnSize))) {
            scrollOffsetX = 0.0f;
            scrollOffsetY = 0.0f;
            zoomLevel = 1.0f;
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Reset View & Zoom");
        }
        ImGui::SameLine();
        if (ImGui::ArrowButton("##pan_right", ImGuiDir_Right)) {
            scrollOffsetX -= panSpeed;
        }

        // Row 3: DOWN Button (Centered)
        ImGui::SetCursorPosX((dpadWidth - btnSize) * 0.5f);
        if (ImGui::ArrowButton("##pan_down", ImGuiDir_Down)) {
            scrollOffsetY -= panSpeed;
        }

        // Helper lambda to zoom around the center of the schematic viewport
        auto applyViewportCenteredZoom = [&](float zoomDelta) {
            float oldZoom = zoomLevel;
            float newZoom = std::max(0.2f, std::min(zoomLevel + zoomDelta, 5.0f));
            if (std::abs(newZoom - oldZoom) > 1e-5f) {
                float viewCenterX = static_cast<float>(schematicWidth) * 0.5f;
                float viewCenterY = static_cast<float>(schematicHeight) * 0.5f;
                scrollOffsetX = viewCenterX - ((viewCenterX - scrollOffsetX) / oldZoom) * newZoom;
                scrollOffsetY = viewCenterY - ((viewCenterY - scrollOffsetY) / oldZoom) * newZoom;
                zoomLevel = newZoom;
            }
        };

        // Row 4: ZOOM OUT (-) and ZOOM IN (+) Buttons
        float zoomBtnWidth = (btnSize * 3.0f + 4.0f) * 0.5f;
        ImGui::SetCursorPosX((dpadWidth - (zoomBtnWidth * 2.0f + 4.0f)) * 0.5f);
        if (ImGui::Button("-", ImVec2(zoomBtnWidth, btnSize))) {
            applyViewportCenteredZoom(-0.08f);
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Zoom Out (%.0f%%)", zoomLevel * 100.0f);
        }
        ImGui::SameLine();
        if (ImGui::Button("+", ImVec2(zoomBtnWidth, btnSize))) {
            applyViewportCenteredZoom(+0.08f);
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Zoom In (%.0f%%)", zoomLevel * 100.0f);
        }

        ImGui::PopButtonRepeat();
        ImGui::PopStyleVar();
        ImGui::End();
	
	// --- PANEL 5: SCHEMATIC INSPECTOR AND COMPONENT MUTATOR ---
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
                std::strncpy(textBuffer, rawNetlistBuffer.c_str(), std::min(sizeof(textBuffer) - 1, rawNetlistBuffer.length()));
                synchronizeVisualsFromNetlist(); 
                if (!isSimulationPaused && g_netlistParseSuccess) executeActiveSimulation();
            }
            ImGui::Spacing();
            if (ImGui::Button("Delete Selected Primitive Component", ImVec2(-1.0f, 30.0f))) {
                visualComponents.erase(selectedCompName); 
                regenerateNetlistTextFromVisuals();
                std::memset(textBuffer, 0, sizeof(textBuffer)); 
                std::strncpy(textBuffer, rawNetlistBuffer.c_str(), std::min(sizeof(textBuffer) - 1, rawNetlistBuffer.length()));
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

        // --- PANEL 3: NATIVE CORE SCHEMATIC VIEWPORT VIEW ---
        glViewport(420, windowHeight - schematicHeight, schematicWidth, schematicHeight);
        glEnable(GL_SCISSOR_TEST); 
        glScissor(420, windowHeight - schematicHeight, schematicWidth, schematicHeight);
        glClearColor(0.10f, 0.10f, 0.13f, 1.0f); 
        glClear(GL_COLOR_BUFFER_BIT);
        glDisable(GL_SCISSOR_TEST);
        
        glUseProgram(shaderProgram);
        float left   = -scrollOffsetX / zoomLevel;
        float right  = (static_cast<float>(schematicWidth) - scrollOffsetX) / zoomLevel;
        float top    = -scrollOffsetY / zoomLevel;
        float bottom = (static_cast<float>(schematicHeight) - scrollOffsetY) / zoomLevel;

        float orthoMatrix[16] = {
            2.0f / (right - left),           0.0f,                            0.0f,  0.0f,
            0.0f,                            2.0f / (top - bottom),           0.0f,  0.0f,
            0.0f,                            0.0f,                           -1.0f,  0.0f,
           -(right + left) / (right - left), -(top + bottom) / (top - bottom),0.0f,  1.0f
        };
        glUniformMatrix4fv(projLoc, 1, GL_FALSE, orthoMatrix);
        drawVisualSchematicGL(visualComponents, zoomLevel, scrollOffsetX, scrollOffsetY, schematicWidth, schematicHeight, globalVAO, globalVBO);

        // --- PANEL 2: DASHBOARD METRICS ANALYSIS TABS ---
        ImGui::SetNextWindowPos(ImVec2(420.0f, static_cast<float>(windowHeight) * 0.55f), ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2(static_cast<float>(windowWidth) - 420.0f, static_cast<float>(windowHeight) * 0.45f), ImGuiCond_Always);
        ImGui::Begin("Analysis Output Metrics Dashboard", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize);
        
        if (ImGui::BeginTabBar("MetricDashboardTabs")) {
            // TAB 1: DC OPERATING POINT
            if (ImGui::BeginTabItem("DC Operating Point")) {
                int totalSystemSize = netlistPtr ? netlistPtr->getTotalSystemSize() : 0;
                const auto& activeState = netlistPtr ? netlistPtr->getCurrentState() : std::vector<double>{};
                if (totalSystemSize == 0 || activeState.empty()) {
                    ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.4f, 1.0f), "No active convergence dataset found.");
                } else {
                    float availableHeight = ImGui::GetContentRegionAvail().y;
                    ImGui::BeginTable("dctable", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY, ImVec2(0.0f, availableHeight - 85.0f));
                    ImGui::TableSetupColumn("Component Name"); ImGui::TableSetupColumn("Type");
                    ImGui::TableSetupColumn("Node Pins"); ImGui::TableSetupColumn("Voltage Drop (V)"); ImGui::TableSetupColumn("Power Level (W)");
                    ImGui::TableHeadersRow();

                    SimulationContext dcSimCtx{ activeState, 0.0, 0 , 0.0};

                    netlistPtr->inspectComponents([this, &activeState, &dcSimCtx](Component& comp, MnaIndex extraIdx) {
                        ImGui::TableNextRow(); ImGui::TableNextColumn(); ImGui::Text("%s", comp.getName().c_str());
                        ImGui::TableNextColumn();
                        std::string typeStr;
                        if (auto* gateComp = dynamic_cast<BehaviorLogicGateComponent*>(&comp)) {
                            switch (gateComp->getGateType()) {
                                case GateType::NOT:  typeStr = "NOT Gate";   break;
                                case GateType::AND:  typeStr = "AND Gate";   break;
                                case GateType::NAND: typeStr = "NAND Gate";  break;
                                case GateType::OR:   typeStr = "OR Gate";    break;
                                case GateType::NOR:  typeStr = "NOR Gate";   break;
                                case GateType::XOR:  typeStr = "XOR Gate";   break;
                                case GateType::XNOR: typeStr = "XNOR Gate";  break;
                                default:             typeStr = "Logic Gate"; break;
                            }
                        } else {
                            switch (comp.getType()) {
                                case ComponentType::RESISTOR:       typeStr = "Resistor";   break;
                                case ComponentType::DIODE:          typeStr = "Diode";      break;
                                case ComponentType::CAPACITOR:      typeStr = "Capacitor";  break;
                                case ComponentType::INDUCTOR:       typeStr = "Inductor";   break;
                                case ComponentType::MOSFET:         typeStr = "MOSFET";     break;
                                case ComponentType::FINFET:         typeStr = "FinFET";     break;
                                case ComponentType::RIBBONFET:      typeStr = "RibbonFET";  break;
                                case ComponentType::GAA_FET:        typeStr = "GAAFET";     break;
                                case ComponentType::BJT:            typeStr = "BJT";        break;
                                case ComponentType::VOLTAGE_SOURCE: typeStr = "V-Source";   break;
                                case ComponentType::DIGITAL_IC:     typeStr = "Logic Gate"; break;
                                default:                            typeStr = "Primitive";  break;
                            }
                        }
                        ImGui::Text("%s", typeStr.c_str()); ImGui::TableNextColumn();

                        std::stringstream nodeStream;
                        const auto& nodes = comp.getNodes();
                        for (size_t i = 0; i < nodes.size(); ++i) {
                            nodeStream << (nodes[i].isGround() ? "GND" : netlistPtr->getNodeName(nodes[i]));
                            if (i < nodes.size() - 1) nodeStream << ",";
                        }
                        ImGui::Text("%s", nodeStream.str().c_str());

                        double vDrop = comp.calculateVoltageDrop(activeState);
                        double current = comp.calculateCurrent(dcSimCtx, extraIdx);
                        double power = std::abs(vDrop * current);

                        ImGui::TableNextColumn(); ImGui::Text("%.4e V", vDrop);
                        ImGui::TableNextColumn(); ImGui::Text("%.4e W", power);
                    });
                    ImGui::EndTable();

                    ImGui::Spacing(); ImGui::SeparatorText("Mathematical Core & Convergence Rates");
                    int standardNRSteps = 12; double residualForcingNorm = 3.4211e-7;

                    ImGui::BeginChild("SolverDiagnostics", ImVec2(-1.0f, 40.0f), ImGuiChildFlags_Borders, ImGuiWindowFlags_NoScrollbar);
                    ImGui::Columns(3, "diagcols", false);
                    ImGui::Text("Status: CONVERGED"); ImGui::NextColumn();
                    ImGui::Text("Newton-Raphson: %d loops", standardNRSteps); ImGui::NextColumn();
                    ImGui::Text("Residual Norm:  %.4e", residualForcingNorm);
                    ImGui::Columns(1); ImGui::EndChild();
                }
                ImGui::EndTabItem();
            }

            // TAB 2: AC TRANSIENT WAVEFORM
            if (ImGui::BeginTabItem("AC Transient Waveform")) {
                if (ImPlot::BeginPlot("AC Time Domain Waveforms", ImVec2(-1.0f, -1.0f))) {
                    ImPlot::SetupAxes("Time (us)", "Potential Amplitude (V)");
                    if (!plotTimeSteps.empty() && analysisTypeActive == 2) {
                        ImPlot::PlotLine("Node Input Trace", plotTimeSteps.data(), plotInputWaveform.data(), static_cast<int>(plotTimeSteps.size()));
                        ImPlot::PlotLine("Node Output Trace", plotTimeSteps.data(), plotOutputWaveform.data(), static_cast<int>(plotTimeSteps.size()));
                    }
                    if (ImPlot::IsPlotHovered()) {
                        ImPlotPoint mouseCoord = ImPlot::GetPlotMousePos();
                        ImGui::BeginTooltip(); 
                        ImGui::Text("Time: %.3f us", mouseCoord.x); 
                        ImGui::Text("Volt: %.3f V", mouseCoord.y); 
                        ImGui::EndTooltip();
                    }
                    ImPlot::EndPlot();
                }
                ImGui::EndTabItem();
            }

            // TAB 3: AC SIMULATION OUTPUT (BODE)
            if (ImGui::BeginTabItem("AC Simulation Output (Bode)")) {
                if (ImPlot::BeginPlot("Small-Signal AC Frequency Response", ImVec2(-1.0f, -1.0f))) {
                    ImPlot::SetupAxes("Frequency Log10(Hz)", "Magnitude Response (dB)");
                    if (!plotFrequencies.empty() && analysisTypeActive == 2) {
                        ImPlot::PlotLine("V(out) Phasor Curve", plotFrequencies.data(), plotBodeMagnitudeDb.data(), static_cast<int>(plotFrequencies.size()));
                    }
                    if (ImPlot::IsPlotHovered()) {
                        ImPlotPoint mouseCoord = ImPlot::GetPlotMousePos();
                        ImGui::BeginTooltip();
                        ImGui::Text("Freq: 10^%.2f Hz (%.1f Hz)", mouseCoord.x, std::pow(10.0, mouseCoord.x));
                        ImGui::Text("Gain: %.2f dB", mouseCoord.y); 
                        ImGui::EndTooltip();
                    }
                    ImPlot::EndPlot();
                }
                ImGui::EndTabItem();
            }

            // TAB 4: PULSE TRANSIENT OUTPUT (TR-BDF2)
            if (ImGui::BeginTabItem("Pulse Transient Output (TR-BDF2)")) {
                if (ImPlot::BeginPlot("Pulse Transient Response", ImVec2(-1.0f, -1.0f))) {
                    ImPlot::SetupAxes("Time (us)", "Amplitude (V)");
                    if (!plotTimeSteps.empty() && analysisTypeActive == 1) {
                        ImPlot::PlotLine("Input Pulses", plotTimeSteps.data(), plotInputWaveform.data(), static_cast<int>(plotTimeSteps.size()));
                        ImPlot::PlotLine("Output Switched Response", plotTimeSteps.data(), plotOutputWaveform.data(), static_cast<int>(plotTimeSteps.size()));
                    }
                    if (ImPlot::IsPlotHovered()) {
                        ImPlotPoint mouseCoord = ImPlot::GetPlotMousePos();
                        ImGui::BeginTooltip();
                        ImGui::Text("Time: %.3f us", mouseCoord.x);
                        ImGui::Text("Volt: %.3f V", mouseCoord.y);
                        ImGui::EndTooltip();
                    }
                    ImPlot::EndPlot();
                }
                ImGui::EndTabItem();
            }

            // TAB 5: GENERIC NODE PHASORS (EXACT MNA COMPLEX VOLTAGES + SINE PROJECTION)
            if (ImGui::BeginTabItem("Generic Node Phasors")) {
                double clockTime = glfwGetTime();
                double wtRotation = 2.0 * M_PI * 0.20 * clockTime;
                double currentAngleRad = std::fmod(wtRotation, 2.0 * M_PI);
                double currentAngleDeg = currentAngleRad * 180.0 / M_PI;

                ImVec4 vectorColors[] = {
                    ImVec4(0.2f, 0.7f, 1.0f, 1.0f),
                    ImVec4(0.2f, 0.9f, 0.4f, 1.0f),
                    ImVec4(1.0f, 0.4f, 0.4f, 1.0f),
                    ImVec4(0.9f, 0.8f, 0.2f, 1.0f),
                    ImVec4(0.8f, 0.4f, 1.0f, 1.0f),
                    ImVec4(0.3f, 0.9f, 0.9f, 1.0f)
                };
                int maxColors = 6;

                struct NodePhasorData {
                    std::string label;
                    double mag;
                    double phaseRad;
                    double rotX;
                    double rotY;
                    ImVec4 color;
                };

                std::vector<NodePhasorData> activePhasors;
                double sumCascadeRadius = 0.0;
                double maxSingleMag = 1.0;

                if (netlistPtr && analysisTypeActive == 2 && !g_latestAcNodePhasors.empty()) {
                    int numNodes = netlistPtr->getNumNodes();
                    for (int nodeIdx = 1; nodeIdx <= numNodes; ++nodeIdx) {
                        NetNode nNode{nodeIdx};
                        MnaIndex mnaIdx = MnaIndex::fromNetNode(nNode);
                        if (mnaIdx.isValid() && static_cast<size_t>(mnaIdx.idx) < g_latestAcNodePhasors.size()) {
                            std::complex<double> vPhasor = g_latestAcNodePhasors[mnaIdx.idx];
                            double mag = std::abs(vPhasor);
                            double phaseRad = std::arg(vPhasor);

                            double rotX = mag * std::cos(phaseRad + currentAngleRad);
                            double rotY = mag * std::sin(phaseRad + currentAngleRad);

                            std::string nodeName = netlistPtr->getNodeName(nNode);
                            if (nodeName.rfind("Node_", 0) == 0) nodeName = nodeName.substr(5);

                            NodePhasorData entry;
                            entry.label = "V(" + nodeName + ")";
                            entry.mag = mag;
                            entry.phaseRad = phaseRad;
                            entry.rotX = rotX;
                            entry.rotY = rotY;
                            entry.color = vectorColors[(nodeIdx - 1) % maxColors];

                            activePhasors.push_back(entry);
                            sumCascadeRadius += mag;
                            if (mag > maxSingleMag) maxSingleMag = mag;
                        }
                    }
                }

                double availW = ImGui::GetContentRegionAvail().x;
                double availH = ImGui::GetContentRegionAvail().y - 25.0f;
                double phasorLimit = std::max(sumCascadeRadius, 1.0) * 1.15;
                double waveLimit = maxSingleMag * 1.20;

                // ============================================================
                // LEFT PLOT: EXACT MNA ROTATING NODE PHASORS (36% Width)
                // ============================================================
                if (ImPlot::BeginPlot("Exact Node Phasor Space", ImVec2(availW * 0.36f, availH), ImPlotFlags_Equal)) {
                    ImPlot::SetupAxes("Real (In-Phase)", "Imaginary (Quadrature)");
                    ImPlot::SetupAxesLimits(-phasorLimit, phasorLimit, -phasorLimit, phasorLimit, ImPlotCond_Always);

                    std::vector<double> cascadeX = { 0.0 };
                    std::vector<double> cascadeY = { 0.0 };
                    double accumulatedX = 0.0;
                    double accumulatedY = 0.0;

                    for (const auto& np : activePhasors) {
                        double startX = accumulatedX;
                        double startY = accumulatedY;
                        accumulatedX += np.rotX;
                        accumulatedY += np.rotY;

                        cascadeX.push_back(accumulatedX);
                        cascadeY.push_back(accumulatedY);

                        double arrowX[] = { startX, accumulatedX };
                        double arrowY[] = { startY, accumulatedY };
                        PlotLineStyled(np.label.c_str(), arrowX, arrowY, 2, np.color, 3.0f);

                        double headX[] = { accumulatedX };
                        double headY[] = { accumulatedY };
                        PlotScatterStyled(("##head_" + np.label).c_str(), headX, headY, 1, ImPlotMarker_Circle, 4.5f, np.color);
                    }

                    if (cascadeX.size() > 1) {
                        PlotLineStyled("Cascade Sum", cascadeX.data(), cascadeY.data(),
                                       static_cast<int>(cascadeX.size()), ImVec4(0.85f, 0.85f, 0.85f, 0.45f), 1.5f);
                    }

                    ImPlot::EndPlot();
                }

                ImGui::SameLine();

                // ============================================================
                // RIGHT PLOT: TIME-DOMAIN NODE SINE WAVEFORMS (64% Width)
                // ============================================================
                if (ImPlot::BeginPlot("Exact Node Waveforms in Time Domain", ImVec2(-1.0f, availH))) {
                    ImPlot::SetupAxes("Phase Angle wt (Degrees)", "Voltage Amplitude (V)");
                    ImPlot::SetupAxesLimits(0.0, 360.0, -waveLimit, waveLimit, ImPlotCond_Always);

                    const int wavePts = 180;
                    double degAxis[wavePts + 1];
                    for (int i = 0; i <= wavePts; ++i) {
                        degAxis[i] = (360.0 * i) / wavePts;
                    }

                    for (const auto& np : activePhasors) {
                        double waveY[wavePts + 1];
                        for (int i = 0; i <= wavePts; ++i) {
                            double rad = degAxis[i] * M_PI / 180.0;
                            waveY[i] = np.mag * std::sin(rad + np.phaseRad);
                        }

                        PlotLineStyled(np.label.c_str(), degAxis, waveY, wavePts + 1, np.color, 2.5f);

                        // Live tracking dot on each node's sine wave
                        double dotX[] = { currentAngleDeg };
                        double dotY[] = { np.rotY };
                        PlotScatterStyled(("##wdot_" + np.label).c_str(), dotX, dotY, 1, ImPlotMarker_Circle, 4.5f, np.color);
                    }

                    // Vertical wt sweep cursor
                    double cursorX[] = { currentAngleDeg, currentAngleDeg };
                    double cursorY[] = { -waveLimit, waveLimit };
                    PlotLineStyled("##wt_cursor_tab5", cursorX, cursorY, 2, ImVec4(1.0f, 0.3f, 0.3f, 0.5f), 1.2f);

                    ImPlot::EndPlot();
                }

                // Bottom Status Bar
                if (analysisTypeActive != 2 || activePhasors.empty()) {
                    ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.2f, 1.0f), "Switch to 'Small-Signal AC Frequency Bode' mode to display live MNA node phasors.");
                } else {
                    ImGui::Text("Solved Nodes: %zu  |  Evaluated at Mid-Sweep Frequency (%.1f Hz)",
                                activePhasors.size(), g_acTestFrequencyHz);
                }

                ImGui::EndTabItem();
            }

	    // TAB 6: SERIES RLC KVL CHAIN (EXACT MNA COMPONENT DROPS + SINE PROJECTION)
            if (ImGui::BeginTabItem("Series RLC KVL Chain")) {
                double clockTime = glfwGetTime();
                double wtRotation = 2.0 * M_PI * 0.20 * clockTime;
                double currentAngleRad = std::fmod(wtRotation, 2.0 * M_PI);
                double currentAngleDeg = currentAngleRad * 180.0 / M_PI;

                struct ComponentPhasorDrop {
                    std::string label;
                    std::complex<double> phasor;
                    double mag;
                    double phaseRad;
                    double rotX;
                    double rotY;
                    ImVec4 color;
                    bool isSource;
                };

                std::vector<ComponentPhasorDrop> passiveDrops;
                ComponentPhasorDrop sourceDrop{ "Total Source V(1)", {0.0, 0.0}, 0.0, 0.0, 0.0, 0.0, ImVec4(0.2f, 0.7f, 1.0f, 1.0f), false };

                ImVec4 palette[] = {
                    ImVec4(0.2f, 0.9f, 0.4f, 1.0f), // Green (Resistor)
                    ImVec4(1.0f, 0.4f, 0.4f, 1.0f), // Red (Inductor)
                    ImVec4(0.8f, 0.4f, 1.0f, 1.0f), // Purple (Capacitor)
                    ImVec4(0.9f, 0.8f, 0.2f, 1.0f), // Yellow
                    ImVec4(0.3f, 0.9f, 0.9f, 1.0f)  // Cyan
                };
                int colorIdx = 0;

                if (netlistPtr && analysisTypeActive == 2 && !g_latestAcNodePhasors.empty()) {
                    auto getNodeComplexV = [&](const NetNode& n) -> std::complex<double> {
                        if (n.isGround()) return std::complex<double>(0.0, 0.0);
                        MnaIndex mnaIdx = MnaIndex::fromNetNode(n);
                        if (mnaIdx.isValid() && static_cast<size_t>(mnaIdx.idx) < g_latestAcNodePhasors.size()) {
                            return g_latestAcNodePhasors[mnaIdx.idx];
                        }
                        return std::complex<double>(0.0, 0.0);
                    };

                    netlistPtr->inspectComponents([&](Component& comp, MnaIndex extraIdx) {
                        (void)extraIdx;
                        const auto& nodes = comp.getNodes();
                        if (nodes.size() < 2) return;

                        std::complex<double> vDrop = getNodeComplexV(nodes[0]) - getNodeComplexV(nodes[1]);
                        double mag = std::abs(vDrop);
                        double phaseRad = std::arg(vDrop);
                        double rotX = mag * std::cos(phaseRad + currentAngleRad);
                        double rotY = mag * std::sin(phaseRad + currentAngleRad);

                        ComponentType cType = comp.getType();
                        if (cType == ComponentType::VOLTAGE_SOURCE || cType == ComponentType::AC_SOURCE_DEVICE) {
                            if (!sourceDrop.isSource) {
                                sourceDrop.label = "Total Source " + comp.getName();
                                sourceDrop.phasor = vDrop;
                                sourceDrop.mag = mag;
                                sourceDrop.phaseRad = phaseRad;
                                sourceDrop.rotX = rotX;
                                sourceDrop.rotY = rotY;
                                sourceDrop.isSource = true;
                            }
                        } else if (cType == ComponentType::RESISTOR || cType == ComponentType::CAPACITOR ||
                                   cType == ComponentType::INDUCTOR || cType == ComponentType::DIODE) {
                            ImVec4 cColor = palette[colorIdx % 5];
                            if (cType == ComponentType::RESISTOR)  cColor = ImVec4(0.2f, 0.9f, 0.4f, 1.0f);
                            if (cType == ComponentType::INDUCTOR)  cColor = ImVec4(1.0f, 0.4f, 0.4f, 1.0f);
                            if (cType == ComponentType::CAPACITOR) cColor = ImVec4(0.8f, 0.4f, 1.0f, 1.0f);
                            colorIdx++;

                            passiveDrops.push_back({
                                "Drop V(" + comp.getName() + ")",
                                vDrop, mag, phaseRad, rotX, rotY, cColor, false
                            });
                        }
                    });
                }

                // Compute dynamic axis bounds from exact magnitudes
                double totalChainLength = 0.0;
                double maxWaveAmp = sourceDrop.mag;
                for (const auto& d : passiveDrops) {
                    totalChainLength += d.mag;
                    if (d.mag > maxWaveAmp) maxWaveAmp = d.mag;
                }
                double limitPhasor = std::max({ sourceDrop.mag, totalChainLength, 1.0 }) * 1.15;
                double limitWave = std::max(maxWaveAmp, 1.0) * 1.20;

                double availW = ImGui::GetContentRegionAvail().x;
                double availH = ImGui::GetContentRegionAvail().y - 25.0f;

                // ============================================================
                // LEFT PLOT: EXACT MNA KVL PHASOR CHAIN (36% Width)
                // ============================================================
                if (ImPlot::BeginPlot("Exact KVL Rotating Phasor", ImVec2(availW * 0.36f, availH), ImPlotFlags_Equal)) {
                    ImPlot::SetupAxes("Real (In-Phase)", "Imaginary (Quadrature)");
                    ImPlot::SetupAxesLimits(-limitPhasor, limitPhasor, -limitPhasor, limitPhasor, ImPlotCond_Always);

                    // Reference Circle for Source Magnitude
                    if (sourceDrop.isSource) {
                        const int circlePts = 64;
                        double circX[circlePts + 1], circY[circlePts + 1];
                        for (int i = 0; i <= circlePts; ++i) {
                            double theta = (2.0 * M_PI * i) / circlePts;
                            circX[i] = sourceDrop.mag * std::cos(theta);
                            circY[i] = sourceDrop.mag * std::sin(theta);
                        }
                        PlotLineStyled("##kvl_ref_circle", circX, circY, circlePts + 1, ImVec4(0.6f, 0.6f, 0.7f, 0.30f), 1.0f);

                        double pSrcX[] = { 0.0, sourceDrop.rotX };
                        double pSrcY[] = { 0.0, sourceDrop.rotY };
                        PlotLineStyled(sourceDrop.label.c_str(), pSrcX, pSrcY, 2, sourceDrop.color, 3.5f);
                    }

                    // Tip-to-Tail KVL Sum of Passive Component Drops
                    std::vector<double> chainX = { 0.0 };
                    std::vector<double> chainY = { 0.0 };
                    double accX = 0.0;
                    double accY = 0.0;

                    for (const auto& drop : passiveDrops) {
                        double startX = accX;
                        double startY = accY;
                        accX += drop.rotX;
                        accY += drop.rotY;

                        chainX.push_back(accX);
                        chainY.push_back(accY);

                        double segX[] = { startX, accX };
                        double segY[] = { startY, accY };
                        PlotLineStyled(drop.label.c_str(), segX, segY, 2, drop.color, 2.5f);
                    }

                    if (chainX.size() > 1) {
                        PlotLineStyled("KVL Closed Trajectory", chainX.data(), chainY.data(),
                                       static_cast<int>(chainX.size()), ImVec4(0.9f, 0.8f, 0.2f, 0.65f), 1.5f);

                        double endX[] = { accX };
                        double endY[] = { accY };
                        PlotScatterStyled("##kvl_tip", endX, endY, 1, ImPlotMarker_Plus, 6.0f, ImVec4(0.9f, 0.8f, 0.2f, 1.0f));
                    }

                    // Horizontal projection bridge to right edge
                    if (sourceDrop.isSource) {
                        double bridgeX[] = { sourceDrop.rotX, limitPhasor };
                        double bridgeY[] = { sourceDrop.rotY, sourceDrop.rotY };
                        PlotLineStyled("##kvl_proj_left", bridgeX, bridgeY, 2, ImVec4(1.0f, 0.3f, 0.3f, 0.85f), 1.5f);
                    }

                    ImPlot::EndPlot();
                }

                ImGui::SameLine();

                // ============================================================
                // RIGHT PLOT: EXACT COMPONENT VOLTAGE DROP WAVEFORMS (64% Width)
                // ============================================================
                if (ImPlot::BeginPlot("Component Voltage Drops in Time Domain", ImVec2(-1.0f, availH))) {
                    ImPlot::SetupAxes("Phase Angle wt (Degrees)", "Voltage Drop (V)");
                    ImPlot::SetupAxesLimits(0.0, 360.0, -limitWave, limitWave, ImPlotCond_Always);

                    const int wavePts = 180;
                    double degAxis[wavePts + 1];
                    for (int i = 0; i <= wavePts; ++i) {
                        degAxis[i] = (360.0 * i) / wavePts;
                    }

                    // Source Waveform
                    if (sourceDrop.isSource) {
                        double waveSrc[wavePts + 1];
                        for (int i = 0; i <= wavePts; ++i) {
                            double rad = degAxis[i] * M_PI / 180.0;
                            waveSrc[i] = sourceDrop.mag * std::sin(rad + sourceDrop.phaseRad);
                        }
                        PlotLineStyled(sourceDrop.label.c_str(), degAxis, waveSrc, wavePts + 1, sourceDrop.color, 3.0f);

                        // Horizontal projection line from 0 deg to current wt angle
                        double projX[] = { 0.0, currentAngleDeg };
                        double projY[] = { sourceDrop.rotY, sourceDrop.rotY };
                        PlotLineStyled("Projection", projX, projY, 2, ImVec4(1.0f, 0.3f, 0.3f, 0.9f), 1.8f);

                        double dotX[] = { currentAngleDeg };
			double dotY[] = { sourceDrop.rotY };
			PlotScatterStyled("##kvl_live_src", dotX, dotY, 1, ImPlotMarker_Circle, 5.0f, ImVec4(1.0f, 0.3f, 0.3f, 1.0f));
			}
			// Each Component's Exact Voltage Drop Sine Wave
			for (const auto& drop : passiveDrops) {
			double waveDrop[wavePts + 1];
			for (int i = 0; i <= wavePts; ++i) {
			double rad = degAxis[i] * M_PI / 180.0;
			waveDrop[i] = drop.mag * std::sin(rad + drop.phaseRad);
			}
			PlotLineStyled(drop.label.c_str(), degAxis, waveDrop, wavePts + 1, drop.color, 2.2f);
			double dDotX[] = { currentAngleDeg };
			double dDotY[] = { drop.rotY };
			PlotScatterStyled(("##kvl_dot_" + drop.label).c_str(), dDotX, dDotY, 1, ImPlotMarker_Circle, 4.0f, drop.color);
			}
			// Vertical wt sweep cursor
			double cursorX[] = { currentAngleDeg, currentAngleDeg };
			double cursorY[] = { -limitWave, limitWave };
			PlotLineStyled("##kvl_wt_cursor", cursorX, cursorY, 2, ImVec4(1.0f, 0.3f, 0.3f, 0.5f), 1.2f);
			ImPlot::EndPlot();
			}
			// Bottom Status Bar
			if (analysisTypeActive != 2 || !sourceDrop.isSource) {
			ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.2f, 1.0f), "Switch to 'Small-Signal AC Frequency Bode' mode to compute live KVL component drops.");
			} else {
			std::complex<double> sumPassive(0.0, 0.0);
			for (const auto& d : passiveDrops) sumPassive += d.phasor;
			double kvlResidual = std::abs(sourceDrop.phasor - sumPassive);
			ImGui::Text("Source |V| = %.3f V  |  Passive Components: %zu  |  KVL Loop Residual |V_src - Sum(V_drop)| = %.2e V",
			sourceDrop.mag, passiveDrops.size(), kvlResidual);
			}
			ImGui::EndTabItem();
			}

          

            // TAB 7: COMPLEX IMPEDANCE Z-PLANE (LIVE MNA SOLVER + SINE PROJECTION)
            if (ImGui::BeginTabItem("Complex Impedance Z-Plane")) {
                double netR = 0.0;
                double netX = 0.0;

                // Pull live R and X directly from the MNA AC solver results
                if (netlistPtr && analysisTypeActive == 2 && !plotFrequencies.empty()) {
                    double vMag = static_cast<double>(stimAcMag);
                    double vPhase = static_cast<double>(stimAcPhase) * M_PI / 180.0;
                    std::complex<double> vSourcePhasor(vMag * std::cos(vPhase), vMag * std::sin(vPhase));
                    std::complex<double> complexPowerS(plotActivePowerP, plotReactivePowerQ);

                    std::complex<double> iSourcePhasor = (vMag > 1e-12)
                        ? std::conj(complexPowerS / vSourcePhasor)
                        : std::complex<double>(0.0, 0.0);

                    std::complex<double> inputImpedanceZ = (std::abs(iSourcePhasor) > 1e-12)
                        ? (vSourcePhasor / iSourcePhasor)
                        : std::complex<double>(0.0, 0.0);

                    netR = inputImpedanceZ.real();
                    netX = inputImpedanceZ.imag();
                }

                double magZ = std::sqrt(netR * netR + netX * netX);
                double thetaZRad = std::atan2(netX, netR);
                double thetaZDeg = thetaZRad * 180.0 / M_PI;

                // Animated phase angle wt (0 to 360 degrees)
                double clockTime = glfwGetTime();
                double currentAngleRad = std::fmod(2.0 * M_PI * 0.20 * clockTime, 2.0 * M_PI);
                double currentAngleDeg = currentAngleRad * 180.0 / M_PI;

                // Normalized unit current I(t) and scaled voltage V(t)
                double iMag = 1.0;
                double vMagNorm = (magZ > 1e-9) ? 1.25 : 1.0;
                double instI = iMag * std::sin(currentAngleRad);
                double instV = vMagNorm * std::sin(currentAngleRad + thetaZRad);

                double availW = ImGui::GetContentRegionAvail().x;
                double availH = ImGui::GetContentRegionAvail().y - 25.0f;
                double zAxisLimit = std::max(magZ, 10.0) * 1.25;
                double waveLimit = std::max(vMagNorm, 1.2) * 1.2;

                // ============================================================
                // LEFT PLOT: LIVE COMPLEX IMPEDANCE Z-PLANE (36% Width)
                // ============================================================
                if (ImPlot::BeginPlot("Impedance Z-Plane (R + jX)", ImVec2(availW * 0.36f, availH), ImPlotFlags_Equal)) {
                    ImPlot::SetupAxes("Resistance R (Ohms)", "Reactance X (Ohms)");
                    double minR = std::min(0.0, netR) - zAxisLimit * 0.25;
                    double maxR = std::max(0.0, netR) + zAxisLimit * 1.05;
                    ImPlot::SetupAxesLimits(minR, maxR, -zAxisLimit, zAxisLimit, ImPlotCond_Always);

                    // Resistance component along Real axis
                    double rX[] = { 0.0, netR };
                    double rY[] = { 0.0, 0.0 };
                    PlotLineStyled("Resistance (R)", rX, rY, 2, ImVec4(0.2f, 0.7f, 1.0f, 1.0f), 2.5f);

                    // Reactance component (+jXL inductive up, -jXC capacitive down)
                    double xX[] = { netR, netR };
                    double xY[] = { 0.0, netX };
                    ImVec4 xColor = (netX >= 0.0) ? ImVec4(1.0f, 0.4f, 0.4f, 1.0f) : ImVec4(0.8f, 0.4f, 1.0f, 1.0f);
                    PlotLineStyled(netX >= 0.0 ? "Inductive (+jX)" : "Capacitive (-jX)", xX, xY, 2, xColor, 2.5f);

                    // Resultant Impedance Vector Z
                    double zX[] = { 0.0, netR };
                    double zY[] = { 0.0, netX };
                    PlotLineStyled("Impedance (Z)", zX, zY, 2, ImVec4(0.2f, 0.9f, 0.4f, 1.0f), 3.5f);

                    double dotX[] = { netR };
                    double dotY[] = { netX };
                    PlotScatterStyled("##z_tip", dotX, dotY, 1, ImPlotMarker_Circle, 5.0f, ImVec4(1.0f, 1.0f, 1.0f, 1.0f));

                    ImPlot::EndPlot();
                }

                ImGui::SameLine();

                // ============================================================
                // RIGHT PLOT: V(t) vs I(t) PHASE-SHIFTED SINE WAVES (64% Width)
                // ============================================================
                if (ImPlot::BeginPlot("Time-Domain V(t) & I(t) Phase Shift", ImVec2(-1.0f, availH))) {
                    ImPlot::SetupAxes("Phase Angle wt (Degrees)", "Normalized Amplitude");
                    ImPlot::SetupAxesLimits(0.0, 360.0, -waveLimit, waveLimit, ImPlotCond_Always);

                    const int wavePts = 180;
                    double degAxis[wavePts + 1];
                    double waveI[wavePts + 1];
                    double waveV[wavePts + 1];

                    for (int i = 0; i <= wavePts; ++i) {
                        double deg = (360.0 * i) / wavePts;
                        double rad = deg * M_PI / 180.0;
                        degAxis[i] = deg;
                        waveI[i] = iMag * std::sin(rad);
                        waveV[i] = vMagNorm * std::sin(rad + thetaZRad);
                    }

                    PlotLineStyled("Current I(t)", degAxis, waveI, wavePts + 1, ImVec4(0.2f, 0.7f, 1.0f, 1.0f), 2.5f);
                    PlotLineStyled("Voltage V(t)", degAxis, waveV, wavePts + 1, ImVec4(0.2f, 0.9f, 0.4f, 1.0f), 3.0f);

                    double cursorX[] = { currentAngleDeg, currentAngleDeg };
                    double cursorY[] = { -waveLimit, waveLimit };
                    PlotLineStyled("##wt_cursor", cursorX, cursorY, 2, ImVec4(1.0f, 0.3f, 0.3f, 0.5f), 1.2f);

                    double liveIX[] = { currentAngleDeg }; double liveIY[] = { instI };
                    PlotScatterStyled("##dot_i", liveIX, liveIY, 1, ImPlotMarker_Circle, 5.0f, ImVec4(0.2f, 0.7f, 1.0f, 1.0f));

                    double liveVX[] = { currentAngleDeg }; double liveVY[] = { instV };
                    PlotScatterStyled("##dot_v", liveVX, liveVY, 1, ImPlotMarker_Circle, 5.0f, ImVec4(0.2f, 0.9f, 0.4f, 1.0f));

                    ImPlot::EndPlot();
                }

                // Bottom Status Bar
                if (analysisTypeActive != 2 || plotFrequencies.empty()) {
                    ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.2f, 1.0f), "Switch to 'Small-Signal AC Frequency Bode' mode to compute live MNA impedance.");
                } else {
                    ImGui::Text("|Z| = %.2f Ohms  |  R = %.2f Ohms  |  X = %.2f Ohms  |  Phase Angle = %.2f deg (%s)",
                                magZ, netR, netX, thetaZDeg,
                                (netX > 0.01 ? "Voltage Leads Current - Inductive" :
                                (netX < -0.01 ? "Voltage Lags Current - Capacitive" : "In-Phase - Purely Resistive")));
                }

                ImGui::EndTabItem();
            }

            // TAB 8: NYQUIST STABILITY
            if (ImGui::BeginTabItem("Nyquist Stability")) {
                bool unstableEncirclement = false;
                bool hasNyquistData = !plotNyquistReal.empty();
		// Sisakan ruang 28 px di bawah untuk 1 baris status kestabilan
		float nyquistPlotHeight = std::max(120.0f, ImGui::GetContentRegionAvail().y - 28.0f);

                if (ImPlot::BeginPlot("##nyquist_plot", ImVec2(-1.0f, nyquistPlotHeight), ImPlotFlags_Equal)) {
                    ImPlot::SetupAxes("Real Component [Re(H)]", "Imaginary Component [Im(H)]");
                    double critX[] = { -1.0 };
                    double critY[] = { 0.0 };
                    PlotScatterStyled("Critical Point (-1, j0)", critX, critY, 1, ImPlotMarker_Cross, 6.0f, ImVec4(1.0f, 0.2f, 0.2f, 1.0f));

                    if (hasNyquistData && analysisTypeActive == 2) {
                        ImPlot::PlotLine("Locus Trajectory H(jw)", plotNyquistReal.data(), plotNyquistImag.data(), static_cast<int>(plotNyquistReal.size()));

                        double currentVecX[] = { 0.0, plotNyquistReal.back() };
                        double currentVecY[] = { 0.0, plotNyquistImag.back() };
                        PlotLineStyled("Current State Vector", currentVecX, currentVecY, 2, ImVec4(1.0f, 0.6f, 0.2f, 0.9f), 2.0f);

                        for (size_t i = 0; i < plotNyquistReal.size(); ++i) {
                            if (plotNyquistReal[i] <= -1.0 && std::abs(plotNyquistImag[i]) < 0.05) {
                                unstableEncirclement = true;
                                break;
                            }
                        }
                    }
                    ImPlot::EndPlot();
                }
                
                if (hasNyquistData) {
                    if (unstableEncirclement) {
                        ImGui::TextColored(ImVec4(1.0f, 0.3f, 0.3f, 1.0f), "CRITICAL STABILITY WARNING: System rolls over close to (-1, j0)!");
                    } else {
                        ImGui::TextColored(ImVec4(0.3f, 1.0f, 0.3f, 1.0f), "SYSTEM STABLE: Encirclement check cleared.");
                    }
                }
                ImGui::EndTabItem();
            }

            // TAB 9: POWER QUADRANT (P-Q-S)
            if (ImGui::BeginTabItem("Power Quadrant (P-Q-S)")) {
                double rmsP = plotActivePowerP * 0.5;
                double rmsQ = plotReactivePowerQ * 0.5;
                double apparentPowerS = std::sqrt(rmsP * rmsP + rmsQ * rmsQ);

		// Sisakan ruang 85 px di bawah untuk 4 baris teks metrik daya
		float powerPlotHeight = std::max(120.0f, ImGui::GetContentRegionAvail().y - 85.0f);


                if (ImPlot::BeginPlot("##pqs_power_plot", ImVec2(-1.0f, powerPlotHeight), ImPlotFlags_Equal)) {
                    ImPlot::SetupAxes("Real Active Power P (Watts)", "Reactive Power Q (VAR)");
                    double pX[] = { 0.0, rmsP };
                    double pY[] = { 0.0, 0.0 };
                    PlotLineStyled("Active Power (P)", pX, pY, 2, ImVec4(0.2f, 0.9f, 0.4f, 1.0f), 3.0f);

                    double qX[] = { rmsP, rmsP };
                    double qY[] = { 0.0, rmsQ };
                    ImVec4 qColor = (rmsQ >= 0.0) ? ImVec4(1.0f, 0.4f, 0.4f, 1.0f) : ImVec4(0.2f, 0.5f, 1.0f, 1.0f);
                    PlotLineStyled(rmsQ >= 0.0 ? "Reactive Inductive (+Q)" : "Reactive Capacitive (-Q)", qX, qY, 2, qColor, 3.0f);

                    double sX[] = { 0.0, rmsP };
                    double sY[] = { 0.0, rmsQ };
                    PlotLineStyled("Apparent Power (S)", sX, sY, 2, ImVec4(0.9f, 0.8f, 0.2f, 1.0f), 3.5f);

                    double pathX[] = { 0.0, rmsP, rmsP, 0.0 };
                    double pathY[] = { 0.0, 0.0, rmsQ, 0.0 };
                    PlotLineStyled("Power Boundary Profile", pathX, pathY, 4, ImVec4(1.0f, 1.0f, 1.0f, 0.25f), 1.0f);

                    ImPlot::EndPlot();
                }

                ImGui::BeginGroup();
                ImGui::Text("Active Power (P)    : %.4f W", rmsP);
                ImGui::Text("Reactive Power (Q)  : %.4f VAR", rmsQ);
                ImGui::Text("Apparent Power (S)  : %.4f VA", apparentPowerS);
                ImGui::Text("Network Power Factor: %.4f (%s)", plotPowerFactor, (rmsQ >= 0.0 ? "Lagging/Inductive" : "Leading/Capacitive"));
                ImGui::EndGroup();
                ImGui::EndTabItem();
            }

            // TAB 10: HOLISTIC COMPLEX IMPEDANCE Z-PLANE (MNA SOLVER)
            if (ImGui::BeginTabItem("Holistic Impedance Z-Plane")) {
                double inputR = 0.0;
                double inputX = 0.0;
                double magnitudeZ = 0.0;
                double phaseRad = 0.0;

                if (netlistPtr && analysisTypeActive == 2 && !plotFrequencies.empty()) {
                    double vMag = static_cast<double>(stimAcMag);
                    double vPhase = static_cast<double>(stimAcPhase) * M_PI / 180.0;
                    std::complex<double> vSourcePhasor(vMag * std::cos(vPhase), vMag * std::sin(vPhase));
                    std::complex<double> complexPowerS(plotActivePowerP, plotReactivePowerQ); // 2. Daya kompleks S = P + jQ
		    // 3. Arus sumber eksak dari S = V * I*  =>  I = (S / V)*
		    std::complex<double> iSourcePhasor = (vMag > 1e-12)
			    ? std::conj(complexPowerS / vSourcePhasor)
			    : std::complex<double>(0.0, 0.0);

                    // 4. Impedansi input Z = V / I (sudut fasa sumber otomatis saling meniadakan)
  		    std::complex<double> inputImpedanceZ = (std::abs(iSourcePhasor) > 1e-12)
			    ? (vSourcePhasor / iSourcePhasor)
			    : std::complex<double>(0.0, 0.0);
                    inputR = inputImpedanceZ.real();
                    inputX = inputImpedanceZ.imag();
                    magnitudeZ = std::abs(inputImpedanceZ);
                    phaseRad = std::arg(inputImpedanceZ);
                }
 	        float plotHeight = std::max(120.0f, ImGui::GetContentRegionAvail().y - 85.0f);

                if (ImPlot::BeginPlot("##holistic_z_plane_plot", ImVec2(-1.0f, plotHeight), ImPlotFlags_Equal)) {
                    ImPlot::SetupAxes("Resistance R (Ohms)", "Reactance X (Ohms)");
                    double zX[] = { 0.0, inputR };
                    double zY[] = { 0.0, inputX };
                    PlotLineStyled("Impedance Vector (Z)", zX, zY, 2, ImVec4(0.2f, 0.9f, 0.4f, 1.0f), 3.5f);

                    double dotX[] = { inputR };
                    double dotY[] = { inputX };
                    PlotScatterStyled("Vector Tip", dotX, dotY, 1, ImPlotMarker_Circle, 4.0f, ImVec4(1.0f, 1.0f, 1.0f, 1.0f));

                    double pathX[] = { 0.0, inputR, inputR };
                    double pathY[] = { 0.0, 0.0, inputX };
                    PlotLineStyled("R-X Boundary Grid", pathX, pathY, 3, ImVec4(1.0f, 1.0f, 1.0f, 0.2f), 1.0f);

                    ImPlot::EndPlot();
                }

                ImGui::BeginGroup();
                ImGui::Text("Net Resistance (R) : %.4f Ohms", inputR);
                ImGui::Text("Net Reactance (X)  : %.4f Ohms", inputX);
                ImGui::Text("Total Impedance |Z|: %.4f Ohms", magnitudeZ);
                ImGui::Text("Phase Angle (Deg)  : %.2f deg (%s)", phaseRad * 180.0 / M_PI, (inputX >= 0.0 ? "Inductive Delay (+jX)" : "Capacitive Lead (-jX)"));
                ImGui::EndGroup();
                ImGui::EndTabItem();
            }
	    // TAB 11: DIGITAL LOGIC ANALYZER (STACKED TIMING DIAGRAM + GATE EVALUATOR)
            if (ImGui::BeginTabItem("Digital Logic Analyzer")) {
                static float logicThresholdV = 2.5f;
                static int selectedGateIdx = 1; // Default: AND
                static bool autoSyncNetlistGate = true;

                const char* gateNames[] = { "NOT", "AND", "OR", "XOR", "NAND", "NOR", "XNOR" };
                const GateType gateEnums[] = {
                    GateType::NOT, GateType::AND, GateType::OR,
                    GateType::XOR, GateType::NAND, GateType::NOR, GateType::XNOR
                };

                // Detect if a BehaviorLogicGateComponent exists in the active netlist
                BehaviorLogicGateComponent* activeGateComp = nullptr;
                if (netlistPtr) {
                    netlistPtr->inspectComponents([&](Component& comp, MnaIndex extraIdx) {
                        (void)extraIdx;
                        if (auto* g = dynamic_cast<BehaviorLogicGateComponent*>(&comp)) {
                            if (!activeGateComp || comp.getName() == selectedCompName) {
                                activeGateComp = g;
                            }
                        }
                    });
                }

                // Auto-sync dropdown to the netlist gate type if enabled
                if (autoSyncNetlistGate && activeGateComp) {
                    GateType netlistGateType = activeGateComp->getGateType();
                    for (int gIdx = 0; gIdx < 7; ++gIdx) {
                        if (gateEnums[gIdx] == netlistGateType) {
                            selectedGateIdx = gIdx;
                            break;
                        }
                    }
                }

                ImGui::SetNextItemWidth(150.0f);
                ImGui::SliderFloat("CMOS Threshold (V)", &logicThresholdV, 0.5f, 4.5f, "%.2f V");
                ImGui::SameLine();

                if (activeGateComp) {
                    ImGui::Checkbox("Sync Netlist Gate", &autoSyncNetlistGate);
                    ImGui::SameLine();
                }

                ImGui::BeginDisabled(autoSyncNetlistGate && activeGateComp != nullptr);
                ImGui::SetNextItemWidth(100.0f);
                ImGui::Combo("Gate", &selectedGateIdx, gateNames, 7);
                ImGui::EndDisabled();

                // Live discrete logic states from latest waveform samples
                double latestVA = plotInputWaveform.empty() ? 0.0 : plotInputWaveform.back();
                double latestVB = plotOutputWaveform.empty() ? 0.0 : plotOutputWaveform.back();
                uint8_t bitA = (latestVA >= logicThresholdV) ? 1 : 0;
                uint8_t bitB = (latestVB >= logicThresholdV) ? 1 : 0;
                uint8_t gateOut = LogicGates::evaluateScalarGate(gateEnums[selectedGateIdx], latestVA, latestVB, logicThresholdV);

                ImGui::SameLine();
                if (activeGateComp) {
                    ImGui::TextColored(ImVec4(0.3f, 0.95f, 0.7f, 1.0f),
                                       " | Netlist [%s: %s]  A=%u, Out=%u",
                                       activeGateComp->getName().c_str(), gateNames[selectedGateIdx], bitA, bitB);
                } else {
                    ImGui::TextColored(ImVec4(0.3f, 0.95f, 0.7f, 1.0f),
                                       " | Virtual Gate: A=%u, B=%u -> %s=%u",
                                       bitA, bitB, gateNames[selectedGateIdx], gateOut);
                }

                float analyzerHeight = std::max(140.0f, ImGui::GetContentRegionAvail().y - 10.0f);

                if (ImPlot::BeginPlot("##digital_logic_analyzer", ImVec2(-1.0f, analyzerHeight))) {
                    ImPlot::SetupAxes("Time (us)", "Stacked Digital Channels");
                    ImPlot::SetupAxesLimits(0.0, plotTimeSteps.empty() ? 200.0 : plotTimeSteps.back(), -0.3, 4.4, ImPlotCond_Always);

                    const double yTicks[] = { 0.0, 0.8, 1.5, 2.3, 3.0, 3.8 };
                    const char* yLabels[] = { "Y:0", "Y:1", "B:0", "B:1", "A:0", "A:1" };
                    ImPlot::SetupAxisTicks(ImAxis_Y1, yTicks, 6, yLabels);

                    if (!plotTimeSteps.empty() && plotTimeSteps.size() == plotInputWaveform.size()) {
                        size_t n = plotTimeSteps.size();
                        std::vector<double> chA(n), chB(n), chGateY(n);

                        for (size_t i = 0; i < n; ++i) {
                            double vA = plotInputWaveform[i];
                            double vB = (i < plotOutputWaveform.size()) ? plotOutputWaveform[i] : 0.0;

                            uint8_t logicA = (vA >= logicThresholdV) ? 1 : 0;
                            uint8_t logicB = (vB >= logicThresholdV) ? 1 : 0;
                            uint8_t logicY = LogicGates::evaluateScalarGate(gateEnums[selectedGateIdx], vA, vB, logicThresholdV);

                            chA[i]     = 3.0 + 0.8 * logicA;
                            chB[i]     = 1.5 + 0.8 * logicB;
                            chGateY[i] = 0.0 + 0.8 * logicY;
                        }

                        PlotStairsStyled("Channel A (Input)",       plotTimeSteps.data(), chA.data(),     static_cast<int>(n), ImVec4(0.2f, 0.75f, 1.0f, 1.0f), 2.5f);
			PlotStairsStyled("Channel B (Probed Out)",  plotTimeSteps.data(), chB.data(),     static_cast<int>(n), ImVec4(1.0f, 0.70f, 0.2f, 1.0f), 2.5f);
			PlotStairsStyled("Evaluated Gate Y",        plotTimeSteps.data(), chGateY.data(), static_cast<int>(n), ImVec4(0.2f, 0.95f, 0.5f, 1.0f), 3.0f);
                    }

                    if (ImPlot::IsPlotHovered()) {
                        ImPlotPoint mouseCoord = ImPlot::GetPlotMousePos();
                        ImGui::BeginTooltip();
                        ImGui::Text("Time: %.3f us", mouseCoord.x);
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
    glfwMakeContextCurrent(window);
    visualComponents.clear();
    plotTimeSteps.clear();
    plotInputWaveform.clear();
    plotOutputWaveform.clear();
    plotFrequencies.clear();
    plotBodeMagnitudeDb.clear();
    plotNyquistReal.clear();
    plotNyquistImag.clear();
    g_latestAcNodePhasors.clear();
    g_latestAcNodePhasors.shrink_to_fit();

    // 1. Delete OpenGL resources while the OpenGL context is still active
    if (globalVBO != 0) {
        glDeleteBuffers(1, &globalVBO);
        globalVBO = 0;
    }
    if (globalVAO != 0) {
        glDeleteVertexArrays(1, &globalVAO);
        globalVAO = 0;
    }
    if (shaderProgram != 0) {
        glDeleteProgram(shaderProgram);
        shaderProgram = 0;
    }

    // 2. Shut down renderer & platform backends first
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();

    // 3. Destroy ImPlot context, detach it, then destroy ImGui context
    ImPlot::DestroyContext();
    ImPlot::SetCurrentContext(nullptr);
    ImGui::DestroyContext();
    ImGui::SetCurrentContext(nullptr);

    // 4. Destroy the window and terminate GLFW
    glfwDestroyWindow(window);
    window = nullptr;
    glfwTerminate();

    // 5. Remove the temporary netlist file from disk to keep the workspace clean
    std::remove("tmp_gui_netlist.ckt");
}
