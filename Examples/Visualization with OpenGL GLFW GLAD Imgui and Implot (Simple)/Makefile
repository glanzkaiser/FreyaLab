# ============================================================================
# COMPILER AND DIRECTORY CONFIGURATION
# ============================================================================
CXX      := g++
CC       := gcc
CXXFLAGS := -std=c++17 -Wall -Wextra -O2
CFLAGS   := -O2
TARGET   := freyalab_dashboard

# Build artifact output directories
OBJ_DIR  := obj

# Paths to external dependency source folders
GLAD_DIR     := ./glad
FREYALAB_DIR := ./freyalab

# Include path directives
INCLUDES := -I. \
            -I$(FREYALAB_DIR) \
            -I$(GLAD_DIR) \

# ============================================================================
# OS SPECIFIC LINKER FLAGS DETECTION (Linux & macOS)
# ============================================================================
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S), Linux)
    LIBS := -lfreyalab -limplot -limgui -lglfw -lGL -ldl -lpthread
endif
ifeq ($(UNAME_S), Darwin) # macOS
    LIBS := -lglfw -framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo
    CXXFLAGS += -D__APPLE__
endif

# ============================================================================
# SOURCE FILES SEGREGATION
# ============================================================================
# Project Main Source Files
APP_SOURCES := main.cpp Visualization.cpp

# External Dependency Source Files
GLAD_SOURCES  := $(GLAD_DIR)/glad.c


# ============================================================================
# OBJECT FILE GENERATION MAPPING
# ============================================================================
OBJECTS := $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(APP_SOURCES)) \
           $(OBJ_DIR)/glad/glad.o

# ============================================================================
# COMPILATION RULES AND TARGET PIELINES
# ============================================================================
all: $(TARGET)

# Main Executable Linkage Stage
$(TARGET): $(OBJECTS)
	@echo "🔗 Linking system objects into target binary: $(TARGET)..."
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET) $(LIBS)
	@echo "✨ Build fully operational and complete."

# Rule: Compile Application Sources
$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Rule: Compile C-based GLAD Loader Functionality
$(OBJ_DIR)/glad/glad.o: $(GLAD_SOURCES)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Housekeeping Cleanup Pipelines
clean:
	@echo " Cleaning workspace object buffers..."
	rm -rf $(OBJ_DIR) $(TARGET)

.PHONY: all clean
