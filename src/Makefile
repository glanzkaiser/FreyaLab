CXX = g++

FREYALAB_DIR = ../src
LIB = libfreyalab.so 
SOURCES = CircuitComponents.cpp BSIM.cpp
SOURCES +=ACCircuit.cpp DCCircuit.cpp FileParser.cpp KLU.cpp math.cpp PulseCircuit.cpp 
SOURCES +=  BehavioralSource.cpp ElectroThermal.cpp NoiseCircuit.cpp LogicGates.cpp Netlist.cpp 


OBJS = $(addsuffix .o, $(basename $(notdir $(SOURCES))))
UNAME_S := $(shell uname -s)

CXXFLAGS = -fPIC
CXXFLAGS += -g -Wall -Wformat

##---------------------------------------------------------------------
## BUILD FLAGS PER PLATFORM
##---------------------------------------------------------------------

ifeq ($(UNAME_S), Linux) #LINUX
	ECHO_MESSAGE = "FreyaLab from GFreya OS"
	
	CFLAGS = $(CXXFLAGS)
endif
##---------------------------------------------------------------------
## BUILD RULES
##---------------------------------------------------------------------

%.o:$(SFREYALAB_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

all: $(LIB)
	@echo Build complete for $(ECHO_MESSAGE)

$(LIB): $(OBJS)
	$(CXX) -shared -o $@ $^ $(CXXFLAGS) 

clean:
	rm -f $(OBJS)
