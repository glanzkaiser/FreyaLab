CFLAGS = -ggdb
DEFINES = -DDEBUGGA
INCLUDES = 
LIBS = -lstdc++ -lsymintegration -lsfml-graphics -lsfml-window -lsfml-system
MAIN = main.o
CC=g++
# Source files
SRC=*.cpp LinearAlgebra/*.cpp GUI/*.cpp circuit_modules/*.cpp

freyalab::	
	$(CC) $(CFLAGS) $(SRC) $(INCLUDES) -o $@ $(EXE) $(LFLAGS) $(LIBS)

clean: 
	rm -f $(MAIN) main

