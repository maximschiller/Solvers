CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=Thruster-Axi.cpp InputOutput.cpp Grid.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=Thruster-Axi

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf $(OBJECTS) $(EXECUTABLE)


g++ Thruster-Axi.cpp InputOutput.cpp Grid.cpp ReCalculate.cpp Memory.cpp InitialConditions.cpp BoundaryConditions.cpp Solve.cpp SetTimeStep.cpp EvalConv.cpp Jacobian.cpp
