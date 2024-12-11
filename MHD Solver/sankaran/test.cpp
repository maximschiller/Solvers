using namespace std;
#include <iostream>

int main()
{
    if (20 > 18) {
  cout << "20 is greater than 18";
}

}


thruster: InputOutput.o Grid.o Memory.o InitialConditions.o Solve.o SetTimeStep.o BoundaryConditions.o EvalConv.o EvalDiss.o TimeMarch.o ReCalculate.o
	g++ InputOutput.o Grid.o Memory.o InitialConditions.o Solve.o SetTimeStep.o BoundaryConditions.o EvalConv.o EvalDiss.o TimeMarch.o ReCalculate.o -o thruster

Thruster-Axi.o: Thruster-Axi.cpp
	g++ -c Thruster-Axi.cpp

InputOutput.o: InputOutput.cpp
	g++ -c InputOutput.cpp

Grid.o: Grid.cpp
	g++ -c Grid.cpp

Memory.o: Memory.cpp 
	g++ -c Memory.cpp

InitialConditions.o: InitialConditions.cpp 
	g++ -c InitialConditions.cpp

Solve.o: Solve.cpp 
	g++ -c Solve.cpp

SetTimeStep.o: SetTimeStep.cpp 
	g++ -c SetTimeStep.cpp

BoundaryConditions.o: BoundaryConditions.cpp
	g++ -c BoundaryConditions.cpp

EvalConv.o: EvalConv.cpp 
	g++ -c EvalConv.cpp

EvalDiss.o: EvalDiss.cpp 
	g++ -c EvalDiss.cpp

TimeMarch.o: TimeMarch.cpp 
	g++ -c TimeMarch.cpp

ReCalculate.o: ReCalculate.cpp 
	g++ -c ReCalculate.cpp

clean:
	rm *.o thruster



// #include "InputOutput.cpp"
// #include "Grid.cpp"
// #include "Memory.cpp"
// #include "InitialConditions.cpp"
// #include "Solve.cpp"
// #include "SetTimeStep.cpp"
// #include "BoundaryConditions.cpp"
// #include "EvalConv.cpp"
// #include "EvalDiss.cpp"
// #include "TimeMarch.cpp"
// #include "ReCalculate.cpp"
