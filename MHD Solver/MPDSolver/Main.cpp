// Outline of solver 

/*
1. Initial Parameters Declaraion
2. Geometry Creation
3. Grid Creation
4. Initial Conditions
5. Computation of Secondary Variables and Transport Coefficients
6. Setting Boundary Conditions
7. Estimating time step
8. Computing Fluxes
9. Updating conserved variables
10. Convergence checks
11. Solver loop of 5-10
12. Outputs
*/

#include "Vars.h"

int main()
{
    // Get geometry input in order to create the thruster and mesh
    // Scripts - IO.cpp - GridGen.cpp
    GetInput(); 
    // Allocate memory spaces for variables
    // Scripts - Memory.cpp
    MemAllocate();
    // Set the initial conditions for all the variables
    // Scripts - Initialise.cpp - Boundaries.cpp - ReCalcVar.cpp
    SetInitial();
    // Solve the governing equations
    // Scripts - TimeStep.cpp - ConvCalc.cpp (Jacobian - EigenSystem) - TimeMarch.cpp - ReCalcVar.cpp - EnergyCalc.cpp - LitCalc.cpp - ReCalcVar.cpp - IO.cpp - Jacobian.cpp - MatrixCalcs.cpp - EigenSystem.cpp
    Solve();

}