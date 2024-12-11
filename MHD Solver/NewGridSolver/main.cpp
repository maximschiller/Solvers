#include "Variables.h"
#include <iostream>

// NOTE for counters - i,I,x,k,Z are in thruster lengthwise direction
//                   - j,J,y,R are all in radial direction
// index accessing is done by (r,z) i.e. (j,i/k)
int i,j; // computational domain variables
int I,J,K; // declaring the counters for cell centres
int Imax, Jmax; // declaring the number of cell centres radially and thruster-length way
int imax=15, jmax=15; // Note that numZ = imax and numR = jmax - this sets the discretisation in the computational space indicating the number of vertices on the grid corners
// imax = numR+1 and jmax = numZ+1
int numR=jmax, numZ=imax; // setting the size of number of divisions radially and length-wise
int Zone;

double RAD = 1., Z = 1.; // domain sizing - need to update - COULD potentially set it as a function of the GMSH file input
double deltaR = RAD/jmax, deltaZ = Z/imax; // This sets the discretisation in the physical space depending on the domain size
double dx = deltaZ, dy = deltaR;

int numSteps=0, SubStep, MaxStep=10; // number of steps
int t; // time counter

// From linear algebra and eigensystem evaluation - used for matrix calculation
int row, col;
int l,k;

int main(){

    // Allocates Memory Spaces to Matrices
    MemAllocate();

    // Gets Input of Required Variables
    GetInput();

    // Domain and Mesh Generation / Evaluation
    DomainInput();
    discretiseGeom();
    gridgen();
        // grideval.cpp //
        cellCentres();
        cellGeom();
        dimensions();
        // gridError();

    // // Sets Initial and Boundary Conditions
    // setInitial();

    // // Solves Governing Equations
    // sysSolver();

    // // Outputs Results
    // outputResults();
}