/* TAKES MAIN INPUTS OF THRUSTER PARAMETERS AND OUTPUTS SOLUTION FILES */
#include "Variables.h"
#include <fstream>
#include <iostream>
#include <string>


double mdot, JTot, TotTime; // variables being taken in by input file
// original was MassFlowRate - Jmax

void GetInput()
{
    //Accepting input from file...
    std::ifstream in("INPUT.DAT");
    in >> mdot >> JTot >> TotTime ;
}

//-------------------------------------
//-------------------------------------
//-------------------------------------
// Might need to put this on a seperate script so it doesnt mess with declaration order
void outputResults(){

}