#ifndef MY_HEADER_FILE_H
#define MY_HEADER_FILE_H

#include <iostream>

// IO.cpp variables - initial parameters
extern double MassFlowRate, Jmax, RAD, 
              Ranode, Rcathode, Z, AnodeLength, CathodeLength, 
              numR, numZ, TotTime;

// GridGen.cpp variables
const int Zones = 4;
extern double deltaR, deltaZ, Rc, Ra, Za, Zc, Rmin[Zones], Rmax[Zones], Zmin[Zones], Zmax[Zones];


// Functions

void GetInput();
void Grid();
void WriteVTK(std::string filename);

#endif