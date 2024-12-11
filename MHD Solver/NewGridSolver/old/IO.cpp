#include "Vars.h"
#include <fstream>
#include <iostream>
#include <string>

double MassFlowRate, Jmax, RAD, 
              Ranode, Rcathode, Z, AnodeLength, CathodeLength, 
              numR, numZ, TotTime;
int Zag = 1;

void GetInput()
{
    //Accepting input from file...
    std::ifstream in("INPUT.DAT");
    in >> MassFlowRate >> Jmax >> RAD >> Ranode >> Rcathode
    >> Z >> AnodeLength >> CathodeLength >> numR >> numZ
    >> TotTime ;
    Grid();
    std::cout << Zag << std::endl;
}


void WriteVTK(std::string filename)
{
    std::ofstream outfile(filename);
    
    // Write header
    outfile << "# vtk DataFile Version 3.0\n";
    outfile << "Paraview Output\n";
    outfile << "ASCII\n";
    outfile << "DATASET RECTILINEAR_GRID\n";
    
    // Write dimensions
    outfile << "DIMENSIONS " << numR+1 << " " << numZ+1 << " 1\n";
    
    // Write X coordinates
    outfile << "X_COORDINATES " << numR+1 << " float\n";
    for(int i = 0; i <= numR; i++)
    {
        double x = i*deltaR;
        outfile << x << " ";
    }
    outfile << "\n";
    
    // Write Y coordinates
    outfile << "Y_COORDINATES " << numZ+1 << " float\n";
    for(int j = 0; j <= numZ; j++)
    {
        double y = j*deltaZ;
        outfile << y << " ";
    }
    outfile << "\n";
    
    // Write Z coordinates
    outfile << "Z_COORDINATES 1 float\n";
    outfile << "0\n";
    
    // Write CELL_DATA
    outfile << "CELL_DATA " << numR*numZ << "\n";
    
    // Write ZONE data
    for(int j = 0; j < numZ; j++)
    {
        for(int i = 0; i < numR; i++)
        {
            // Determine which zone this cell is in
            int zone = 0;
            for(int z = 1; z <= Zones; z++)
            {
                if(i >= Rmin[z] && i <= Rmax[z] &&
                   j >= Zmin[z] && j <= Zmax[z])
                {
                    zone = z;
                    break;
                }
            }
            
            // Write zone number as scalar value
            outfile << zone << " ";
        }
        outfile << "\n";
    }
    
    outfile.close();
}
