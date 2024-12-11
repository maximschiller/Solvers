/* GENERATING THE MESH */
#include "Vars.h"
#include <math.h>

double deltaR, deltaZ, Rc, Ra, Za, Zc, Rmin[Zones], Rmax[Zones], Zmin[Zones], Zmax[Zones];

void Grid()
{
    /*
    In the cell-centered scheme, there is one face more than the center.
    In the channel:
    Faces: j = Rcathode, ... , Ra
    k = 0, ... , EL
    Centers:J = Rcathode+1, ... , Ra
    K = 1, ... , EL
    In the plume:
    Faces: j = 0, ... , numR
    k = EL, ..., numZ
    Centers:J = 1, ... , numR
    K = EL+1, ... , numZ
    */

    deltaR = RAD/double(numR);
    deltaZ = Z/double(numZ);

    //These 0.5s are added for the following reason:
    //to get a integer value of, for e.g, 3 from 2.99999999 and 3.00000001,
    //using int() would give 2 and 3 respectively.
    //To get 3 from both, add 0.5 to them and then int() them.
    Rc = int((Rcathode/deltaR)+0.5);
    Ra = int((Ranode/deltaR)+0.5);
    Za = int((AnodeLength/deltaZ)+0.5);
    Zc = int((CathodeLength/deltaZ)+0.5);

    //Corners of ZONE1
    Rmin[1] = Rc;
    Rmax[1] = Ra;
    Zmin[1] = 0;
    if(Zc <= Za)
    Zmax[1] = Zc;
    else
    Zmax[1] = Za;

    //Corners of ZONE2
    if(Zc <= Za)
    Rmin[2] = 0;
    else
    Rmin[2] = Rc;
    if(Zc <= Za)
    Rmax[2] = Ra;
    else
    Rmax[2] = numR;
    if(Zc <= Za)
    Zmin[2] = Zc;
    else
    Zmin[2] = Za;
    if(Zc <= Za)
    Zmax[2] = Za;
    else
    Zmax[2] = Zc;
    
    //Corners of ZONE3
    Rmin[3] = 0;
    Rmax[3] = numR;
    if(Zc <= Za)
    Zmin[3] = Za;
    else
    Zmin[3] = Zc;
    Zmax[3] = numZ;

    WriteVTK("output.vtk");
}
//-------------------------------------