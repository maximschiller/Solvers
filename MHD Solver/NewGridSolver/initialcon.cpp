/* CALCULATES THE INITIAL CONDITIONS OF THE SYSTEM */
#include "Variables.h"
#include <math.h>
#include <iostream>

// Initial condition variables
FILE *ConvergeDiff = nullptr;
FILE *ConvergeVal = nullptr;
FILE *Converge = nullptr;
double RhoIn, VzIn, peIn, TeIn, phIn, ThIn; // declaring the inlet value variables
double Jback=0.0, Vappl = 10.0; //Initial guesses.

void setInitial(){

    // Call three functions to set convective and dissipative boundary conditions
    // calls from boundarycon.cpp and calcconv.cpp
    // ConvBound();
    // EvalParam();
    // DissBound();

    // Can use J by setting it to be the integer that loops through cell centers
    // note that cell center array counts from 0 to imax-1 and 0 to jmax-1 (always hae one less cell center than vertex)

    // also need to compute the values of Rc and Ra for the given grid
    // These are all the fluid properties at the inlet of the thruster
    for(J=Rc+1; J<Ra; J++)
    {
        aSq[J][0] = gamConst*p[J][0]/Rho[J][0];
        CmSq[J][0] = aSq[J][0]+(Bt[J][0]*Bt[J][0]
        /(Mu*Rho[J][0]));
        Cf[J][0] = sqrt(CmSq[J][0]);
        //Now, Cf = Cfz. So Cfz is no longer computed
        //Largest eigenvalue
        SpecRadZ[J][0] = fabs(Vz[J][0]+Cf[J][0]);
    }
    
    // determining the size of the grid in the y direction
    int rows = jmax/2;
    
    // Initial values at inlet
    // Setting the value for these variables at the midpoint between Rc and Ra
    RhoIn= Rho[Rc+rows/2][1];
    VzIn = Vz[Rc+rows/2][1];
    peIn = pe[Rc+rows/2][1];
    TeIn = Te[Rc+rows/2][1];
    phIn = ph[Rc+rows/2][1];
    ThIn = Th[Rc+rows/2][1];

    // Print convergence data for initial conditions
    
    ConvergeDiff= fopen("ConvergeDiff.DAT","wb");
    // // fprintf(ConvergeDiff, "%s \t %s \t %s \t %s \t %s \t %s \t %s \t
    // // %s \t %s \t %s \t %s \t %s \t %s \n",
    // // "numSteps", "R", "Z", "u0Dmax", "R", "Z", "u1Dmax",
    // // "R", "Z", "u2Dmax","R", "Z", "u3Dmax", "R", "Z",
    // // "u4Dmax","u0Davg","u1Davg","u2Davg","u3Davg",
    // // "u4Davg");

    // fprintf(ConvergeDiff, "%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n",
    //     "numSteps", "R", "Z", "u0Dmax", "u1Dmax", "u2Dmax", "u3Dmax", "u4Dmax",
    //     "R", "Z", "u0Dmax", "u1Dmax", "u2Dmax", "u3Dmax", "u4Dmax",
    //     "R", "Z", "u0Davg", "u1Davg", "u2Davg", "u3Davg", "u4Davg");


    ConvergeVal = fopen("ConvergeVal.DAT","wb");
    // // fprintf(ConvergeVal, "%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n",
    // // "numSteps", "R", "Z", "u0Max", "R", "Z", "u1Max",
    // // "R", "Z", "u2Max","R", "Z", "u3Max", "R", "Z",
    // // "u4Max","u0Avg","u1Avg","u2Avg","u3Avg","u4Avg");
    
    // fprintf(ConvergeVal, "%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n",
    //     "numSteps", "R", "Z", "u0Max", "u1Max", "u2Max", "u3Max", "u4Max",
    //     "R", "Z", "u0Max", "u1Max", "u2Max", "u3Max", "u4Max",
    //     "R", "Z", "u0Avg", "u1Avg", "u2Avg", "u3Avg", "u4Avg");


    Converge = fopen("Converge.DAT","wb");
    // // fprintf(Converge, "%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n",
    // // "numSteps", "R", "Z", "DensDmax", "R", "Z",
    // // "KEDmax", "R", "Z", "BtDmax", "R", "Z",
    // // "EtotDmax", "R", "Z", "EthDmax", "DensDavg",
    // // "KEDavg","BtDavg", "EtotDavg", "EthDavg");

    // fprintf(Converge, "%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n",
    //     "numSteps", "R", "Z", "DensDmax", "KEDmax", "BtDmax", "EtotDmax", "EthDmax",
    //     "R", "Z", "DensDmax", "KEDmax", "BtDmax", "EtotDmax", "EthDmax",
    //     "R", "Z", "DensDavg", "KEDavg", "BtDavg", "EtotDavg", "EthDavg");


}
