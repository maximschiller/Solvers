//---------------------------------------------------------
//InitialConditions.CPP; Written by Kameshwaran Sankaran
//---------------------------------------------------------
#include "Thruster-Axi.h"
#include <math.h>

// Variable definitions
double peIn, phIn, TeIn, ThIn;

void SetInitial()
{
    /*
    //This is in case of starting after purely fluid flow.
    ifstream InRho("RhoInitial.DAT");
    ifstream InPress("PressInitial.DAT");
    ifstream InVr("VrInitial.DAT");
    ifstream InVz("VzInitial.DAT");
    for(Zone = 1; Zone <= 3; Zone++)
    {
    for(J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
    {
    for(K=Zmin[Zone]+1; K<=Zmax[Zone]; K++)
    {
    InRho>>Rho[J][K];
    InPress>>p[J][K];
    InVr>>Vr[J][K];
    InVz>>Vz[J][K];
    gam[J][K]= 5.0/3.0;
    pe[J][K] = 0.5*p[J][K];
    ph[J][K] = pe[J][K];
    Th[J][K] = ph[J][K]/((Rho[J][K]/mAr)*kBoltz);
    Te[J][K] = Th[J][K];
    }
    }
    }
    */
    // //This is in case of restarting after some period of MHD flow.
    // std::ifstream LoadNo("NoOut.DAT");
    // std::ifstream LoadNe("NeOut.DAT");
    // std::ifstream LoadTe("TeOut.DAT");
    // std::ifstream LoadTh("ThOut.DAT");
    // std::ifstream LoadPh("PhOut.DAT");
    // std::ifstream LoadGam("GamOut.DAT");
    // std::ifstream LoadVr("VrOut.DAT");
    // std::ifstream LoadVz("VzOut.DAT");
    // std::ifstream LoadBt("BtOut.DAT");

    // for(Zone = 1; Zone <= 3; Zone++)
    // {
    //     for(J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
    //     {
    //         for(K=Zmin[Zone]+1; K<=Zmax[Zone]; K++)
    //         {
    //             LoadNo>>n[J][K];
    //             LoadNe>>ne[J][K];
    //             LoadTe>>Te[J][K];
    //             LoadTh>>Th[J][K];
    //             LoadPh>>ph[J][K];
    //             LoadGam>>gam[J][K];
    //             LoadVr>>Vr[J][K];
    //             LoadVz>>Vz[J][K];
    //             LoadBt>>Bt[J][K];
    //             Rho[J][K] = n[J][K]*mAr;
    //             pe[J][K]= ne[J][K]*kBoltz*Te[J][K];
    //             VSq[J][K] = (Vr[J][K]*Vr[J][K])
    //             +(Vz[J][K]*Vz[J][K]);
    //             BSq[J][K] = Bt[J][K]*Bt[J][K];
    //             p[J][K] = pe[J][K] + ph[J][K];
    //             KE[J][K] = 0.5*Rho[J][K]*VSq[J][K];
    //             Etherm[J][K]= p[J][K]/(gam[J][K]-1);
    //             E[J][K] = Etherm[J][K]+KE[J][K]
    //             +(0.5*BSq[J][K]/Mu);
    //             E[J][K] = (p[J][K]/(gam[J][K]-1))
    //             +(0.5*Rho[J][K]*VSq[J][K])
    //             +(0.5*BSq[J][K]/Mu);
    //         }
    //     }
    // }

    // std::ifstream LoadVolt("Voltage.DAT");
    // LoadVolt>>Vappl;
    // Saha();
    // for(Zone = 1; Zone <= 3; Zone++)
    // {
    //     for(J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
    //     {
    //         for(K=Zmin[Zone]+1; K<=Zmax[Zone]; K++)
    //         {
    //             u[0][J][K] = Rho[J][K];
    //             u[1][J][K] = Rho[J][K]*Vr[J][K];
    //             u[2][J][K] = Rho[J][K]*Vz[J][K];
    //             u[3][J][K] = Bt[J][K];
    //             u[4][J][K] = E[J][K];
    //         }
    //     }
    // }

    ConvBound();
    EvalParam();
    DissBound();
    for(J=Rc+1; J<=Ra; J++)
    {
        aSq[J][0] = gamConst*p[J][0]/Rho[J][0];
        CmSq[J][0] = aSq[J][0]+(Bt[J][0]*Bt[J][0]
        /(Mu*Rho[J][0]));
        Cf[J][0] = sqrt(CmSq[J][0]);
        //Now, Cf = Cfz. So Cfz is no longer computed
        //Largest eigenvalue
        SpecRadZ[J][0] = fabs(Vz[J][0]+Cf[J][0]);
    }
    
    RhoIn= Rho[Rc+25][1];
    VzIn = Vz[Rc+25][1];
    peIn = pe[Rc+25][1];
    TeIn = Te[Rc+25][1];
    phIn = ph[Rc+25][1];
    ThIn = Th[Rc+25][1];

    ConvergeDiff= fopen("ConvergeDiff.DAT","wb");
    // fprintf(ConvergeDiff, "%s \t %s \t %s \t %s \t %s \t %s \t %s \t
    // %s \t %s \t %s \t %s \t %s \t %s \n",
    // "numSteps", "R", "Z", "u0Dmax", "R", "Z", "u1Dmax",
    // "R", "Z", "u2Dmax","R", "Z", "u3Dmax", "R", "Z",
    // "u4Dmax","u0Davg","u1Davg","u2Davg","u3Davg",
    // "u4Davg");

    fprintf(ConvergeDiff, "%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n",
        "numSteps", "R", "Z", "u0Dmax", "u1Dmax", "u2Dmax", "u3Dmax", "u4Dmax",
        "R", "Z", "u0Dmax", "u1Dmax", "u2Dmax", "u3Dmax", "u4Dmax",
        "R", "Z", "u0Davg", "u1Davg", "u2Davg", "u3Davg", "u4Davg");


    ConvergeVal = fopen("ConvergeVal.DAT","wb");
    // fprintf(ConvergeVal, "%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n",
    // "numSteps", "R", "Z", "u0Max", "R", "Z", "u1Max",
    // "R", "Z", "u2Max","R", "Z", "u3Max", "R", "Z",
    // "u4Max","u0Avg","u1Avg","u2Avg","u3Avg","u4Avg");
    
    fprintf(ConvergeVal, "%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n",
        "numSteps", "R", "Z", "u0Max", "u1Max", "u2Max", "u3Max", "u4Max",
        "R", "Z", "u0Max", "u1Max", "u2Max", "u3Max", "u4Max",
        "R", "Z", "u0Avg", "u1Avg", "u2Avg", "u3Avg", "u4Avg");


    Converge = fopen("Converge.DAT","wb");
    // fprintf(Converge, "%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n",
    // "numSteps", "R", "Z", "DensDmax", "R", "Z",
    // "KEDmax", "R", "Z", "BtDmax", "R", "Z",
    // "EtotDmax", "R", "Z", "EthDmax", "DensDavg",
    // "KEDavg","BtDavg", "EtotDavg", "EthDavg");

    fprintf(Converge, "%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n",
        "numSteps", "R", "Z", "DensDmax", "KEDmax", "BtDmax", "EtotDmax", "EthDmax",
        "R", "Z", "DensDmax", "KEDmax", "BtDmax", "EtotDmax", "EthDmax",
        "R", "Z", "DensDavg", "KEDavg", "BtDavg", "EtotDavg", "EthDavg");

}//End of function
//-------------------------------------