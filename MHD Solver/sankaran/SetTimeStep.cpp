//---------------------------------------------------------
//SetTimeStep.CPP; Written by Kameshwaran Sankaran
//---------------------------------------------------------
#include "Thruster-Axi.h"

// Variable Definitions
double Lmax, deltaTconv, deltaTtherm, deltaTresist;

void EvalTimeStep()
{
    //The time step is chosen to satisfy the CFL criterion.
    //For the hyperbolic part, deltaT = C * (GridSize/Max.wave speed)
    Lmax = SpecRadZ[Rc+25][0];//Initial guess.
    //SpecRadR and SpecRadZ are at centers. So, the face values are
    //obtained by averaging them.
    for(Zone = 1; Zone <= 3; Zone++)
    {                               
        for(j=J=Rmin[Zone]+1; j<=Rmax[Zone]-1; j++,J++)
        {
            for(K=Zmin[Zone]+1; K<=Zmax[Zone]; K++)
            {
                if (0.5*(SpecRadR[J][K]+SpecRadR[J+1][K]) > Lmax)
                {
                    Lmax = 0.5*(SpecRadR[J][K]+SpecRadR[J+1][K]);
                }
            }
        }
        
        for(J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
        {
            for(k=K=Zmin[Zone+1]; k<=Zmax[Zone]; k++,K++)
            {
                if (0.5*(SpecRadZ[J][K]+SpecRadZ[J][K+1]) > Lmax)
                {
                    Lmax = 0.5*(SpecRadZ[J][K]+SpecRadZ[J][K+1]);
                }
            }
        }
    }//end for all zones.

    deltaTconv = 0.6*deltaR/Lmax;
    deltaT = deltaTconv;
    //For thermal conduction, VNSA dictates that
    //deltaT <= 0.25*n[J][K]*kBoltz*deltaR*deltaR/kTherm[J][K].
    deltaTtherm = 5.0e-9;//Initial guess.
    for(Zone = 1; Zone <= 3; Zone++)
    {
        for (J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
        {
            for(K=Zmin[Zone]+1; K<=Zmax[Zone]; K++)
            {
                if(0.2*0.25*n[J][K]*kBoltz*deltaR*deltaR/kTherm[J][K]
                < deltaTtherm)
                    deltaTtherm = 0.2*0.25*n[J][K]*kBoltz*deltaR*deltaR
                    /kTherm[J][K];
            }
        }
    }
   
    if (deltaTtherm < deltaT)
    {
        deltaT = deltaTtherm;
    }
    
    //For resistive diffusion, VNSA dictates that
    //deltaT <= 0.25*Mu*deltaR*deltaR/Res[J][K].
    deltaTresist = 5.0e-10;//Initial guess.
    for(Zone = 1; Zone <= 3; Zone++)
    {
        for (J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
        {
            for(K=Zmin[Zone]+1; K<=Zmax[Zone]; K++)
            {
                if(0.2*0.25*Mu*deltaR*deltaR/Res[J][K] < deltaTresist)
                    deltaTresist =0.2*0.25*Mu*deltaR*deltaR/Res[J][K];
            }
        }
    }

    if (deltaTresist < deltaT)
    {
        deltaT = deltaTresist;
    }
    
    //This is if substepping is used.
    // if(deltaTtherm > deltaTresist)
    // MaxStep = int((deltaTconv/deltaTtherm)/10);
    // else
    // MaxStep = int((deltaTconv/deltaTresist)/10);
    // if(MaxStep < 1)
    MaxStep = 1;
    // if(MaxStep > 5)
    // MaxStep = 5;
}
//-------------------------------------