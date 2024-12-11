/* EVALUATES THE CONVECTIVE VALUES AT EACH TIME STEP
    USES THE HYPERBOLIC GODUNOV SCHEME */
#include "Variables.h"
#include <math.h>


// Currently uses a zoning scheme which is not required since the interior values of the domain are the grids which are being considered
void EvalConv()
{
    //IMPORTANT:
    //Refer to notes about indexing scheme for the fluxes.
    //This is the cell center value of convective fluxes. The
    //cell face values are given by Hr and Hz.
    EvalF();
    //This is the cell face value of numerical dissipation.
    EvalNumDiss();
    //While computing SourceConv, remember:
    //The ’r’ is not from the variable, but from the upper flux
    //surface, and the term is from the lower flux surface.
    //From the implementation in EvalF, SourceR evaluates the
    //term at ’J’ = ’j-0.5’. The value at ’j’ is needed.
    for( i=0; i<5; i++)
    {
        for(Zone = 1; Zone <= 3; Zone++)
        {
            for(j=J=Rmin[Zone]+1; j<=Rmax[Zone]; j++,J++)
            {
                for (K=Zmin[Zone]+1; K<=Zmax[Zone]; K++)
                {
                Hr[i][j][K] = (0.5*(Fr[i][J+1][K]+Fr[i][J][K]))
                - Dr[i][j][K];
                SourceConv[i][j][K] = (0.5*(SourceR[i][J+1][K]
                +SourceR[i][J][K]))/((j+1)*deltaR);
                }
            }

            for(J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
            {
                for (k=K=Zmin[Zone]+1; k<=Zmax[Zone]; k++,K++)
                {
                    Hz[i][J][k] = (0.5*(Fz[i][J][K+1]+Fz[i][J][K]))
                    - Dz[i][J][k];
                }
            }
        }
    }

    ConvBound();//Calculates the convective terms at the boundaries.
    for(i=0; i<5; i++)
    {
        for(Zone = 1; Zone <= 3; Zone++)
        {
            for (j=J=Rmin[Zone]+1; J<=Rmax[Zone]; j++,J++)
            {
                for(k=K=Zmin[Zone]+1; K<=Zmax[Zone]; k++,K++)
                {
                    ConvFlux[i][J][K] = ((Hr[i][j][K]-Hr[i][j-1][K])
                    /deltaR)+SourceConv[i][j-1][K]
                    + ((Hz[i][J][k]-Hz[i][J][k-1])
                    /deltaZ);
                }
            }
        }
    }
}

//-------------------------------------
void EvalF()
{
    //The variables are stored at cell centers. So, ’J’ refers to
    //’j-0.5’ and ’K’ refers to ’k-0.5’.
    for(Zone = 1; Zone <= 3; Zone++)
    {
        for (J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
        {
            for (K=Zmin[Zone]+1; K<=Zmax[Zone]; K++)
            {
                //Fluxes in the "r" direction
                Fr[0][J][K] = Rho[J][K]*Vr[J][K];
                Fr[1][J][K] = (Rho[J][K]*Vr[J][K]*Vr[J][K])
                +p[J][K]+(0.5*Bt[J][K]*Bt[J][K]/Mu);
                Fr[2][J][K] = Rho[J][K]*Vr[J][K]*Vz[J][K];
                Fr[3][J][K] = Bt[J][K]*Vr[J][K];
                Fr[4][J][K] = Vr[J][K]*(E[J][K]+P[J][K]);
                //Source terms due to the 1/r dependence.
                SourceR[0][J][K] = Rho[J][K]*Vr[J][K];
                SourceR[1][J][K] = (Rho[J][K]*Vr[J][K]*Vr[J][K])
                +(Bt[J][K]*Bt[J][K]/Mu);
                SourceR[2][J][K] = Rho[J][K]*Vr[J][K]*Vz[J][K];
                SourceR[3][J][K] = 0.0;
                SourceR[4][J][K] = Vr[J][K]*(E[J][K]+P[J][K]);
                //Fluxes in the "z" direction
                Fz[0][J][K] = Rho[J][K]*Vz[J][K];
                Fz[1][J][K] = Rho[J][K]*Vr[J][K]*Vz[J][K];
                Fz[2][J][K] = (Rho[J][K]*Vz[J][K]*Vz[J][K])
                +p[J][K]+(0.5*Bt[J][K]*Bt[J][K]/Mu);
                Fz[3][J][K] = Bt[J][K]*Vz[J][K];
                Fz[4][J][K] = Vz[J][K]*(E[J][K]+P[J][K]);
            }
        }
    }
}

//-------------------------------------
void EvalNumDiss()
{
    //Numerical dissipation.
    EvalLim();
    //The numerical dissipation through boundaries should be zero.
    //Thus, the values at boundaries are overwritten
    //by 0 in BoundaryConditions.

    // ExtraMemory++;//First time it is called, it is =1. Then it is >1.
    // if(ExtraMemory ==1)
    //     AllocMemoryEXTRA();

    for(i=0; i<5;i++)
    {
        for(Zone=1; Zone<=3; Zone++)
        {
                for (j=Rmin[Zone]+1; j<=Rmax[Zone]-1; j++)
            {
                for (k=Zmin[Zone]+1; k<=Zmax[Zone]; k++)
                {
                    RdelU[i][j][K] = u[i][J+1][K] - u[i][J][K];
                }
            }

            for (j=Rmin[Zone]+1; j<=Rmax[Zone]; j++)
            {
                for (k=Zmin[Zone]+1; k<=Zmax[Zone]; k++)
                {
                    ZdelU[i][J][k] = u[i][J][K+1] - u[i][J][K];
                }
            }
        }
    }//end of ’i’ loop.

    //CHAMBER
    EvalA(1);
    MatTimesVec(AbsA, RdelU, Dr, 1);
    MatTimesVec(AbsB, ZdelU, Dz, 1);
    //PLUME
    EvalA(2);
    MatTimesVec(AbsA, RdelU, Dr, 2);
    MatTimesVec(AbsB, ZdelU, Dz, 2);

    /****************************
    for(i=0; i<5;i++)
    {
    for(Zone = 1; Zone <= 3; Zone++)
    {
    for(j=J=Rmin[Zone]+1; j<=Rmax[Zone]-1; j++,J++)
    {
    for (K=Zmin[Zone]+1; K<=Zmax[Zone]; K++)
    {
    RdelU[i][j][K] = u[i][J+1][K] - u[i][J][K];
    Dr[i][j][K]=0.5*0.5*(SpecRadR[J][K]
    +SpecRadR[J+1][K])
    *((RdelU[i][j][K])-Lr[i][j][K]);
    }
    }
    for(J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
    {
    for (k=K=Zmin[Zone]+1; k<=Zmax[Zone]; k++,K++)
    {
    ZdelU[i][J][k]= u[i][J][K+1] - u[i][J][K];
    Dz[i][J][k]=0.5*0.5*(SpecRadZ[J][K]
    +SpecRadZ[J][K+1])
    *((ZdelU[i][J][k])-Lz[i][J][k]);
    }
    }
    }
    }
    ****************/
}

//-------------------------------------
void EvalLim()
{
    //This stupid compiler doesn’t have min(), max() functions.
    EvalS();

    //Alpha-bee
    double a, b, aN, bN, d, q;
    q=0.3;

    for (i=0; i<5; i++)
    {
        for(Zone = 1; Zone <= 3; Zone++)
        {
            for (j=Rmin[Zone]+1; j<=Rmax[Zone]-2; j++)
            {
                for(k=Zmin[Zone]+1; k<=Zmax[Zone]-2; k++)
                {
                    //"r" direction
                    a = q*fabs(u[i][j+2][k] - u[i][j+1][k]);
                    b = fabs(u[i][j][k] - u[i][j-1][k]);
                    if (a<b)
                    aN = a;
                    else
                    aN = b;
                    a = fabs(u[i][j+2][k] - u[i][j+1][k]);
                    b = q*fabs(u[i][j][k] - u[i][j-1][k]);
                    if (a<b)
                    bN = a;
                    else
                    bN = b;
                    if (aN > bN)
                    d = aN;
                    else
                    d = bN;
                    Lr[i][j][k] = Sr[i][j][k]*d;

                    //"z" direction
                    a = q*fabs(u[i][j][k+2] - u[i][j][k+1]);
                    b = fabs(u[i][j][k] - u[i][j][k-1]);
                    if (a<b)
                    aN = a;
                    else
                    aN = b;
                    a = fabs(u[i][j][k+2] - u[i][j][k+1]);
                    b = q*fabs(u[i][j][k] - u[i][j][k-1]);
                    if (a<b)
                    bN = a;
                    else
                    bN = b;
                    if (aN > bN)
                    d = aN;
                    else
                    d = bN;
                    Lz[i][j][k] = Sz[i][j][k]*d;
                }
            }
        }//end for all zones.
    }//end for all "i".
}

//-------------------------------------
void EvalS()
{
    double t1, t2;
    for( i=0; i<5; i++)
    {
        for(Zone = 1; Zone <= 3; Zone++)
        {
            for( j=Rmin[Zone]+1; j<=Rmax[Zone]-2; j++)
            {
                for(k=Zmin[Zone]+1; k<=Zmax[Zone]-2; k++)
                {
                    //"r" direction
                    t1 = EvalSign(u[i][j+2][k] - u[i][j+1][k]);
                    t2 = EvalSign(u[i][j][k] - u[i][j-1][k]);
                    Sr[i][j][k] = 0.5*(t1+t2);
                    //"z" direction
                    t1 = EvalSign(u[i][j][k+2] - u[i][j][k+1]);
                    t2 = EvalSign(u[i][j][k] - u[i][j][k-1]);
                    Sz[i][j][k] = 0.5*(t1+t2);
                }
            }
        }
    }
}

//-------------------------------------
double EvalSign(double x)
{
    double sign;
    if (x == 0)
    sign = 0.0;
    else
    sign = x/fabs(x);
    return sign;
}
//-------------------------------------