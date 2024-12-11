/* RECALCULATES ALL THE REQUIRED VARIABLES */
#include "Variables.h"
#include <math.h>
#include<fstream>
#include<iostream>

// // Variable Declarations

int Count; // counter for convergence checks
int In1, In2, In3, In4, In5, In2A, In2B, In3A, In3B, In4A, In4B;
double K0, K1, K2; //These are coefficients in the equation of state.
double DensD=0.0, KED=0.0, BtD=0.0, EtotD=0.0, EthD=0.0, DensDavg=0.0, KEDavg=0.0, BtDavg=0.0,EtotDavg=0.0, EthDavg=0.0;
double DensMaxR, DensMaxZ, KEmaxR, KEmaxZ, BtMaxR, BtMaxZ, EtotMaxR, EtotMaxZ, EthMaxR, EthMaxZ;
double u0Dmax=0.0, u1Dmax=0.0, u2Dmax=0.0, u3Dmax=0.0, u4Dmax=0.0, u0Davg=0.0, u1Davg=0.0, u2Davg=0.0, u3Davg=0.0, u4Davg=0.0,
       u0DMaxR, u0DMaxZ, u1DMaxR, u1DMaxZ, u2DMaxR, u2DMaxZ, u3DMaxR, u3DMaxZ, u4DMaxR, u4DMaxZ, u0Max=0.0, u1Max=0.0,
       u2Max=0.0, u3Max=0.0, u4Max=0.0, u0Avg=0.0, u1Avg=0.0, u2Avg=0.0, u3Avg=0.0, u4Avg=0.0, u0MaxR, u0MaxZ, u1MaxR,
       u1MaxZ, u2MaxR, u2MaxZ, u3MaxR, u3MaxZ, u4MaxR, u4MaxZ;
double Tinlet=0.0, Tupstream=0.0, Texit=0.0, Thrust=0.0,MaeckerT=0.0, Isp;

void ReCalculate()
{
    //Variables for convergence checks.
    u0Dmax = 0.0;
    u1Dmax = 0.0;
    u2Dmax = 0.0;
    u3Dmax = 0.0;
    u4Dmax = 0.0;
    u0Davg = 0.0;
    u1Davg = 0.0;
    u2Davg = 0.0;
    u3Davg = 0.0;
    u4Davg = 0.0;
    DensD = 0.0;
    KED = 0.0;
    BtD = 0.0;
    EtotD = 0.0;
    EthD = 0.0;
    DensDavg= 0.0;
    KEDavg = 0.0;
    BtDavg = 0.0;
    EtotDavg= 0.0;
    EthDavg = 0.0;
    u0Max = 0.0;
    u1Max = 0.0;
    u2Max = 0.0;
    u3Max = 0.0;
    u4Max = 0.0;
    u0Avg = 0.0;
    u1Avg = 0.0;
    u2Avg = 0.0;
    u3Avg = 0.0;
    u4Avg = 0.0;
    Count = 0;
    
    //The variables are stored at cell centers.
    //So, ’J’ refers to ’j-0.5’ and ’K’ refers to ’k-0.5’.
    for(Zone = 1; Zone <= 3; Zone++)
    {
        for( J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
        {
            for (K=Zmin[Zone]+1; K<=Zmax[Zone]; K++)
            {
                //Storing some old variables
                KEold[J][K] = KE[J][K];
                EthOld[J][K]= Etherm[J][K];
                //Recalculating new variables
                Rho[J][K] = u[0][J][K];
                Vr[J][K] = u[1][J][K]/u[0][J][K];
                Vz[J][K] = u[2][J][K]/u[0][J][K];
                Bt[J][K] = u[3][J][K];
                E[J][K] = u[4][J][K];
                //Other necessary variables
                VSq[J][K] = (Vr[J][K]*Vr[J][K])+(Vz[J][K]*Vz[J][K]);
                BSq[J][K] = Bt[J][K]*Bt[J][K];
                KE[J][K] = 0.5*Rho[J][K]*VSq[J][K];
                Etherm[J][K]= E[J][K]-KE[J][K]-(0.5*BSq[J][K]/Mu);
                // p[J][K] = (gam[J][K]-1)*Etherm[J][K];
                p[J][K] = (gam[J][K]-1)*(E[J][K]-(0.5*Rho[J][K]
                *VSq[J][K])-(0.5*BSq[J][K]/Mu));
                n[J][K] = Rho[J][K]/mAr;
                Jencl[J][K]= 2*PI*(J-0.5)*deltaR*Bt[J][K]/Mu;
                Potential[J][K] = Potential[J-1][K]
                +((0.5*(DissipZ[3][J][K]
                +DissipZ[3][J][K-1]))-
                (Vz[J][K]*Bt[J][K]));
            }
        }
    }

    //Convergence checks
    for(Zone = 1; Zone <= 3; Zone++)
    {
        for( J=Rmin[Zone]+1; J<=Rmax[Zone]-1; J++)
        {
        for (K=Zmin[Zone]+1; K<=Zmax[Zone]-1; K++)
        {
            if(fabs(u[0][J][K]-uOld[0][J][K])/u[0][J][K]>u0Dmax)
            {
                u0Dmax = fabs(u[0][J][K]-uOld[0][J][K])
                /u[0][J][K];
                u0DMaxR = (J-1)*deltaR;
                u0DMaxZ = (K-1)*deltaZ;
            }

            if(fabs(u[1][J][K]-uOld[1][J][K])/u[1][J][K]>u1Dmax)
            {
                u1Dmax = fabs(u[1][J][K]-uOld[1][J][K])
                /u[1][J][K];
                u1DMaxR = (J-1)*deltaR;
                u1DMaxZ = (K-1)*deltaZ;
            }

            if(fabs(u[2][J][K]-uOld[2][J][K])/u[2][J][K]>u2Dmax)
            {
                u2Dmax = fabs(u[2][J][K]-uOld[2][J][K])
                /u[2][J][K];
                u2DMaxR = (J-1)*deltaR;
                u2DMaxZ = (K-1)*deltaZ;
            }

            if(fabs((u[3][J][K]-uOld[3][J][K])/u[3][J][K])>u3Dmax)
            {
                u3Dmax = fabs((u[3][J][K]-uOld[3][J][K])
                /u[3][J][K]);
                u3DMaxR = (J-1)*deltaR;
                u3DMaxZ = (K-1)*deltaZ;
            }

            if(fabs(u[4][J][K]-uOld[4][J][K])/u[4][J][K]>u4Dmax)
            {
                u4Dmax = fabs(u[4][J][K]-uOld[4][J][K])
                /u[4][J][K];
                u4DMaxR = (J-1)*deltaR;
                u4DMaxZ = (K-1)*deltaZ;
            }

            u0Davg += fabs(u[0][J][K]-uOld[0][J][K])/u[0][J][K];
            u1Davg += fabs(u[1][J][K]-uOld[1][J][K])/u[1][J][K];
            u2Davg += fabs(u[2][J][K]-uOld[2][J][K])/u[2][J][K];
            u3Davg += fabs((u[3][J][K]-uOld[3][J][K])/u[3][J][K]);
            u4Davg += fabs(u[4][J][K]-uOld[4][J][K])/u[4][J][K];
            
            if(u[0][J][K] > u0Max)
            {
                u0Max = u[0][J][K];
                u0MaxR = (J-1)*deltaR;
                u0MaxZ = (K-1)*deltaZ;
            }

            if(u[1][J][K] > u1Max)
            {
                u1Max = u[1][J][K];
                u1MaxR = (J-1)*deltaR;
                u1MaxZ = (K-1)*deltaZ;
            }

            if(u[2][J][K] > u2Max)
            {
                u2Max = u[2][J][K];
                u2MaxR = (J-1)*deltaR;
                u2MaxZ = (K-1)*deltaZ;
            }

            if(fabs(u[3][J][K]) > fabs(u3Max))//Bt < 0.
            {
                u3Max = fabs(u[3][J][K]);
                u3MaxR = (J-1)*deltaR;
                u3MaxZ = (K-1)*deltaZ;
            }

            if(u[4][J][K] > u4Max)
            {
                u4Max = u[4][J][K];
                u4MaxR = (J-1)*deltaR;
                u4MaxZ = (K-1)*deltaZ;
            }

            u0Avg += u[0][J][K];
            u1Avg += u[1][J][K];
            u2Avg += u[2][J][K];
            u3Avg += fabs(u[3][J][K]);
            u4Avg += u[4][J][K];

            //Old convergence checks
            if(fabs(Rho[J][K]-uOld[0][J][K])/Rho[J][K] > DensD)
            {
                DensD =fabs(Rho[J][K]-uOld[0][J][K])/Rho[J][K];
                DensMaxR=(J-1)*deltaR;
                DensMaxZ=(K-1)*deltaZ;
            }

            if(fabs(KE[J][K]-KEold[J][K])/KE[J][K] > KED)
            {
                KED =fabs(KE[J][K]-KEold[J][K])/KE[J][K];
                KEmaxR =(J-1)*deltaR;
                KEmaxZ =(K-1)*deltaZ;
            }

            if(fabs((u[3][J][K]-uOld[3][J][K])/u[3][J][K])>BtD)
            {
                BtD =fabs((u[3][J][K]-uOld[3][J][K])
                /u[3][J][K]);
                BtMaxR =(J-1)*deltaR;
                BtMaxZ =(K-1)*deltaZ;
            }

            if(fabs(E[J][K]-uOld[4][J][K])/E[J][K] > EtotD)
            {
                EtotD =fabs(E[J][K]-uOld[4][J][K])/E[J][K];
                EtotMaxR=(J-1)*deltaR;
                EtotMaxZ=(K-1)*deltaZ;
            }

            if(fabs(Etherm[J][K]-EthOld[J][K])/Etherm[J][K]>EthD)
            {
                EthD =fabs(Etherm[J][K]-EthOld[J][K])
                /Etherm[J][K];
                EthMaxR =(J-1)*deltaR;
                EthMaxZ =(K-1)*deltaZ;
            }

            DensDavg+= fabs(Rho[J][K]-uOld[0][J][K])/Rho[J][K];
            KEDavg += fabs(KE[J][K]-KEold[J][K])/KE[J][K];
            BtDavg += fabs((u[3][J][K]-uOld[3][J][K])
            /u[3][J][K]);
            EtotDavg+= fabs(E[J][K]-uOld[4][J][K])/E[J][K];
            EthDavg += fabs(Etherm[J][K]-EthOld[J][K])
            /Etherm[J][K];
            Count++;
            }
        }
    }

    //Convergence checks
    u0Davg /= Count;
    u1Davg /= Count;
    u2Davg /= Count;
    u3Davg /= Count;
    u4Davg /= Count;
    u0Avg /= Count;
    u1Avg /= Count;
    u2Avg /= Count;
    u3Avg /= Count;
    u4Avg /= Count;

    //Old convergence checks
    DensDavg /= Count;
    KEDavg /= Count;
    BtDavg /= Count;
    EtotDavg /= Count;
    EthDavg /= Count;

    if(numSteps%50 == 0)
    {
        fprintf(ConvergeDiff, "%i \t %f \t %f \t %e \t %f \t %f \t %e \t %f \t %f \t %e \t %f \t %f \t %e \t %f \t %f \t %e \t %e \t %e \t %e \t %e \t %e \n",
        numSteps, u0DMaxR, u0DMaxZ, u0Dmax, u1MaxR, u1MaxZ,
        u1Dmax, u2DMaxR, u2DMaxZ, u2Dmax, u3DMaxR, u3DMaxZ,
        u3Dmax, u4DMaxR, u4DMaxZ, u4Dmax, u0Davg, u1Davg,
        u2Davg, u3Davg, u4Davg);

        fprintf(ConvergeVal, "%i \t %f \t %f \t %e \t %f \t %f \t %e \t %f \t %f \t %e \t %f \t %f \t %e \t %f \t %f \t %e \t %e \t %e \t %e \t %e \t %e \n",
        numSteps, u0MaxR, u0MaxZ, u0Max, u1MaxR, u1MaxZ,
        u1Max, u2MaxR, u2MaxZ, u2Max, u3MaxR, u3MaxZ,
        u3Max, u4MaxR, u4MaxZ, u4Max, u0Avg, u1Avg,
        u2Avg, u3Avg, u4Avg);

        fprintf(Converge, "%i \t %f \t %f \t %e \t %f \t %f \t %e \t %f \t %f \t %e \t %f \t %f \t %e \t %f \t %f \t %e \t %e \t %e \t %e \t %e \t %e \n",
        numSteps, DensMaxR, DensMaxZ, DensD, KEmaxR,
        KEmaxZ, KED, BtMaxR, BtMaxZ, BtD, EtotMaxR,
        EtotMaxZ, EtotD, EthMaxR, EthMaxZ, EthD,
        DensDavg, KEDavg, BtDavg, EtotDavg, EthDavg);
    }

    for (J=Rc+1; J<=Ra; J++)
    {
        Jback += Jencl[J][1];
    }

    //This is the total current flowing in the channel.
    Jback /= Ra-Rc;
}
//-------------------------------------
void EvalParam()//This is called every convective time scale.
{
    for(Zone = 1; Zone <= 3; Zone++)
    {
        for (j=J=Rmin[Zone]+2; J<=Rmax[Zone]-1; j++,J++)
        {
            for(k=K=Zmin[Zone]+2; K<=Zmax[Zone]-1; k++,K++)
            {
                jDens[J][K] = sqrt(pow(0.5*(jr[J][k]
                +jr[J][k-1]),2.0)+
                pow(0.5*(jz[j][K]
                +jz[j-1][K]),2.0));
            }
        }

        for (J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
        {
            for (K=Zmin[Zone]+1; K<=Zmax[Zone]; K++)
            {
                P[J][K] = p[J][K]+(0.5*BSq[J][K]/Mu);
                aSq[J][K] = gam[J][K]*p[J][K]/Rho[J][K];
                CmSq[J][K]= aSq[J][K]+(BSq[J][K]/(Mu*Rho[J][K]));
                Cf[J][K] = sqrt(CmSq[J][K]);//Cf = Cfz.

                //Largest eigenvalue
                if (Vr[J][K]<0)
                SpecRadR[J][K] = fabs(Vr[J][K]-Cf[J][K]);
                else
                SpecRadR[J][K] = fabs(Vr[J][K]+Cf[J][K]);
                if (Vz[J][K]<0)
                SpecRadZ[J][K] = fabs(Vz[J][K]-Cf[J][K]);
                else
                SpecRadZ[J][K] = fabs(Vz[J][K]+Cf[J][K]);
                CoulombLog[J][K] = log(sqrt(Eps0*kBoltz*Te[J][K]
                /(q*q*n[J][K]))/
                (q*q/(12*PI*Eps0*kBoltz
                *Te[J][K])));
                EIcollFreq[J][K] = 3.63312216e-06*n[J][K]
                *CoulombLog[J][K]
                *pow(Te[J][K],-1.5);
                if(EIcollFreq[J][K] < 1.0e9)
                EIcollFreq[J][K]= 1.0e9;
                ElecGyro[J][K] = q*sqrt(BSq[J][K])/mEl;
                ElecHall[J][K] = ElecGyro[J][K]
                /EIcollFreq[J][K];
                Vde[J][K] = jDens[J][K]/(q*n[J][K]);
                Vti[J][K] = sqrt(2*kBoltz*Th[J][K]/mAr);

                if(Vde[J][K]/Vti[J][K] > 1.5)
                {
                    AnomFreqEl[J][K]
                    = EIcollFreq[J][K]*(((0.192)
                    +(3.33e-2*ElecHall[J][K])
                    +(0.212*pow(ElecHall[J][K],2.0))
                    +(-8.27e-5*pow(ElecHall[J][K],3.0)))
                    +((Th[J][K]/Te[J][K])*((1.23e-3)
                    +(-1.58e-2*ElecHall[J][K])
                    +(-7.89e-3*pow(ElecHall[J][K],2.0)))));
                }
                else
                {
                    AnomFreqEl[J][K]= 0.0;
                }

                TotCollFreq[J][K]=EIcollFreq[J][K] + AnomFreqEl[J][K];
                Res[J][K] =mEl*TotCollFreq[J][K]/(n[J][K]*q*q);

                //This is so that the near vacuum regions
                //do not set time step constraint.
                if(Res[J][K] > 7.5e-4)
                Res[J][K] = 7.5e-4;
                kTherm[J][K] = 3.20*kBoltz*kBoltz*n[J][K]*Te[J][K]/
                (mEl*EIcollFreq[J][K]);
                kIon[J][K] = 2.84586e-13*pow(Th[J][K],2.5)
                /log(1.239e7*sqrt(pow(Th[J][K],3.0)
                /n[J][K]));
                if(kTherm[J][K] > 20.0)   
                kTherm[J][K] = 20.0;
            }
        }
    }
    
    //end for all zones.
    //At the inlet, the value of SpecRadZ is needed.
    for (J=Rc+1; J<=Ra; J++)
    {
        aSq[J][0] = gam[J][K]*p[J][0]/Rho[J][0];
        CmSq[J][0] = aSq[J][0]+(BSq[J][0]/(Mu*Rho[J][0]));
        Cf[J][0] = sqrt(CmSq[J][0]);
        SpecRadZ[J][0] = Vz[J][0]+Cf[J][0];
    }
}

//-------------------------------------
void EquationOfState()
{
    //Calculates the temperature based on the ratio of pressure to density.
    //The coefficients are obtained by a fit to the data
    //from a table of partition functions.
    for(Zone = 1; Zone <= 3; Zone++)
    {
        for (J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
        {
            for (K=Zmin[Zone]+1; K<=Zmax[Zone]+1; K++)
            {
                PR[J][K] = ph[J][K]/Rho[J][K];

                //The model used above is only good for temperatures above ˜5000 K.
                if(PR[J][K] <= 2994000)
                In1 = 1;
                else
                In1 = 0;
                if(PR[J][K] > 2994000)
                In2B = 1;
                else
                In2B = 0;
                if(PR[J][K] <= 5507000)
                In2A = 1;
                else
                In2A = 0;
                In2 = In2A * In2B;
                if(PR[J][K] > 5507000)
                In3B =1;
                else
                In3B = 0;
                if(PR[J][K] <= 9827000)
                In3A = 1;
                else
                In3A = 0;
                In3 = In3A * In3B;
                if(PR[J][K] > 9827000)
                In4B = 1;
                else
                In4B = 0;
                if(PR[J][K] <= 18290000)
                In4A = 1;
                else
                In4A = 0;
                In4 = In4A * In4B;
                if(PR[J][K] > 18290000)
                In5 = 1;
                else
                In5 = 0;

                K0 = (In2*7935) + (In4*12460) + (In5*14820);
                K1 = (In1*0.00599) + (In2*0.00119) + (In3*0.00317)
                + (In4*0.00094) + (In5*0.000811);
                K2 = (-7.18e-10 * In1) + (-9.79e-11 * In3);

                if(PR[J][K] < 1.0e6)
                {
                    Th[J][K] = ph[J][K]/((Rho[J][K]/mAr)*kBoltz);
                    gam[J][K]= 5.0/3.0;
                }
                else
                {
                    Th[J][K] = ( ((K2*PR[J][K])+K1)*PR[J][K] ) + K0;
                }
            }
        }
    }//end for all zones.
}

//-------------------------------------
void TerminalChars()
{
    Tinlet = RhoIn*VzIn*VzIn*PI*(pow(Ranode[0],2.0)-pow(Rcathode[0],2.0));
    for(K= Za+1; K<=numZ; K++){
            Tupstream += 2*PI*RAD*Rho[numR][K]*Vr[numR][K]*Vz[numR][K]*deltaZ;
    }
    for(J=1; J<= numR; J++){
            Texit += 2*PI*Rho[J][numZ]*Vz[J][numZ]*Vz[J][numZ]*(J-0.5)*deltaR*deltaR;
    }
    Thrust = Texit+Tupstream-Tinlet;

    //The coefficient 0.15 was obtained from Villani.
    MaeckerT = 1.0e-7*(log(Ranode[0]/Rcathode[0]) + 0.15)*JTot*JTot; // ?? check which values of Ranode and Rcathode should be using
    Isp = (Thrust/mdot)/go;
}

//-------------------------------------
void EvalGamma()
{
    for(Zone = 1; Zone <= 3; Zone++)
    {
        for (J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
        {
            for (K=Zmin[Zone]+1; K<=Zmax[Zone]; K++)
            {
                gamOld[J][K] = gam[J][K];
                if(Th[J][K] < 8000.0)
                {
                    gam[J][K] = 5.0/3.0;
                }
                else
                {
                    if(Th[J][K] < 13000.0)
                    {
                        gam[J][K] = 1.11217+(0.529956*exp(-pow
                        ((Th[J][K]-8050.61)/1318.85,2.0)));
                    }
                    else
                    {
                        if(Th[J][K] < 40000.0)
                        {
                            gam[J][K] = 1.1054+(0.0252666*exp(
                            -pow((Th[J][K]-15142.8)/2394.06,2.0)));
                        }
                        //Above 40000K, gamma is almost a const ˜ 1.1.
                        else
                        {
                            gam[J][K] = 1.10;
                        }
                    }
                }

                //Relaxation for gamma.
                gam[J][K] = (0.99999*gamOld[J][K])+(0.00001*gam[J][K]);
                // EintH[J][K] = ph[J][K]/(gam[J][K]-1);
                // E[J][K] = EintEl[J][K]+EintH[J][K]+(0.5*Rho[J][K]*VSq[J][K])+(0.5*BSq[J][K]/Mu);
                // u[4][J][K] = E[J][K];
            }//end of ’k’
        }//end of ’j’
    }//end for all zones.
}
//-------------------------------------
void Saha()
{
    //From the given Te & n, it computes ne, ni, nii, niii and nA.
    for(Zone = 1; Zone <= 3; Zone++)
    {
    for (J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
        {
        for (K=Zmin[Zone]+1; K<=Zmax[Zone]; K++)
        {
            //Normalizing
            double h=exp(1.5*log(Te[J][K]))*2.41500819e21/n[J][K];

            //Local variables
            double K1 = 11. * h * exp(-1.82892573e5/Te[J][K]);
            double K2 = 3.72 * h * exp(-3.20294100e5/Te[J][K]);
            double K3 = 1.6 * h * exp(-4.74638720e5/Te[J][K]);
            double K4 = 2.1 * h * exp(-6.92810064e5/Te[J][K]);
            double K5 = 1.04 * h * exp(-8.72685373e5/Te[J][K]);
            double K6 = 0.48 * h * exp(-1.05836311e6/Te[J][K]);
            double K12 = K1*K2;
            double K123 = K12*K3;
            double K1234 = K123*K4;
            double K12345 = K1234*K5;
            double K123456 = K12345*K6;
            double a[8];
            a[7] = 1.;
            a[6] = K1;
            a[5] = K12 - K1;
            a[4] = K123 - 2.*K12;
            a[3] = K1234 - 3.*K123;
            a[2] = K12345 - 4.*K1234;
            a[1] = K123456 - 5.*K12345;
            a[0] = -6.*K123456;
            double neN=6.; //Normalized initial guess for ne=6*n
            double ne_old;

            do
            {
                double poly=a[7];
                double dp=0.;
                for(int i=6;i>=0;i--)
                {
                    dp=dp*neN+poly;
                    poly=poly*neN+a[i];
                }
                ne_old=neN;
                neN = neN - poly/dp; // Newton-Raphson
            }

            while(fabs(neN-ne_old)/neN > 1.e-6);
            double ne6=neN*neN*neN;
            ne6=ne6*ne6;
            nA[J][K] =n[J][K]*ne6/(K123456+neN*(K12345+neN
            *(K1234+neN*(K123+neN*(K12+neN*(neN+K1))))));
            ni[J][K] =nA[J][K]*K1/neN;
            nii[J][K] =ni[J][K]*K2/neN;
            niii[J][K]=nii[J][K]*K3/neN;
            //Ensuring that there is at least 1 particle of every species.
            nA[J][K] = (nA[J][K]<1.) ? 1. : nA[J][K];
            ni[J][K] = (ni[J][K]<1.) ? 1. : ni[J][K];
            nii[J][K] = (nii[J][K]<1.) ? 1. : nii[J][K];
            niii[J][K] = (niii[J][K]<1.) ? 1. : niii[J][K];
            //Saving the old value before updating ne.
            neOld[J][K] = ne[J][K];
            ne[J][K] = neN*n[J][K]; //De-normalizing.
            if(ne[J][K]/n[J][K] < 0.05)
            ne[J][K]= 0.05*n[J][K];
            //Relaxation for ne.
            ne[J][K] = (0.999*neOld[J][K])+(0.001*ne[J][K]);
            pe[J][K] = ne[J][K]*kBoltz*Te[J][K];
            //Freeze everything else
            // EintEl[J][K]= pe[J][K]/(gamConst-1);
            ph[J][K] = n[J][K]*kBoltz*Th[J][K];
            // EintH[J][K] = ph[J][K]/(gam[J][K]-1);
            // E[J][K] = EintEl[J][K] + EintH[J][K]
            // + (0.5*BSq[J][K]/Mu)+ (0.5*Rho[J][K]*VSq[J][K]);
            // u[4][J][K] = E[J][K];
            ph[J][K] = p[J][K] - pe[J][K];
            // EintH[J][K] = E[J][K] - ((0.5*BSq[J][K]/Mu)
            // + (0.5*Rho[J][K]*VSq[J][K])) - EintEl[J][K];
            // ph[J][K] = (gam[J][K]-1)*EintH[J][K];
            }//end of ’K’
        }//end of ’J’
    }//end for all Zones.
}
//-------------------------------------
void Transport()
{
    //Obtained from Joerg Heiermann.
    Zeff[J][K] = ne[J][K]/n[J][K];
    Ce[J][K] = sqrt(2*kBoltz*Te[J][K]/mEl);
    //Computes collision cross sections of electrons with
    //each of the heavy species
    CGvosQ[J][K]= 2.19320509e-10 * log(1. + 1.53478873e14
    *Te[J][K]*Te[J][K]*Te[J][K]/
    (ne[J][K]*Zeff[J][K]*Zeff[J][K]*(Zeff[J][K]+1.)));
    //Computes the collision frequencies of electrons with
    //each of the heavy species
    //Neutrals to be included later: Qea = 4.0e-20.
    Qei[J][K] = CGvosQ[J][K]/(Te[J][K]*Te[J][K]);
    Qeii[J][K] = CGvosQ[J][K]*2*2/(Te[J][K]*Te[J][K]);
    Qeiii[J][K] = CGvosQ[J][K]*3*3/(Te[J][K]*Te[J][K]);
    nuei[J][K] = ni[J][J]*Qei[J][K]*Ce[J][K];
    nueii[J][K] = nii[J][J]*Qeii[J][K]*Ce[J][K];
    nueiii[J][K]= niii[J][J]*Qeiii[J][K]*Ce[J][K];
    //Computes resistivity
    double Sum=0.;
    Sum += nuei[J][K]+nueii[J][K]+nueiii[J][K];
    Res[J][K] = Sum/2.11355633e-8*ne[J][K];
}
//-------------------------------------