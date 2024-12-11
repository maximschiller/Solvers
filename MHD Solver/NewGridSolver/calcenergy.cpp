/* CALCULATES ENERGY OF THE SYSTEM */
#include "Variables.h"

void EvalEnergy()
{
    //Equivalent to "EvalF".
    for(Zone = 1; Zone <= 3; Zone++)
    {
        for (J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
        {
            for (K=Zmin[Zone]+1; K<=Zmax[Zone]; K++)
            {
                FrEn[J][K] = Vr[J][K]*(EintEl[J][K]+pe[J][K]);
                SrEn[J][K] = FrEn[J][K];
                FzEn[J][K] = Vz[J][K]*(EintEl[J][K]+pe[J][K]);
            }
        }
    }

    //Equivalent to "EvalNumDiss".
    for(Zone = 1; Zone <= 3; Zone++)
    {
        for(j=J=Rmin[Zone]+1; j<=Rmax[Zone]-1; j++,J++)
        {
            for (K=Zmin[Zone]+1; K<=Zmax[Zone]; K++)
            {
                DrEn[j][K]=0.5*0.5*(SpecRadR[J][K]+SpecRadR[J+1][K])
                *(EintEl[J+1][K]-EintEl[J][K]);
            }
        }

        for(J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
        {
            for (k=K=Zmin[Zone]+1; k<=Zmax[Zone]; k++,K++)
            {
                DzEn[J][k]=0.5*0.5*(SpecRadZ[J][K]+SpecRadZ[J][K+1])
                *(EintEl[J][K+1]-EintEl[J][K]);
            }
        }
    }
    //Equivalent to "EvalConv".
    for(Zone = 1; Zone <= 3; Zone++)
    {
        for(j=J=Rmin[Zone]+1; j<=Rmax[Zone]; j++,J++)
        {
            for (K=Zmin[Zone]+1; K<=Zmax[Zone]; K++)
            {
                EnFluxR[j][K] = (0.5*(FrEn[J+1][K]+FrEn[J][K]))
                - DrEn[j][K];
                EnSourceR[j][K] = (0.5*(SrEn[J+1][K]+SrEn[J][K]))
                /((j+1)*deltaR);
            }
        }

        for(J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
        {
            for (k=K=Zmin[Zone]+1; k<=Zmax[Zone]; k++,K++)
            {
                EnFluxZ[J][k] = (0.5*(FzEn[J][K+1]+FzEn[J][K]))
                - DzEn[J][k];
            }
        }
    }

    //Equivalent to "TimeMarch".
    for(Zone = 1; Zone <= 3; Zone++)
    {
        for (j=J=Rmin[Zone]+1; J<=Rmax[Zone]; j++,J++)
        {
            for(k=K=Zmin[Zone]+1; K<=Zmax[Zone]; k++,K++)
            {
                //First term is the compressional energy,
                //the second term is the Ohmic heating, the third term is the
                //e-i energy exchange, the fourth term is the convective flux,
                //the fifth term is the thermal conduction.
                if(J==Rmax[Zone])
                {
                    pe[J+1][K] = pe[J][K];
                    jz[j][K] = 0.5*Bt[J][K]/Mu;
                }

                if(J==Rmin[Zone]+1)
                {
                    pe[J-1][K] = pe[J][K];
                    jz[j-1][K] = 0.5*Bt[J][K]/Mu;
                }

                if((K==Zmax[Zone])&&(Zone==3))
                pe[J][K+1] = pe[J][K];

                if(K==Zmin[Zone]+1)
                {
                    if(Zone==1)
                    {
                        pe[J][K-1] = peIn;
                        jr[J][k-1] = 0.0;
                    }

                    if((Zone==2)&&(J>=Ra))
                    {
                        jr[J][k-1] = 0.5*Bt[J][K]/Mu;
                        pe[J][K-1] = pe[J][K];
                    }

                    if((Zone==3)&&(J<=Rc))
                    {
                        jr[J][k-1] = 0.5*Bt[J][K]/Mu;
                        pe[J][K-1] = pe[J][K];
                    }
                }

                ElecComp[J][K] = (Vr[J][K]*(0.5*(pe[J+1][K]
                -pe[J-1][K])
                /deltaR))
                +(Vz[J][K]*(0.5*(pe[J][K+1]
                -pe[J][K-1])
                /deltaZ));
                Ohmic[J][K] = (0.5*((-jr[J][k]*Er[J][k])
                +(-jr[J][k-1]*Er[J][k-1])))
                +(0.5*((-jz[j][K]*Ez[j][K])
                +(-jz[j-1][K]*Ez[j-1][K])));

                //Tab allignments are intentionally offset.
                Exchange[J][K] = 3*EIcollFreq[J][K]*n[J][K]*mEl
                *kBoltz*(Te[J][K]-Th[J][K])/mAr;
                EintEl[J][K] += deltaT*(ElecComp[J][K]+Ohmic[J][K]-Exchange[J][K]
                -(((EnFluxR[j][K]-EnFluxR[j-1][K])/deltaR)
                +EnSourceR[j-1][K]+((EnFluxZ[J][k]-
                EnFluxZ[J][k-1])/deltaZ))
                +(((ElCondR[j][K]-ElCondR[j-1][K])/deltaR)
                +((ElCondZ[J][k]-ElCondZ[J][k-1])/deltaZ)));

                // pe[J][K] = (gamConst-1)*EintEl[J][K];
                // ph[J][K] = p[J][K] - pe[J][K];
                Te[J][K] = (gamConst-1)*EintEl[J][K]/(ne[J][K]*kBoltz);
            }
        }
    }
}
//-------------------------------------