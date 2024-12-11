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
    Tinlet = RhoIn*VzIn*VzIn*PI*(pow(Ranode,2.0)-pow(Rcathode,2.0));
    for(K= Za+1; K<=numZ; K++)
    Tupstream += 2*PI*RAD*Rho[numR][K]*Vr[numR][K]*Vz[numR][K]
    *deltaZ;
    for(J=1; J<= numR; J++)
    Texit += 2*PI*Rho[J][numZ]*Vz[J][numZ]*Vz[J][numZ]
    *(J-0.5)*deltaR*deltaR;
    Thrust = Texit+Tupstream-Tinlet;
    //The coefficient 0.15 was obtained from Villani.
    MaeckerT = 1.0e-7*(log(Ranode/Rcathode) + 0.15)*Jmax*Jmax;
    Isp = (Thrust/MassFlowRate)/go;
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