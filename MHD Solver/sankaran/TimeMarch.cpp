//---------------------------------------------------------
//TimeMarch.CPP; Written by Kameshwaran Sankaran
//---------------------------------------------------------
void TimeMarch()
{
    //The variables are stored at cell centers. So, ’J’ refers
    //to ’j-0.5’ and ’K’ refers to ’k-0.5’.
    //For fluxes, ’j’ and ’k’ refer to correct values.
    for(SubStep = 1; SubStep <= MaxStep; SubStep++)
    {
        EvalDiss();//Computes physical dissipation.
        for(i=0; i<5; i++)
        {
            for(Zone = 1; Zone <= 3; Zone++)
            {
                for (j=J=Rmin[Zone]+1; J<=Rmax[Zone]; j++,J++)
                {
                    for(k=K=Zmin[Zone]+1; K<=Zmax[Zone]; k++,K++)
                    {
                        uOld[i][J][K] = u[i][J][K];
                        u[i][J][K]+=deltaT*(((DissipR[i][j][k]
                        -DissipR[i][j-1][k])
                        /deltaR)
                        +SourceDiss[i][j-1][k]
                        +((DissipZ[i][j][k]
                        -DissipZ[i][j][k-1])
                        /deltaZ)
                        - ConvFlux[i][J][K] );
                    }
                }
            }
        }
        //The primary variables are computed from the conservation variables.
        ReCalculate();
        //This accounts for the thermal nonequilibrium between
        //electrons and ions.
        EvalEnergy();
        Saha();
        //This calculates the ion temperature.
        EquationOfState();
        EvalGamma();
        t += deltaT;
        numSteps++;
    }//end of "SubStep" loop.
}
//-------------------------------------
