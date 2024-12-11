//---------------------------------------------------------
//Solve.CPP; Written by Kameshwaran Sankaran
//---------------------------------------------------------
#include "Thruster-Axi.h"

void Solve()
{
    if(t >= TotTime)
    //This is only to deal with the case TotTime <= 0.0.
    {
        ReCalculate();
        WriteFile();
    }

    while(t < TotTime)
    {
        //This is updated only at convective time scales.
        EvalTimeStep();
        //Convective fluxes are calculated only at convective time scales.
        EvalConv();
        //Calculates "EvalDiss()", updates "u[]", and increments time.
        TimeMarch();
        //These variables are computed only at convective time scales.
        EvalParam();
        //Verifies if it is time to write to a file, and then does so.
        WriteFile();
    }//finished computing for "TotTime".
    
    //Closing the files for verifying convergence.
    fclose(ConvergeDiff);
    fclose(ConvergeVal);
    fclose(Converge);
}
//-------------------------------------