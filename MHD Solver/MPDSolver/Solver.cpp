/* CALLING ALL ELEMENTS OF THE SOLVER */
#include "Vars.h"

void Solve()
{

    while (t < TotTime)
    {
        // Calculate the time step required
        // This is updated only at convective time scales
        SetTimeStep();
        //Convective fluxes are calculated only at convective time scales.
        EvalConv();
        //Calculates "EvalDiss()", updates "u[]", and increments time.
        TimeMarch();
        //These variables are computed only at convective time scales.
        EvalParam();
        //Verifies if it is time to write to a file, and then does so.
        WriteFile();
    }//finished computing for "TotTime".


}