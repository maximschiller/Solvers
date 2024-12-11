/* Evaluates the jacobian of the system */
#include "Variables.h"
#include <math.h>

double **rho, **ur, **uz, **bt, **hOld, **h, **BtAvg;
double **Lplus[5][5], **Lminus[5][5], **LplusZ[5][5], **LminusZ[5][5];
double **Lambda[5][5], **LambdaZ[5][5], **AbsLambda[5][5], **AbsLambdaZ[5][5];
double **Aplus[5][5], **Aminus[5][5], **Bplus[5][5], **Bminus[5][5];
double **Rp[5][5], **Lp[5][5], **RpZ[5][5], **LpZ[5][5];
double **R[5][5], **Rinv[5][5], **RZ[5][5], **RinvZ[5][5];
double **M[5][5], **Minv[5][5], **AbsA[5][5], **AbsB[5][5];
double **temp[5][5];

void EvalA(int Zone)
{
    //---------------------------------------------------------
    //Jacobian.CPP; Written by Kameshwaran Sankaran
    //---------------------------------------------------------
    /*************************************************
    Performing calculations in the "r" direction
    **************************************************/
    //"AverageR" calculates the values of the variables at j+(1/2).
    //Since they are stored in the same array, the old values
    //are over-written.
    AverageR(Zone);
    //M= dU/dW, where U=conservative variables, W=primitive variables
    EvalM(Zone);
    EvalMinv(Zone);
    //Calculating the diagonal matrix of eigenvalues, Plus & Minus
    EvalLambdaR(Zone);
    //Calculating the riht eigenvalue matrix and its inverse
    EvalVectorsR(Zone);
    MulMatrix(R, Lplus, temp, Zone);
    //Aplus has non-negative eigenvalues
    MulMatrix(temp, Rinv, Aplus, Zone);
    MulMatrix(R, Lminus, temp, Zone);
    //Aminus has non-positive eigenvalues
    MulMatrix(temp, Rinv, Aminus, Zone);
    SubMatrix(Aplus, Aminus, AbsA, Zone);
    //As mentioned in the beginning, the values at node-points
    //were over-written by values at half-points. "ReCalculate"
    //resets it back to node-points. This also calls the function
    //"EvalVar", thus resetting the speeds to the node-points.
    ReCalculate();
    /**********************************************************
    Performing calculations in the "z" direction
    **********************************************************/
    //"AverageZ" calculates the values of the variables at k+(1/2).
    //Since they are stored in the same array, the old values
    //are over-written.
    AverageZ(Zone);
    //M= dU/dW, where U=conservative variables, W=primitive variables
    EvalM(Zone);
    EvalMinv(Zone);
    //Calculating the diagonal matrix of eigenvalues, Plus & Minus
    EvalLambdaZ(Zone);
    //Calculating the right eigenvalue matrix and its inverse
    EvalVectorsZ(Zone);
    MulMatrix(RZ, LplusZ, temp, Zone);
    //Aplus has non-negative eigenvalues
    MulMatrix(temp, RinvZ, Bplus, Zone);
    MulMatrix(RZ, LminusZ, temp, Zone);
    //Aminus has non-positive eigenvalues
    MulMatrix(temp, RinvZ, Bminus, Zone);
    SubMatrix(Bplus, Bminus, AbsB, Zone);
    //As mentioned in the beginning, the values at node-points
    //were over-written by values at half-points. "ReCalculate"
    //resets it back to node-points. This also calls the function
    //"EvalVar", thus resetting the speeds to the node-points.
    ReCalculate();
}

//-------------------------------------
void AverageR(int Zone)
{
    //Setting temporary variables to the value at node-points
    for (j=Rmin[Zone]+1; j<=Rmax[Zone]; j++)
    {
        for (k=Zmin[Zone]+1; k<=Zmax[Zone]; k++)
        {
            rho[j][k] = Rho[j][k];
            ur[j][k] = Vr[j][k];
            uz[j][k] = Vz[j][k];
            hOld[j][k] = h[j][k];
            bt[j][k] = Bt[j][k];
        }
    }

    //This is Roe averaging. So far, the values of the variables are at
    //node-points.Once this averaging is done, the values are at j+(1/2).
    for (j=Rmin[Zone]+1; j<=Rmax[Zone]; j++)
    {
        for (k=Zmin[Zone]+1; k<=Zmax[Zone]; k++)
        {
            Rho[j][k] = pow(rho[j+1][k]*rho[j][k],0.5);
            Vr[j][k] = ((pow(rho[j+1][k],0.5)*ur[j+1][k])
            +(pow(rho[j][k],0.5)*ur[j][k]))/
            (pow(rho[j+1][k],0.5)+pow(rho[j][k],0.5));
            Vz[j][k] = ((pow(rho[j+1][k],0.5)*uz[j+1][k])
            +(pow(rho[j][k],0.5)*uz[j][k]))/
            (pow(rho[j+1][k],0.5)+pow(rho[j][k],0.5));
            h[j][k] = ((hOld[j+1][k]*pow(rho[j+1][k], 0.5))
            +(hOld[j][k]*pow(rho[j][k], 0.5)))/
            (pow(rho[j+1][k],0.5)+pow(rho[j][k],0.5));
            BtAvg[j][k] = ((bt[j+1][k]*pow(rho[j][k], 0.5))
            +(bt[j][k]*pow(rho[j+1][k], 0.5)))/2;
            p[j][k] = ((**gam-1)/(**gam))*Rho[j][k]*( h[j][k]
            -0.5*((Vr[j][k]*Vr[j][k])+(Vz[j][k]+Vz[j][k]))
            - (BtAvg[j][k]*BtAvg[j][k]/(Mu*rho[j][k])) );
        }//end of "k"
    }//end of "j"

    //The function "EvalVectors" needs values of speeds at j+(1/2).
    //So, call "EvalVar"
    // EvalVar(Zone);
    //To restore the values at node-points, "ReCalculate"
    //function is called by "EvalA".
}
//-------------------------------------
void AverageZ(int Zone)
{
    //Setting temporary variables to the value at node-points
    for (j=Rmin[Zone]+1; j<=Rmax[Zone]; j++)
    {
        for (k=Zmin[Zone]+1; k<=Zmax[Zone]; k++)
        {
            rho[j][k] = Rho[j][k];
            ur[j][k] = Vr[j][k];
            uz[j][k] = Vz[j][k];
            hOld[j][k] = h[j][k];
            bt[j][k] = Bt[j][k];
        }
    }

    //This is Roe averaging. So far, the values of the variables are at
    //node-points.Once this averaging is done, the values are at k+(1/2).
    for (j=Rmin[Zone]+1; j<=Rmax[Zone]; j++)
    {
        for (k=Zmin[Zone]+1; k<=Zmax[Zone]; k++)
        {
            Rho[j][k] = pow(rho[j][k+1]*rho[j][k],0.5);
            Vr[j][k] = ((pow(rho[j][k+1],0.5)*ur[j][k+1])
            +(pow(rho[j][k],0.5)*ur[j][k]))/
            (pow(rho[j][k+1],0.5)+pow(rho[j][k],0.5));
            Vz[j][k] = ((pow(rho[j][k+1],0.5)*uz[j][k+1])
            +(pow(rho[j][k],0.5)*uz[j][k]))/
            (pow(rho[j][k+1],0.5)+pow(rho[j][k],0.5));
            h[j][k] = ((hOld[j][k+1]*pow(rho[j][k+1], 0.5))
            +(hOld[j][k]*pow(rho[j][k], 0.5)))/
            (pow(rho[j][k+1],0.5)+pow(rho[j][k],0.5));
            BtAvg[j][k] = ((bt[j][k+1]*pow(rho[j][k], 0.5))
            +(bt[j][k]*pow(rho[j][k+1], 0.5)))/2;
            p[j][k] = ((**gam-1)/(**gam))*Rho[j][k]*( h[j][k]
            -0.5*((Vr[j][k]*Vr[j][k])+(Vz[j][k]+Vz[j][k]))
            - (BtAvg[j][k]*BtAvg[j][k]/(Mu*rho[j][k])) );
        }//end of "k"
    }//end of "j"
    
    //The function "EvalVectors" needs values of speeds at k+(1/2).
    //So, call "EvalVar"
    // EvalVar(Zone);
    // To restore the values at node-points, "ReCalculate" function
    //is called by "EvalA".
}

//-------------------------------------
void EvalM(int Zone)
{
    //"M" is the Jacobian = dU/dW, i.e., the transformation matrix between
    //primitive and conservation variables.
    for (j=Rmin[Zone]+1; j<=Rmax[Zone]; j++)
    {
        for (k=Zmin[Zone]+1; k<=Zmax[Zone]; k++)
        {
            M[0][0][j][k] = 1.0;
            M[0][1][j][k] = 0.0;
            M[0][2][j][k] = 0.0;
            M[0][3][j][k] = 0.0;
            M[0][4][j][k] = 0.0;
            M[1][0][j][k] = Vr[j][k];
            M[1][1][j][k] = Rho[j][k];
            M[1][2][j][k] = 0.0;
            M[1][3][j][k] = 0.0;
            M[1][4][j][k] = 0.0;
            M[2][0][j][k] = Vz[j][k];
            M[2][1][j][k] = 0.0;
            M[2][2][j][k] = Rho[j][k];
            M[2][3][j][k] = 0.0;
            M[2][4][j][k] = 0.0;
            M[3][0][j][k] = 0.0;
            M[3][1][j][k] = 0.0;
            M[3][2][j][k] = 0.0;
            M[3][3][j][k] = 1;
            M[3][4][j][k] = 0.0;
            M[4][0][j][k] = 0.5*((Vr[j][k]*Vr[j][k])+(Vz[j][k]*Vz[j][k]));
            M[4][1][j][k] = Rho[j][k]*Vr[j][k];
            M[4][2][j][k] = Rho[j][k]*Vz[j][k];
            M[4][3][j][k] = Bt[j][k]/Mu;
            M[4][4][j][k] = 1/(**gam-1);
        }//end of k
    }//end of j
}

//-------------------------------------
void EvalMinv(int Zone)
{
    //"Minv" is the inverse of "M".
    for (j=Rmin[Zone]+1; j<=Rmax[Zone]; j++)
    {
        for (k=Zmin[Zone]+1; k<=Zmax[Zone]; k++)
        {
            Minv[0][0][j][k] = 1.0;
            Minv[0][1][j][k] = 0.0;
            Minv[0][2][j][k] = 0.0;
            Minv[0][3][j][k] = 0.0;
            Minv[0][4][j][k] = 0.0;
            Minv[1][0][j][k] = -Vr[j][k]/Rho[j][k];
            Minv[1][1][j][k] = 1/Rho[j][k];
            Minv[1][2][j][k] = 0.0;
            Minv[1][3][j][k] = 0.0;
            Minv[1][4][j][k] = 0.0;
            Minv[2][0][j][k] = -Vz[j][k]/Rho[j][k];
            Minv[2][1][j][k] = 0.0;
            Minv[2][2][j][k] = 1/Rho[j][k];
            Minv[2][3][j][k] = 0.0;
            Minv[2][4][j][k] = 0.0;
            Minv[3][0][j][k] = 0.0;
            Minv[3][1][j][k] = 0.0;
            Minv[3][2][j][k] = 0.0;
            Minv[3][3][j][k] = 1.0;
            Minv[3][4][j][k] = 0.0;
            Minv[4][0][j][k] = (**gam-1)*0.5*( (Vr[j][k]*Vr[j][k])
            +(Vz[j][k]*Vz[j][k]) );
            Minv[4][1][j][k] = -(**gam-1)*Vr[j][k];
            Minv[4][2][j][k] = -(**gam-1)*Vz[j][k];
            Minv[4][3][j][k] = -(**gam-1)*Bt[j][k]/Mu;
            Minv[4][4][j][k] = **gam-1;
        }//end of k
    }//end of j
}
//-------------------------------------