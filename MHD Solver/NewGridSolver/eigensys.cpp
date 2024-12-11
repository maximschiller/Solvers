/* SCRIPT WHICH EVALUATES THE EIGENSYSTEM OF THE HYPERBOLIC (CONVECTION) EQUATIONS */
#include "Variables.h"
#include <math.h>

double **alfF, **alfS, **betaT, **ZalfF, **ZalfS, **ZbetaT;
double **L1, **L2, **L3, **L4, **L5, **L1m, **L2m, **L3m, **L4m, **L5m, **L1p, **L2p, **L3p, **L4p, **L5p,
        **ZL1, **ZL2, **ZL3, **ZL4, **ZL5, **ZL1m, **ZL2m, **ZL3m, **ZL4m, **ZL5m, **ZL1p, **ZL2p, **ZL3p, **ZL4p, **ZL5p;
// int row, col;

//These are vectors are every point.
double **R1[5], **R2[5], **R3[5], **R4[5], **R5[5], **Left1[5], **Left2[5], **Left3[5], **Left4[5], **Left5[5], 
        **R1z[5], **R2z[5], **R3z[5], **R4z[5], **R5z[5],  **Left1z[5], **Left2z[5], **Left3z[5], **Left4z[5], **Left5z[5];

void EvalLambdaR(int Zone)
{
    for (j=Rmin[Zone]+1; j<=Rmax[Zone]; j++)
    {
        for (k=Zmin[Zone]+1; k<=Zmax[Zone]; k++)
        {
            //These are the eigenvalues of the Jacobian
            L1[j][k] = Vr[j][k]-Cf[j][k];
            L2[j][k] = Vr[j][k];
            L3[j][k] = Vr[j][k];
            L4[j][k] = Vr[j][k];
            L5[j][k] = Vr[j][k]+Cf[j][k];

            //Splitting into positive values...
            L1p[j][k] = (L1[j][k] + fabs(L1[j][k]))/2;
            L2p[j][k] = (L2[j][k] + fabs(L2[j][k]))/2;
            L3p[j][k] = (L3[j][k] + fabs(L3[j][k]))/2;
            L4p[j][k] = (L4[j][k] + fabs(L4[j][k]))/2;
            L5p[j][k] = (L5[j][k] + fabs(L5[j][k]))/2;

            //Splitting into negative values...
            L1m[j][k] = (L1[j][k] - fabs(L1[j][k]))/2;
            L2m[j][k] = (L2[j][k] - fabs(L2[j][k]))/2;
            L3m[j][k] = (L3[j][k] - fabs(L3[j][k]))/2;
            L4m[j][k] = (L4[j][k] - fabs(L4[j][k]))/2;
            L5m[j][k] = (L5[j][k] - fabs(L5[j][k]))/2;

            //Creating the diagonal matrix with positive eigenvalues...
            for (row =0; row<5; row++)
            {
                for (col=0; col<5; col++)
                {
                    Lplus[row][col][j][k] = 0.0;
                }
            }

            Lplus[0][0][j][k] = L1p[j][k];
            Lplus[1][1][j][k] = L2p[j][k];
            Lplus[2][2][j][k] = L3p[j][k];
            Lplus[3][3][j][k] = L4p[j][k];
            Lplus[4][4][j][k] = L5p[j][k];

            //Creating the diagonal matrix with negative eigenvalues...
            for (row =0; row<5; row++)
            {
                for (col =0; col<5; col++)
                {
                    Lminus[row][col][j][k] = 0.0;
                }
            }

            Lminus[0][0][j][k] = L1m[j][k];
            Lminus[1][1][j][k] = L2m[j][k];
            Lminus[2][2][j][k] = L3m[j][k];
            Lminus[3][3][j][k] = L4m[j][k];
            Lminus[4][4][j][k] = L5m[j][k];
        }//end of "k" loop
    }//end of "j" loop
    //The final matrix is created by adding the positive and negative.
    //NOTE: This is not used.
    AddMatrix(Lplus, Lminus, Lambda, Zone);
    SubMatrix(Lplus, Lminus, AbsLambda, Zone);
}

//-------------------------------------
void EvalLambdaZ(int Zone)
{
    for (j=Rmin[Zone]+1; j<=Rmax[Zone]; j++)
    {
        for (k=Zmin[Zone]+1; k<=Zmax[Zone]; k++)
        {
            //These are the eigenvalues of the Jacobian
            ZL1[j][k] = Vz[j][k]-Cf[j][k];
            ZL2[j][k] = Vz[j][k];
            ZL3[j][k] = Vz[j][k];
            ZL4[j][k] = Vz[j][k];
            ZL5[j][k] = Vz[j][k]+Cf[j][k];

            //Splitting into positive values...
            ZL1p[j][k] = (ZL1[j][k] + fabs(ZL1[j][k]))/2;
            ZL2p[j][k] = (ZL2[j][k] + fabs(ZL2[j][k]))/2;
            ZL3p[j][k] = (ZL3[j][k] + fabs(ZL3[j][k]))/2;
            ZL4p[j][k] = (ZL4[j][k] + fabs(ZL4[j][k]))/2;
            ZL5p[j][k] = (ZL5[j][k] + fabs(ZL5[j][k]))/2;

            //Splitting into negative values...
            ZL1m[j][k] = (ZL1[j][k] - fabs(ZL1[j][k]))/2;
            ZL2m[j][k] = (ZL2[j][k] - fabs(ZL2[j][k]))/2;
            ZL3m[j][k] = (ZL3[j][k] - fabs(ZL3[j][k]))/2;
            ZL4m[j][k] = (ZL4[j][k] - fabs(ZL4[j][k]))/2;
            ZL5m[j][k] = (ZL5[j][k] - fabs(ZL5[j][k]))/2;

            //Creating the diagonal matrix with positive eigenvalues...
            for (row =0; row<5; row++)
            {
                for (col=0; col<5; col++)
                {
                    LplusZ[row][col][j][k] = 0.0;
                }
            }

            LplusZ[0][0][j][k] = ZL1p[j][k];
            LplusZ[1][1][j][k] = ZL2p[j][k];
            LplusZ[2][2][j][k] = ZL3p[j][k];
            LplusZ[3][3][j][k] = ZL4p[j][k];
            LplusZ[4][4][j][k] = ZL5p[j][k];

            //Creating the diagonal matrix with negative eigenvalues...
            for (row =0; row<5; row++)
            {
                for (col =0; col<5; col++)
                {
                    LminusZ[row][col][j][k] = 0.0;
                }
            }

            LminusZ[0][0][j][k] = ZL1m[j][k];
            LminusZ[1][1][j][k] = ZL2m[j][k];
            LminusZ[2][2][j][k] = ZL3m[j][k];
            LminusZ[3][3][j][k] = ZL4m[j][k];
            LminusZ[4][4][j][k] = ZL5m[j][k];
        }//end of "k" loop
    }//end of "j" loop
    //The final matrix is created by adding the positive and negative.
    //NOTE: This is not used.
    AddMatrix(LplusZ, LminusZ, LambdaZ, Zone);
    SubMatrix(LplusZ, LminusZ, AbsLambdaZ, Zone);
}

//-------------------------------------
void EvalVectorsR(int Zone)
{
    for (j=Rmin[Zone]+1; j<=Rmax[Zone]; j++)
    {
        for (k=Zmin[Zone]+1; k<=Zmax[Zone]; k++)
        {
            //Definitions:
            alfF[j][k] = sqrt(aSq[j][k]/(Cf[j][k]*Cf[j][k]));
            alfS[j][k] = sqrt(((Cf[j][k]*Cf[j][k]) - aSq[j][k])
            /(Cf[j][k]*Cf[j][k]));
            betaT[j][k] = Bt[j][k]/sqrt(Bt[j][k]*Bt[j][k]);

            //These are the right eigenvectors, R1 through R5.
            R1[0][j][k] = Rho[j][k]*alfF[j][k];
            R1[1][j][k] = -alfF[j][k]*Cf[j][k];
            R1[2][j][k] = 0.0;
            R1[3][j][k] = alfS[j][k]*sqrt(Mu*Rho[j][k]*aSq[j][k])
            *betaT[j][k];
            R1[4][j][k] = alfF[j][k]*gam[j][k]*p[j][k];

            //------------------
            R2[0][j][k] = 0.0;
            R2[1][j][k] = 0.0;
            R2[2][j][k] = betaT[j][k];
            R2[3][j][k] = 0.0;
            R2[4][j][k] = 0.0;

            //------------------
            R3[0][j][k] = 1.0;
            R3[1][j][k] = 0.0;
            R3[2][j][k] = 0.0;
            R3[3][j][k] = 0.0;
            R3[4][j][k] = 0.0;

            //------------------
            R4[0][j][k] = Root2*Rho[j][k]*alfS[j][k];
            R4[1][j][k] = 0.0;
            R4[2][j][k] = 0.0;
            R4[3][j][k] = -Root2*alfF[j][k]*betaT[j][k]
            *sqrt(Mu*Rho[j][k]*aSq[j][k]);
            R4[4][j][k] = Root2*alfS[j][k]*gam[j][k]*p[j][k];

            //------------------
            R5[0][j][k] = Rho[j][k]*alfF[j][k];
            R5[1][j][k] = alfF[j][k]*Cf[j][k];
            R5[2][j][k] = 0.0;
            R5[3][j][k] = alfS[j][k]*betaT[j][k]
            *sqrt(Mu*Rho[j][k]*aSq[j][k]);
            R5[4][j][k] = alfF[j][k]*gam[j][k]*p[j][k];

            //------------------
            //Creating the matrix of right eigenvectors.
            for (row=0; row<5; row++)
            {
            Rp[row][0][j][k] = R1[row][j][k];
            Rp[row][1][j][k] = R2[row][j][k];
            Rp[row][2][j][k] = R3[row][j][k];
            Rp[row][3][j][k] = R4[row][j][k];
            Rp[row][4][j][k] = R5[row][j][k];
            }

            //------------------
            //Creating the inverse of the matrix.
            Left1[0][j][k] = 0.0;
            Left1[1][j][k] = -alfF[j][k]*Cf[j][k]/(2*aSq[j][k]);
            Left1[2][j][k] = 0.0;
            Left1[3][j][k] = alfS[j][k]*betaT[j][k]
            /(2*sqrt(Mu*Rho[j][k]*aSq[j][k]));
            Left1[4][j][k] = alfF[j][k]/(2*Rho[j][k]*aSq[j][k]);

            //------------------
            Left2[0][j][k] = 0.0;
            Left2[1][j][k] = 0.0;
            Left2[2][j][k] = betaT[j][k]/Root2;
            Left2[3][j][k] = 0.0;
            Left2[4][j][k] = 0.0;

            //------------------
            Left3[0][j][k] = 1.0;
            Left3[1][j][k] = 0.0;
            Left3[2][j][k] = 0.0;
            Left3[3][j][k] = 0.0;
            Left3[4][j][k] = -1/aSq[j][k];

            //------------------
            Left4[0][j][k] = 0.0;
            Left4[1][j][k] = 0.0;
            Left4[2][j][k] = 0.0;
            Left4[3][j][k] = -alfF[j][k]*betaT[j][k]
            /sqrt(2*Mu*Rho[j][k]*aSq[j][k]);
            Left4[4][j][k] = alfS[j][k]/(Root2*Rho[j][k]*aSq[j][k]);

            //------------------
            Left5[0][j][k] = 0.0;
            Left5[1][j][k] = alfF[j][k]*Cf[j][k]/(2*aSq[j][k]);
            Left5[2][j][k] = 0.0;
            Left5[3][j][k] = alfS[j][k]*betaT[j][k]
            /(2*sqrt(Mu*Rho[j][k]*aSq[j][k]));
            Left5[4][j][k] = alfF[j][k]/(2*Rho[j][k]*aSq[j][k]);

            //------------------
            //Creating the matrix of left eigenvectors.
            for (col=0; col<5; col++)
            {
                Lp[0][col][j][k] = Left1[col][j][k];
                Lp[1][col][j][k] = Left2[col][j][k];
                Lp[2][col][j][k] = Left3[col][j][k];
                Lp[3][col][j][k] = Left4[col][j][k];
                Lp[4][col][j][k] = Left5[col][j][k];
            }
        }//end of "k" loop
    }//end of "j" loop
    //------------------
    //Obtaining the eigenvectors for the conservative form:
    MulMatrix(M, Rp, R, Zone);
    MulMatrix(Lp, Minv, Rinv, Zone);
}
//-------------------------------------
void EvalVectorsZ(int Zone)
{
    for (j=Rmin[Zone]+1; j<=Rmax[Zone]; j++)
    {
        for (k=Zmin[Zone]+1; k<=Zmax[Zone]; k++)
        {
            //Definitions:
            ZalfF[j][k] = sqrt(aSq[j][k]/(Cf[j][k]*Cf[j][k]));
            ZalfS[j][k] = sqrt(((Cf[j][k]*Cf[j][k]) - aSq[j][k])
            /(Cf[j][k]*Cf[j][k]));
            ZbetaT[j][k] = Bt[j][k]/sqrt(Bt[j][k]*Bt[j][k]);

            //These are the right eigenvectors, R1z through R5z.
            R1z[0][j][k] = Rho[j][k]*ZalfF[j][k];
            R1z[1][j][k] = 0.0;
            R1z[2][j][k] = -Cf[j][k]*ZalfF[j][k];
            R1z[3][j][k] = ZalfS[j][k]*ZbetaT[j][k]
            *sqrt(Mu*Rho[j][k]*aSq[j][k]);
            R1z[4][j][k] = Rho[j][k]*aSq[j][k]*ZalfF[j][k];

            //------------------
            R2z[0][j][k] = 0.0;
            R2z[1][j][k] = -ZbetaT[j][k]/Root2;
            R2z[2][j][k] = 0.0;
            R2z[3][j][k] = 0.0;
            R2z[4][j][k] = 0.0;

            //------------------
            R3z[0][j][k] = 1.0;
            R3z[1][j][k] = 0.0;
            R3z[2][j][k] = 0.0;
            R3z[3][j][k] = 0.0;
            R3z[4][j][k] = 0.0;

            //------------------
            R4z[0][j][k] = Root2*Rho[j][k]*ZalfS[j][k];
            R4z[1][j][k] = 0.0;
            R4z[2][j][k] = 0.0;
            R4z[3][j][k] = -Root2*ZalfF[j][k]*ZbetaT[j][k]
            *sqrt(Mu*Rho[j][k]*aSq[j][k]);
            R4z[4][j][k] = Root2*Rho[j][k]*aSq[j][k]*ZalfS[j][k];

            //------------------
            R5z[0][j][k] = Rho[j][k]*ZalfF[j][k];
            R5z[1][j][k] = 0.0;
            R5z[2][j][k] = Cf[j][k]*ZalfF[j][k];
            R5z[3][j][k] = ZalfS[j][k]*ZbetaT[j][k]
            *sqrt(Mu*Rho[j][k]*aSq[j][k]);
            R5z[4][j][k] = Rho[j][k]*aSq[j][k]*ZalfF[j][k];

            //------------------
            //Creating the matrix of right eigenvectors.
            for (row=0; row<5; row++)
            {
                RpZ[row][0][j][k] = R1z[row][j][k];
                RpZ[row][1][j][k] = R2z[row][j][k];
                RpZ[row][2][j][k] = R3z[row][j][k];
                RpZ[row][3][j][k] = R4z[row][j][k];
                RpZ[row][4][j][k] = R5z[row][j][k];
            }

            //------------------
            //Creating the inverse of the matrix.
            Left1z[0][j][k] = 0.0;
            Left1z[1][j][k] = 0.0;
            Left1z[2][j][k] = -Cf[j][k]*ZalfF[j][k]/(2*aSq[j][k]);
            Left1z[3][j][k] = ZalfS[j][k]*ZbetaT[j][k]
            /(2*sqrt(Mu*Rho[j][k]*aSq[j][k]));
            Left1z[4][j][k] = ZalfF[j][k]/(2*Rho[j][k]*aSq[j][k]);

            //------------------
            Left2z[0][j][k] = 0.0;
            Left2z[1][j][k] = -ZbetaT[j][k]/Root2;
            Left2z[2][j][k] = 0.0;
            Left2z[3][j][k] = 0.0;
            Left2z[4][j][k] = 0.0;

            //------------------
            Left3z[0][j][k] = 1.0;
            Left3z[1][j][k] = 0.0;
            Left3z[2][j][k] = 0.0;
            Left3z[3][j][k] = 0.0;
            Left3z[4][j][k] = -1/aSq[j][k];

            //------------------
            Left4z[0][j][k] = 0.0;
            Left4z[1][j][k] = 0.0;
            Left4z[2][j][k] = 0.0;
            Left4z[3][j][k] = -ZalfF[j][k]*ZbetaT[j][k]
            /(Root2*sqrt(Mu*Rho[j][k]*aSq[j][k]));
            Left4z[4][j][k] = ZalfS[j][k]/(Root2*Rho[j][k]*aSq[j][k]);

            //------------------
            Left5z[0][j][k] = 0.0;
            Left5z[1][j][k] = 0.0;
            Left5z[2][j][k] = Cf[j][k]*ZalfF[j][k]/(2*aSq[j][k]);
            Left5z[3][j][k] = ZalfS[j][k]*ZbetaT[j][k]
            /(2*sqrt(Mu*Rho[j][k]*aSq[j][k]));
            Left5z[4][j][k] = ZalfF[j][k]/(2*Rho[j][k]*aSq[j][k]);

            //------------------
            //Creating the matrix of left eigenvectors.
            for (col=0; col<5; col++)
            {
                LpZ[0][col][j][k] = Left1z[col][j][k];
                LpZ[1][col][j][k] = Left2z[col][j][k];
                LpZ[2][col][j][k] = Left3z[col][j][k];
                LpZ[3][col][j][k] = Left4z[col][j][k];
                LpZ[4][col][j][k] = Left5z[col][j][k];
            }
        }//end of "k" loop
    }//end of "j" loop
    //------------------
    //Obtaining the eigenvectors for the conservative form:
    MulMatrix(M, RpZ, RZ, Zone);
    MulMatrix(LpZ, Minv, RinvZ, Zone);
}
//-------------------------------------