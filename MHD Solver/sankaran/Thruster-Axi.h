//---------------------------------------------------------
//Thruster-Axi.h; Written by Kameshwaran Sankaran
//---------------------------------------------------------
#ifndef MY_HEADER_FILE_H
#define MY_HEADER_FILE_H

// #include <iostream>
// #include <iomanip>
// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>
// #include <complex>
#include <cstdio>

/////---------------------------------------------------------/////

// FILE *ConvergeDiff, *ConvergeVal;
// FILE *Converge;

extern FILE *ConvergeDiff;
extern FILE *ConvergeVal;
extern FILE *Converge;

/////---------------------------------------------------------/////

const double PI=3.14159265359;
const double Root2=1.4142136;

//Acceleration due to gravity at Earth, in m/sˆ2
const double go=9.81;
//Permeability of free space
const double Mu=1.256637061436e-6;
//Permittivity of free space
const double Eps0=8.854187817e-12;
//Ideal ratio of specific heats for Argon
const double gamConst=1.666666666;
//Boltzmann’s constant in Joules/Kelvin
const double kBoltz=1.380658e-23;
//Gas constant for Argon in J/Kg.K
const double Rargon= 208.1333519883;
//Mass of Argon ion in Kg
const double mAr=6.63352599096e-26;
//Mass of election in Kg
const double mEl=9.1093897e-31;
//Charage of electron in Coulombs
const double q=1.60217733e-19;
//1 electron volt in K
//Resistivity of tungsten in Ohm.m
const double eV=11604.44751705;
const double ResTungst=5.6e-8;
//Resistivity of copper in Ohm.m
const double ResCopper=1.7e-8;
//Thermal conductivity of tungsten in W/m.K
const double ThermTungst=174.0;
//Thermal conductivity of copper in W/m.K
const double ThermCopper=401.0;

/////---------------------------------------------------------/////

// int Ra, Rc, Za, Zc, numR=0, numZ=0, numT=0,
// row, col, i=0, j=0, k=0, l=0, m=0, J=0, K=0, numSteps=0, SubStep;
// int MaxStep=10;

extern int Ra, Rc, Za, Zc, numR, numZ, numT,
    row, col, i, j, k, l, m, J, K, numSteps, SubStep;
extern int MaxStep;

/////---------------------------------------------------------/////


//These are variables for the equation of state.
// int In1, In2, In3, In4, In5, In2A, In2B, In3A, In3B, In4A, In4B;
// int Zone, Count;
// int ExtraMemory =0;
// double AnodeLength,CathodeLength,Ranode,Rcathode,MassFlowRate,Jmax;
// double RAD, Z, deltaR, deltaZ, TotTime, deltaT=1.0e-8,
// deltaTconv, deltaTtherm, deltaTresist, Lmax;
// double t = 0.0;
// double eps = 1.0e-15;//A small number, used to estimate cut offs.
// double Jback=0.0, Vappl = 10.0;//Initial guesses.

extern int In1, In2, In3, In4, In5, In2A, In2B, In3A, In3B, In4A, In4B;
extern int Zone, Count;
extern int ExtraMemory;
extern double AnodeLength, CathodeLength, Ranode, Rcathode, MassFlowRate, Jmax, RAD, Z;
extern double deltaR, deltaZ, TotTime, deltaT, deltaTconv, deltaTtherm, deltaTresist, Lmax;
extern double t, eps, Jback, Vappl;

/////---------------------------------------------------------/////

//Guess for the temperature of the electrodes, in K.
// double ElecTemp = 3.0e3;
// double K0, K1, K2;//These are coefficients in the equation of state.
// double RhoIn, VzIn, peIn, phIn, TeIn, ThIn;
// double Tinlet=0.0, Tupstream=0.0, Texit=0.0, Thrust=0.0,
// MaeckerT=0.0, Isp;
// int Rmin[4], Rmax[4], Zmin[4], Zmax[4]; //Vertices of the zone.

// MyHeader.h
extern double ElecTemp;
extern double K0, K1, K2;
extern double RhoIn, VzIn, peIn, phIn, TeIn, ThIn;
extern double Tinlet, Tupstream, Texit, Thrust, MaeckerT, Isp;
extern int Rmin[4], Rmax[4], Zmin[4], Zmax[4];

/////---------------------------------------------------------/////


//These are scalars at every point.
extern double **Vr, **Vz, **Cf, **SpecRadR, **SpecRadZ,
**VSq, **aSq, **CmSq;
extern double **Bt, **BSq;
extern double **Er, **Ez, **ErH, **EzH, **EFr, **EFz, **jr, **jz,
**jDens, **Jencl, **Potential;
extern double **n, **Rho, **p, **E, **P;
extern double **Res, **EIcollFreq, **ElecGyro, **ElecHall,
**CoulombLog, **kTherm, **ElCondR, **ElCondZ,
**kIon, **IonCondR, **IonCondZ, **ThermCondR,
**ThermCondZ;
extern double **Vde, **Vti, **AnomFreqEl, **TotCollFreq;
extern double **PR, **gam, **gamOld;

//These are vectors at every point.
extern double **u[5], **ConvFlux[5],
**Hr[5], **Fr[5], **Dr[5], **SourceConv[5], **SourceR[5], **Lr[5], **Sr[5],
**Hz[5], **Fz[5], **Dz[5], **Lz[5], **Sz[5],
**DissFlux[5], **DissipR[5], **SourceDiss[5], **DissipZ[5];

//These are from "EvalD()".
extern double **RdelU[5], **ZdelU[5];

/////---------------------------------------------------------/////

//Variables for the two-temperature model.
// double **EintEl, **EintH, **ne, **neOld, **ni, **nii, **niii, **nA,
// **pe, **ph, **Te, **Th, **Zeff, **Qei, **Qeii, **Qeiii,
// **nuei, **nueii, **nueiii, **CGvosQ, **Ce, **FrEn, **FzEn,
// **DrEn, **DzEn, **SrEn, **ElecComp, **Ohmic, **Exchange,
// **EnFluxR, **EnFluxZ, **EnSourceR;
extern double **EintEl, **EintH, **ne, **neOld, **ni, **nii, **niii, **nA,
**pe, **ph, **Te, **Th, **Zeff, **Qei, **Qeii, **Qeiii,
**nuei, **nueii, **nueiii, **CGvosQ, **Ce, **FrEn, **FzEn,
**DrEn, **DzEn, **SrEn, **ElecComp, **Ohmic, **Exchange,
**EnFluxR, **EnFluxZ, **EnSourceR;

/////---------------------------------------------------------/////


//Convergence checks
// double **uOld[5];
// double u0Dmax=0.0, u1Dmax=0.0, u2Dmax=0.0, u3Dmax=0.0, u4Dmax=0.0,
// u0Davg=0.0, u1Davg=0.0, u2Davg=0.0, u3Davg=0.0, u4Davg=0.0,
// u0DMaxR, u0DMaxZ, u1DMaxR, u1DMaxZ, u2DMaxR, u2DMaxZ,
// u3DMaxR, u3DMaxZ, u4DMaxR, u4DMaxZ, u0Max=0.0, u1Max=0.0,
// u2Max=0.0, u3Max=0.0, u4Max=0.0, u0Avg=0.0, u1Avg=0.0,
// u2Avg=0.0, u3Avg=0.0, u4Avg=0.0, u0MaxR, u0MaxZ, u1MaxR,
// u1MaxZ, u2MaxR, u2MaxZ, u3MaxR, u3MaxZ, u4MaxR, u4MaxZ;
// double **KEold, **KE, **Etherm, **EthOld;
// double DensD=0.0, KED=0.0, BtD=0.0, EtotD=0.0, EthD=0.0,
// DensDavg=0.0, KEDavg=0.0, BtDavg=0.0,EtotDavg=0.0,
// EthDavg=0.0, DensMaxR, DensMaxZ, KEmaxR, KEmaxZ,
// BtMaxR, BtMaxZ, EtotMaxR, EtotMaxZ, EthMaxR, EthMaxZ;

extern double **uOld[5];
extern double u0Dmax, u1Dmax, u2Dmax, u3Dmax, u4Dmax,
            u0Davg, u1Davg, u2Davg, u3Davg, u4Davg,
            u0DMaxR, u0DMaxZ, u1DMaxR, u1DMaxZ, u2DMaxR, u2DMaxZ,
            u3DMaxR, u3DMaxZ, u4DMaxR, u4DMaxZ, u0Max, u1Max,
            u2Max, u3Max, u4Max, u0Avg, u1Avg, u2Avg, u3Avg, u4Avg,
            u0MaxR, u0MaxZ, u1MaxR, u1MaxZ, u2MaxR, u2MaxZ,
            u3MaxR, u3MaxZ, u4MaxR, u4MaxZ;
extern double **KEold, **KE, **Etherm, **EthOld;
extern double DensD, KED, BtD, EtotD, EthD,
            DensDavg, KEDavg, BtDavg, EtotDavg, EthDavg,
            DensMaxR, DensMaxZ, KEmaxR, KEmaxZ,
            BtMaxR, BtMaxZ, EtotMaxR, EtotMaxZ, EthMaxR, EthMaxZ;

/////---------------------------------------------------------/////


//These are required by CharSplit().
extern double **rho, **ur, **uz, **bt, **hOld, **h, **BtAvg;
extern double **alfF, **alfS, **betaT, **ZalfF, **ZalfS, **ZbetaT;
extern double **L1, **L2, **L3, **L4, **L5, **L1m, **L2m, **L3m, **L4m,
**L5m, **L1p, **L2p, **L3p, **L4p, **L5p, **ZL1, **ZL2,
**ZL3, **ZL4, **ZL5, **ZL1m, **ZL2m, **ZL3m, **ZL4m,
**ZL5m, **ZL1p, **ZL2p, **ZL3p, **ZL4p, **ZL5p;

//These are vectors are every point.
extern double **R1[5], **R2[5], **R3[5], **R4[5], **R5[5], **Left1[5],
**Left2[5], **Left3[5], **Left4[5], **Left5[5],
**R1z[5], **R2z[5], **R3z[5], **R4z[5], **R5z[5], **Left1z[5],
**Left2z[5], **Left3z[5], **Left4z[5], **Left5z[5];

//These are matrices at every point.
extern double **Lplus[5][5], **Lminus[5][5], **LplusZ[5][5], **LminusZ[5][5];
extern double **Lambda[5][5], **LambdaZ[5][5], **AbsLambda[5][5],
**AbsLambdaZ[5][5];
extern double **Aplus[5][5], **Aminus[5][5], **Bplus[5][5], **Bminus[5][5];
extern double **Rp[5][5], **Lp[5][5], **RpZ[5][5], **LpZ[5][5];
extern double **R[5][5], **Rinv[5][5], **RZ[5][5], **RinvZ[5][5];
extern double **M[5][5], **Minv[5][5], **AbsA[5][5], **AbsB[5][5];
extern double **temp[5][5];


//FUNCTION DECLARATIONS...
double EvalSign(double x);
void ConvBound();
void DissBound();
void EquationOfState();
void EvalConv();
void EvalDiss();
void EvalEnergy();
void EvalF();
void EvalGamma();
void EvalLim();
void EvalNumDiss();
void EvalParam();
void EvalS();
void EvalTimeStep();
void EvalVar();
void GetInput();
void Grid();
void MemAllocate();
void ReCalculate();
void Saha();
void SetInitial();
void Solve();
void TerminalChars();
void TimeMarch();
void Transport();
void WriteFile();



//Memory functions
double **AllocMatrix(int RowBegin,int RowEnd,int ColBegin,int ColEnd);


void ArtVisc();
void CharSplit();
void AllocMemoryEXTRA();
void EvalLambdaR(int Zone);
void EvalLambdaZ(int Zone);
void EvalVectorsR(int Zone);
void EvalVectorsZ(int Zone);
void EvalA(int Zone);
void AverageR(int Zone);
void AverageZ(int Zone);
void EvalM(int Zone);
void EvalMinv(int Zone);
void MulMatrix(double **MatA[5][5], double **MatB[5][5], double **MatC[5][5], int Zone);
void AddMatrix(double **MatA[5][5], double **MatB[5][5], double **MatC[5][5], int Zone);
void SubMatrix(double **MatA[5][5], double **MatB[5][5], double **MatC[5][5], int Zone);
void MatTimesVec(double **MatA[5][5], double **VecB[5], double **VecC[5], int Zone);
void AddVector(double **VecA[5], double **VecB[5], double **VecC[5], int Zone);
void SubVector(double **VecA[5], double **VecB[5], double **VecC[5], int Zone);


#endif