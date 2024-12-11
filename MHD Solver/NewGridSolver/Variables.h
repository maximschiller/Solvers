#ifndef MY_HEADER_FILE_H
#define MY_HEADER_FILE_H

#include <cstdio>

// Basic variables declaration
const double PI = 3.1415926535897932384626433;
const double Root2=1.4142136;
const double go=9.81; //Acceleration due to gravity at Earth, in m/sˆ2
const double Mu=1.256637061436e-6; //Permeability of free space
const double Eps0=8.854187817e-12; //Permittivity of free space
const double gamConst=1.666666666; //Ideal ratio of specific heats for Argon
const double kBoltz=1.380658e-23; //Boltzmann’s constant in Joules/Kelvin
const double Rargon= 208.1333519883; //Gas constant for Argon in J/Kg.K
const double mAr=6.63352599096e-26; //Mass of Argon ion in Kg
const double mEl=9.1093897e-31; //Mass of election in Kg
const double q=1.60217733e-19; //Charage of electron in Coulombs
const double eV=11604.44751705; //1 electron volt in K

// Material Properties
const double ResTungst=5.6e-8; //Resistivity of tungsten in Ohm.m
const double ResCopper=1.7e-8; //Resistivity of copper in Ohm.m
const double ThermTungst=174.0; //Thermal conductivity of tungsten in W/m.K
const double ThermCopper=401.0; //Thermal conductivity of copper in W/m.K

// Counters
extern int i,j; // counters for cell vertexes
extern int I,J; // counters for cell centres
extern int k,K; // thruster-length counters
extern int imax,jmax; // number of cell vertexes
extern int Imax, Jmax; // maximum number of cell-centres
extern int numR, numZ; // number of discretisations in radial and thruster-length directions
extern int Zone; // zone number counter
// Used for linear algebra and eigensystem calcs
extern int row, col;
extern int l;

// Main variables
extern double RAD, Z;
extern double deltaR, deltaZ;
extern double dx, dy;

// Time counters
extern int numSteps, SubStep, MaxStep; // number of steps
extern int t;


// Residuals
extern double eps;


/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/


// Grid Generating variables declaration
// domain.cpp
extern int iter;
extern double **x, **y, **Pdom, **Qdom; // domain.cpp matrices

// grideval.cpp
extern double **xc, **yc; // centroid position arrays
extern double **totVol; // Total volume of each cell
extern double **fAn, **fAe, **fAs, **fAw; // face areas
extern double **nvni, **nvnj, **nvei, **nvej, **nvsi, **nvsj, **nvwi, **nvwj; // normal vectors
extern double **ezetani, **ezetanj, **ezetaei, **ezetaej, **ezetasi, **ezetasj, **ezetawi, **ezetawj; // cell centre to cell centre vector
extern double **eetani, **eetanj, **eetaei, **eetaej, **eetasi, **eetasj, **eetawi, **eetawj; // face area vector
extern double **thetan, **thetae, **thetas, **thetaw;
extern double *Ranode, *Rcathode; // radius of anode/cathode at different z-direction positions
extern int Ra, Rc, Za, Zc; // Grid positioning of the anode and cathode in radial and thruster-length direction
extern int Rmin[4], Rmax[4], Zmin[4], Zmax[4]; //Vertices of the zone.

// io.cpp
extern double mdot, JTot, TotTime;

// memory.cpp
extern double **Rho, **Vr, **Vz, **VSq, **Bt, **BSq, **E; // Conservation Variables
extern double **n, **p, **P, **Er, **Ez, **ErH, **EzH, **EFr, **EFz, **jr, **jz, **jDens, **Jencl, **Potential; // Primitive Variables
extern double **EintEl, **EintH, **ne, **neOld, **ni, **nii, **niii, **nA, **pe, **ph, **Te, **Th, **FrEn, **FzEn,
        **DrEn, **DzEn, **SrEn, **ElecComp, **Ohmic, **Exchange, **EnFluxR, **EnFluxZ, **EnSourceR; // Variables from the two temperature model
extern double **Res, **EIcollFreq, **AnomFreqEl, **TotCollFreq, **ElecGyro, **ElecHall, **CoulombLog, **kTherm,
        **kIon, **IonCondR, **IonCondZ, **ElCondR, **ElCondZ, **ThermCondR, **ThermCondZ; //Dissipation variables
extern double **Zeff, **Qei, **Qeii, **Qeiii, **nuei, **nueii, **nueiii, **CGvosQ, **Ce; //Transport
extern double **aSq, **CmSq, **Cf, **SpecRadR, **SpecRadZ; // Speeds
extern double **Vde, **Vti, **AnomFreqEl; //For anomalous transport
extern double **PR, **gam, **gamOld; //For the equation of state
extern double **KEold, **KE, **Etherm, **EthOld; //For convergence checks
extern double **u[5], **uOld[5], **ConvFlux[5], **Dr[5], **Fr[5], **Hr[5], **SourceR[5], **SourceConv[5], **Lr[5], **Sr[5], **Dz[5], **Fz[5],
                **Hz[5], **Lz[5], **Sz[5], **RdelU[5], **ZdelU[5], **DissFlux[5], **DissipR[5], **SourceDiss[5], **DissipZ[5]; // Vector variables

// initialcon.cpp
extern double RhoIn, VzIn, peIn, TeIn, phIn, ThIn; // declaring the inlet value variables
extern double Jback, Vappl; //Initial guesses.
extern FILE *ConvergeDiff, *ConvergeVal, *Converge ;

// boundarycon.cpp


// recalc.cpp


// solver.cpp


// timestep.cpp
extern double deltaT;


// calcconv.cpp


// timemarch.cpp


// calcdiss.cpp


// calcenergy.cpp
void EvalEnergy();


// eigensys.cpp


// mjacobian.cpp
extern double **rho, **ur, **uz, **bt, **hOld, **h, **BtAvg;
extern double **Lplus[5][5], **Lminus[5][5], **LplusZ[5][5], **LminusZ[5][5];
extern double **Lambda[5][5], **LambdaZ[5][5], **AbsLambda[5][5], **AbsLambdaZ[5][5];
extern double **Aplus[5][5], **Aminus[5][5], **Bplus[5][5], **Bminus[5][5];
extern double **Rp[5][5], **Lp[5][5], **RpZ[5][5], **LpZ[5][5];
extern double **R[5][5], **Rinv[5][5], **RZ[5][5], **RinvZ[5][5];
extern double **M[5][5], **Minv[5][5], **AbsA[5][5], **AbsB[5][5];
extern double **temp[5][5];



// mlinearalgebra.cpp



/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

// Function declarations
// domain.cpp
void DomainInput();
void discretiseGeom();
void testgeom();
void gridgen();

// grideval.cpp
void cellCentres();
void cellGeom();
void gridError();
void dimensions();

// io.cpp
void GetInput();
void outputResults();

// memory.cpp
void MemAllocate();
double **AllocMatrix(int RowBegin,int RowEnd,int ColBegin,int ColEnd);

// initialcon.cpp
void setInitial();

// boundarycon.cpp
void ConvBound();
void DissBound();

// recalc.cpp
void ReCalculate();
void EvalParam();
void EquationOfState();
void TerminalChars();
void EvalGamma();
void Saha();
void Transport();


// solver.cpp
void sysSolver();

// timestep.cpp
void EvalTimeStep();

// calcconv.cpp
void EvalConv();
void EvalF();
void EvalNumDiss();
void EvalLim();
void EvalS();
double EvalSign(double x);


// timemarch.cpp
void TimeMarch();


// calcdiss.cpp
void EvalDiss();

// calcenergy.cpp
void EvalEnergy();


// eigensys.cpp
void EvalLambdaR(int Zone);
void EvalLambdaZ(int Zone);
void EvalVectorsR(int Zone);
void EvalVectorsZ(int Zone);


// mjacobian.cpp
void EvalA(int Zone);
void AverageR(int Zone);
void AverageZ(int Zone);
void EvalM(int Zone);
void EvalMinv(int Zone);


// mlinearalgebra.cpp
void MulMatrix(double **MatA[5][5], double **MatB[5][5], double **MatC[5][5], int Zone);
void AddMatrix(double **MatA[5][5], double **MatB[5][5], double **MatC[5][5], int Zone);
void SubMatrix(double **MatA[5][5], double **MatB[5][5], double **MatC[5][5], int Zone);
void MatTimesVec(double **MatA[5][5], double **VecB[5], double **VecC[5], int Zone);
void AddVector(double **VecA[5], double **VecB[5], double **VecC[5], int Zone);
void SubVector(double **VecA[5], double **VecB[5], double **VecC[5], int Zone);



#endif