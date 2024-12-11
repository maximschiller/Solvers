/* INITIALISES MEMORY FOR MATRICES */
#include "Variables.h"
#include <iostream>

// Domain and Grid Generation variables
double **x, **y, **Pdom, **Qdom; // domain.cpp matrices
double **xc, **yc; // centroid position arrays
double **totVol; // Total volume of each cell
double **fAn, **fAe, **fAs, **fAw; // face areas
double **nvni, **nvnj, **nvei, **nvej, **nvsi, **nvsj, **nvwi, **nvwj; // normal vectors
double **ezetani, **ezetanj, **ezetaei, **ezetaej, **ezetasi, **ezetasj, **ezetawi, **ezetawj; // cell centre to cell centre vector
double **eetani, **eetanj, **eetaei, **eetaej, **eetasi, **eetasj, **eetawi, **eetawj; // face area vector
double **thetan, **thetae, **thetas, **thetaw;
double *Ranode, *Rcathode; // radius of anode/cathode at different z-direction positions

// Solver computation variables
double **Rho, **Vr, **Vz, **VSq, **Bt, **BSq, **E; // Conservation Variables
double **n, **p, **P, **Er, **Ez, **ErH, **EzH, **EFr, **EFz, **jr, **jz, **jDens, **Jencl, **Potential; // Primitive Variables
double **EintEl, **EintH, **ne, **neOld, **ni, **nii, **niii, **nA, **pe, **ph, **Te, **Th, **FrEn, **FzEn,
        **DrEn, **DzEn, **SrEn, **ElecComp, **Ohmic, **Exchange, **EnFluxR, **EnFluxZ, **EnSourceR; // Variables from the two temperature model
double **Res, **EIcollFreq, **TotCollFreq, **ElecGyro, **ElecHall, **CoulombLog, **kTherm,
        **kIon, **IonCondR, **IonCondZ, **ElCondR, **ElCondZ, **ThermCondR, **ThermCondZ; //Dissipation variables
double **Zeff, **Qei, **Qeii, **Qeiii, **nuei, **nueii, **nueiii, **CGvosQ, **Ce; //Transport
double **aSq, **CmSq, **Cf, **SpecRadR, **SpecRadZ; // Speeds
double **Vde, **Vti, **AnomFreqEl; //For anomalous transport
double **PR, **gam, **gamOld; //For the equation of state
double **KEold, **KE, **Etherm, **EthOld; //For convergence checks
double **u[5], **uOld[5], **ConvFlux[5], **Dr[5], **Fr[5], **Hr[5], **SourceR[5], **SourceConv[5], **Lr[5], **Sr[5], **Dz[5], **Fz[5],
        **Hz[5], **Lz[5], **Sz[5], **RdelU[5], **ZdelU[5], **DissFlux[5], **DissipR[5], **SourceDiss[5], **DissipZ[5]; // Vector variables


void MemAllocate(){
 // Allocating dimension to domain and grid variables
    x =AllocMatrix(0, numR+1, 0, numZ+1);
    y =AllocMatrix(0, numR+1, 0, numZ+1);
    Pdom =AllocMatrix(0, numR+1, 0, numZ+1);
    Qdom =AllocMatrix(0, numR+1, 0, numZ+1);

    // Centre position arrays
    xc =AllocMatrix(0, numR+1, 0, numZ+1);
    yc =AllocMatrix(0, numR+1, 0, numZ+1);
    totVol =AllocMatrix(0, numR+1, 0, numZ+1);
    
    // Non-orthogonal face areas and vectors
    fAn =AllocMatrix(0, numR+1, 0, numZ+1);
    fAe =AllocMatrix(0, numR+1, 0, numZ+1);
    fAs =AllocMatrix(0, numR+1, 0, numZ+1);
    fAw =AllocMatrix(0, numR+1, 0, numZ+1);
    nvni =AllocMatrix(0, numR+1, 0, numZ+1);
    nvnj =AllocMatrix(0, numR+1, 0, numZ+1);
    nvei =AllocMatrix(0, numR+1, 0, numZ+1);
    nvej =AllocMatrix(0, numR+1, 0, numZ+1);
    nvsi =AllocMatrix(0, numR+1, 0, numZ+1);
    nvsj =AllocMatrix(0, numR+1, 0, numZ+1);
    nvwi =AllocMatrix(0, numR+1, 0, numZ+1);
    nvwj =AllocMatrix(0, numR+1, 0, numZ+1);
    ezetani =AllocMatrix(0, numR+1, 0, numZ+1);
    ezetanj =AllocMatrix(0, numR+1, 0, numZ+1);
    ezetaei =AllocMatrix(0, numR+1, 0, numZ+1);
    ezetaej =AllocMatrix(0, numR+1, 0, numZ+1);
    ezetasi =AllocMatrix(0, numR+1, 0, numZ+1);
    ezetasj =AllocMatrix(0, numR+1, 0, numZ+1);
    ezetawi =AllocMatrix(0, numR+1, 0, numZ+1);
    ezetawj =AllocMatrix(0, numR+1, 0, numZ+1);
    eetani =AllocMatrix(0, numR+1, 0, numZ+1);
    eetanj =AllocMatrix(0, numR+1, 0, numZ+1);
    eetaei =AllocMatrix(0, numR+1, 0, numZ+1);
    eetaej =AllocMatrix(0, numR+1, 0, numZ+1);
    eetasi =AllocMatrix(0, numR+1, 0, numZ+1);
    eetasj =AllocMatrix(0, numR+1, 0, numZ+1);
    eetawi =AllocMatrix(0, numR+1, 0, numZ+1);
    eetawj =AllocMatrix(0, numR+1, 0, numZ+1);
    thetan =AllocMatrix(0, numR+1, 0, numZ+1);
    thetae =AllocMatrix(0, numR+1, 0, numZ+1);
    thetas =AllocMatrix(0, numR+1, 0, numZ+1);
    thetaw =AllocMatrix(0, numR+1, 0, numZ+1);

    Ranode = new double[numZ+1];
    Rcathode = new double[numZ+1];

 //Allocating the scalars.
    //Conservation variables
    Rho =AllocMatrix(0, numR+1, 0, numZ+1);
    Vr =AllocMatrix(0, numR+1, 0, numZ+1);
    Vz =AllocMatrix(0, numR+1, 0, numZ+1);
    VSq =AllocMatrix(0, numR+1, 0, numZ+1);
    Bt =AllocMatrix(0, numR+1, 0, numZ+1);
    BSq =AllocMatrix(0, numR+1, 0, numZ+1);
    E =AllocMatrix(0, numR+1, 0, numZ+1);

    //Primitive variables
    n =AllocMatrix(0, numR+1, 0, numZ+1);
    p =AllocMatrix(0, numR+1, 0, numZ+1);
    P =AllocMatrix(0, numR+1, 0, numZ+1);
    Er =AllocMatrix(0, numR+1, 0, numZ+1);
    Ez =AllocMatrix(0, numR+1, 0, numZ+1);
    ErH =AllocMatrix(0, numR+1, 0, numZ+1);
    EzH =AllocMatrix(0, numR+1, 0, numZ+1);
    EFr =AllocMatrix(0, numR+1, 0, numZ+1);
    EFz =AllocMatrix(0, numR+1, 0, numZ+1);
    jr =AllocMatrix(0, numR+1, 0, numZ+1);
    jz =AllocMatrix(0, numR+1, 0, numZ+1);
    jDens =AllocMatrix(0, numR+1, 0, numZ+1);
    Jencl =AllocMatrix(0, numR+1, 0, numZ+1);
    Potential=AllocMatrix(0, numR+1, 0, numZ+1);

    //Variables from the two-temperature model.
    EintEl =AllocMatrix(0, numR+1, 0, numZ+1);
    EintH =AllocMatrix(0, numR+1, 0, numZ+1);
    ne =AllocMatrix(0, numR+1, 0, numZ+1);
    neOld =AllocMatrix(0, numR+1, 0, numZ+1);
    ni =AllocMatrix(0, numR+1, 0, numZ+1);
    nii =AllocMatrix(0, numR+1, 0, numZ+1);
    niii=AllocMatrix(0, numR+1, 0, numZ+1);
    nA =AllocMatrix(0, numR+1, 0, numZ+1);
    pe =AllocMatrix(0, numR+1, 0, numZ+1);
    ph =AllocMatrix(0, numR+1, 0, numZ+1);
    Te =AllocMatrix(0, numR+1, 0, numZ+1);
    Th =AllocMatrix(0, numR+1, 0, numZ+1);
    FrEn =AllocMatrix(0, numR+1, 0, numZ+1);
    FzEn =AllocMatrix(0, numR+1, 0, numZ+1);
    DrEn =AllocMatrix(0, numR+1, 0, numZ+1);
    DzEn =AllocMatrix(0, numR+1, 0, numZ+1);
    SrEn =AllocMatrix(0, numR+1, 0, numZ+1);
    ElecComp=AllocMatrix(0, numR+1, 0, numZ+1);
    Ohmic =AllocMatrix(0, numR+1, 0, numZ+1);
    Exchange=AllocMatrix(0, numR+1, 0, numZ+1);
    EnFluxR =AllocMatrix(0, numR+1, 0, numZ+1);
    EnFluxZ =AllocMatrix(0, numR+1, 0, numZ+1);
    EnSourceR =AllocMatrix(0, numR+1, 0, numZ+1);

    //Dissipation variables
    Res =AllocMatrix(0, numR+1, 0, numZ+1);
    EIcollFreq =AllocMatrix(0, numR+1, 0, numZ+1);
    // AnomFreqEl =AllocMatrix(0, numR+1, 0, numZ+1);
    TotCollFreq =AllocMatrix(0, numR+1, 0, numZ+1);
    ElecGyro =AllocMatrix(0, numR+1, 0, numZ+1);
    ElecHall =AllocMatrix(0, numR+1, 0, numZ+1);
    CoulombLog =AllocMatrix(0, numR+1, 0, numZ+1);
    kTherm =AllocMatrix(0, numR+1, 0, numZ+1);
    kIon =AllocMatrix(0, numR+1, 0, numZ+1);
    IonCondR =AllocMatrix(0, numR+1, 0, numZ+1);
    IonCondZ =AllocMatrix(0, numR+1, 0, numZ+1);
    ElCondR =AllocMatrix(0, numR+1, 0, numZ+1);
    ElCondZ =AllocMatrix(0, numR+1, 0, numZ+1);
    ThermCondR =AllocMatrix(0, numR+1, 0, numZ+1);
    ThermCondZ =AllocMatrix(0, numR+1, 0, numZ+1);

    //Transport
    Zeff = AllocMatrix(0, numR+1, 0, numZ+1);
    Qei = AllocMatrix(0, numR+1, 0, numZ+1);
    Qeii = AllocMatrix(0, numR+1, 0, numZ+1);
    Qeiii = AllocMatrix(0, numR+1, 0, numZ+1);
    nuei = AllocMatrix(0, numR+1, 0, numZ+1);
    nueii = AllocMatrix(0, numR+1, 0, numZ+1);
    nueiii = AllocMatrix(0, numR+1, 0, numZ+1);
    CGvosQ = AllocMatrix(0, numR+1, 0, numZ+1);
    Ce = AllocMatrix(0, numR+1, 0, numZ+1);

    //Speeds
    aSq =AllocMatrix(0, numR+1, 0, numZ+1);
    CmSq=AllocMatrix(0, numR+1, 0, numZ+1);
    Cf =AllocMatrix(0, numR+1, 0, numZ+1);
    SpecRadR =AllocMatrix(0, numR+1, 0, numZ+1);
    SpecRadZ =AllocMatrix(0, numR+1, 0, numZ+1);

    //For anomalous transport
    Vde =AllocMatrix(0, numR+1, 0, numZ+1);
    Vti =AllocMatrix(0, numR+1, 0, numZ+1);
    AnomFreqEl =AllocMatrix(0, numR+1, 0, numZ+1);

    //For the equation of state
    PR =AllocMatrix(0, numR+1, 0, numZ+1);
    gam =AllocMatrix(0, numR+1, 0, numZ+1);
    gamOld =AllocMatrix(0, numR+1, 0, numZ+1);

    //For convergence checks
    KEold = AllocMatrix(0, numR+1, 0, numZ+1);
    KE = AllocMatrix(0, numR+1, 0, numZ+1);
    Etherm = AllocMatrix(0, numR+1, 0, numZ+1);
    EthOld = AllocMatrix(0, numR+1, 0, numZ+1);

    //Allocating the vectors.
    for (i=0; i<5; i++)
    {
        u[i] =AllocMatrix(0, numR+1, 0, numZ+1);
        uOld[i] =AllocMatrix(0, numR+1, 0, numZ+1);
        ConvFlux[i] =AllocMatrix(0, numR+1, 0, numZ+1);
        Dr[i] =AllocMatrix(0, numR+1, 0, numZ+1);
        Fr[i] =AllocMatrix(0, numR+1, 0, numZ+1);
        Hr[i] =AllocMatrix(0, numR+1, 0, numZ+1);
        SourceR[i] =AllocMatrix(0, numR+1, 0, numZ+1);
        SourceConv[i] =AllocMatrix(0, numR+1, 0, numZ+1);
        Lr[i] =AllocMatrix(0, numR+1, 0, numZ+1);
        Sr[i] =AllocMatrix(0, numR+1, 0, numZ+1);
        Dz[i] =AllocMatrix(0, numR+1, 0, numZ+1);
        Fz[i] =AllocMatrix(0, numR+1, 0, numZ+1);
        Hz[i] =AllocMatrix(0, numR+1, 0, numZ+1);
        Lz[i] =AllocMatrix(0, numR+1, 0, numZ+1);
        Sz[i] =AllocMatrix(0, numR+1, 0, numZ+1);
        RdelU[i] = AllocMatrix(0, numR+1, 0, numZ+1);
        ZdelU[i] = AllocMatrix(0, numR+1, 0, numZ+1);
        DissFlux[i] =AllocMatrix(0, numR+1, 0, numZ+1);
        DissipR[i] =AllocMatrix(0, numR+1, 0, numZ+1);
        SourceDiss[i] =AllocMatrix(0, numR+1, 0, numZ+1);
        DissipZ[i] =AllocMatrix(0, numR+1, 0, numZ+1);
    }        
}

//-------------------------------------
double **AllocMatrix(int RowBegin,int RowEnd,int ColBegin,int ColEnd)
{
    int i;
    int k,l; // Counters for matrix building
    double **Mat;

    //Creating the blank matrix
    // Mat= new double*[ColEnd-ColBegin+1];
    // Mat -= ColBegin;
    Mat= new double*[RowEnd-RowBegin+1];
    Mat -= RowBegin;
    for(i=RowBegin;i<=RowEnd;i++)
    {
        // Mat[i]= new double[RowEnd-RowBegin+1];
        // Mat[i] -= RowBegin;
        Mat[i]= new double[ColEnd-ColBegin+1];
        Mat[i] -= ColBegin;
    }

    //Initializing
    for(k=RowBegin;k<=RowEnd;k++)
    {
    for(l=ColBegin;l<=ColEnd;l++)
        {
        Mat[k][l]=0.0;
        }
    }
    return Mat;
}
//-------------------------------------