/*ALLOCATING MEMORY TO VARIABLES*/

#include "Vars.h"

void MemAllocate()
{
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
    AnomFreqEl =AllocMatrix(0, numR+1, 0, numZ+1);
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
