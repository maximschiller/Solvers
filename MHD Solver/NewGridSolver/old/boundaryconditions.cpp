/* SCRIPT WHICH INITIALISES THE BOUNDARY CONDITIONS */
#include "Variables.h"
#include <math.h>

// Variable Definitions
double ElecTemp = 3.0e3; //Guess for the temperature of the electrodes, in K.
double eps = 1.0e-15; //A small number, used to estimate cut offs.

// Currently uses an orthogonal grid - need to rewrite this so that it follows the boundary curves of the actual domain
void ConvBound()
{
    //These are the boundary conditions for convective fluxes
    //Backplate conditions. These variables are defined at
    //an imaginary point K=0. Fluid enters the domain from
    //a porous backplate at sonic conditions.
    for(J=Rc[0]+1; J<=Ra[0]; J++)
    {
        //There is no real need to compute the variables here,
        //but it makes the code for specifying the fluxes neater.
        Te[J][0] = 15000.0;
        Th[J][0] = 15000.0;
        Vz[J][0] = sqrt(gamConst*kBoltz*(Te[J][0]+Th[J][0])/mAr);
        VSq[J][0]= Vz[J][0]*Vz[J][0];
        Rho[J][0]= MassFlowRate/(Vz[J][0]*PI*(pow(Ranode,2.0)
        -pow(Rcathode,2.0)));
        n[J][0] = Rho[J][0]/mAr;
        ne[J][0] = 0.999*n[J][0];
        pe[J][0] = ne[J][0]*kBoltz*Te[J][0];
        ph[J][0] = n[J][0]*kBoltz*Th[J][0];

        //For the given Th & pe, this the value of gamma is 1.13.
        gam[J][0]= 1.13;
        p[J][0] = pe[J][0]+ph[J][0];

        //The field diffuses instantaneously inside an insulator.
        Bt[J][0] = Bt[J][1];
        BSq[J][0]= Bt[J][0]*Bt[J][0];

        //Now, the boundary condition for the fluxes. The point k=0
        //corresponds to the physical backplate boundary.
        Hz[0][J][0] = Rho[J][0]*Vz[J][0];
        Hz[1][J][0] = 0.0;
        Hz[2][J][0] = (Rho[J][0]*VSq[J][0])+p[J][0]
        +(0.5*BSq[J][0]/Mu);

        //The inductive drop is also included in the dissipative part,
        //becuase it is evaluated more frequently.
        Hz[3][J][0] = 0.0;
        E[J][0] = (p[J][0]/(gam[J][0]-1))
        +(0.5*Rho[J][0]*VSq[J][0])
        +(0.5*BSq[J][0]/Mu);
        Hz[4][J][0] = Vz[J][0]*(E[J][0]+p[J][0]
        +(0.5*BSq[J][0]/Mu));
    }

    //Anode inner.
    for(K=1; K<=Za; K++)
    {
        Hr[0][Ra][K] = 0.0;
        // Hr[1][Ra][K] = P[Ra][K];
        Hr[1][Ra][K] = Hr[1][Ra-1][K]
        -(deltaR*SourceConv[1][Ra-1][K]);
        Hr[2][Ra][K] = 0.0;
        Hr[3][Ra][K] = 0.0;
        Hr[4][Ra][K] = 0.0;
    }

    //Anode tip.
    for(J=Ra+1; J<=numR; J++)
    {
        Hz[0][J][Za] = 0.0;
        Hz[1][J][Za] = 0.0;
        Hz[2][J][Za] = P[J][Za+1];
        Hz[3][J][Za] = 0.0;
        Hz[4][J][Za] = 0.0;
        EnFluxZ[J][Za] = 0.0;
    }

    //Upper freestream.
    for(K=Za+1; K<=numZ; K++)
    {
        Rho[numR+1][K] = Rho[numR][K];
        n[numR+1][K] = n[numR][K];
        ne[numR+1][K] = ne[numR][K];
        Vr[numR+1][K] = Vr[numR][K];
        Vz[numR+1][K] = Vz[numR][K];
        Bt[numR+1][K] = 0.0;
        E[numR+1][K] = E[numR][K];
        VSq[numR+1][K] = (Vr[numR+1][K]*Vr[numR+1][K])
        +(Vz[numR+1][K]*Vz[numR+1][K]);
        BSq[numR+1][K] = 0.0;
        gam[numR+1][K] = gam[numR][K];
        p[numR+1][K] = p[numR][K];
        P[numR+1][K] = p[numR+1][K];
        EintEl[numR+1][K]= EintEl[numR][K];
        pe[numR+1][K] = pe[numR][K];
        Te[numR+1][K] = Te[numR][K];
        EintH[numR+1][K]= EintH[numR][K];
        ph[numR+1][K] = ph[numR][K];
        Th[numR+1][K] = Th[numR][K];
        Fr[0][numR+1][K] = Rho[numR+1][K]*Vr[numR+1][K];
        Fr[1][numR+1][K] = (Rho[numR+1][K]*Vr[numR+1][K]
        *Vr[numR+1][K])+P[numR+1][K];
        Fr[2][numR+1][K] = Rho[numR+1][K]*Vr[numR+1][K]
        *Vz[numR+1][K];
        Fr[3][numR+1][K] = Bt[numR+1][K]*Vr[numR+1][K];
        Fr[4][numR+1][K] = Vr[numR+1][K]*(E[numR+1][K]
        +P[numR+1][K]);
    }

    //Exit.
    for(J=numR; J>=1; J--)
    {
        Rho[J][numZ+1] = Rho[J][numZ];
        n[J][numZ+1] = n[J][numZ];
        ne[J][numZ+1] = ne[J][numZ];
        Vr[J][numZ+1] = Vr[J][numZ];
        Vz[J][numZ+1] = Vz[J][numZ];
        Bt[J][numZ+1] = 0.0;
        Jencl[J][numZ+1]= 0.0;
        E[J][numZ+1] = E[J][numZ];
        gam[J][numZ+1] = gam[J][numZ];
        EintEl[J][numZ+1]=EintEl[J][numZ];
        pe[J][numZ+1] = pe[J][numZ];
        Te[J][numZ+1] = Te[J][numZ];
        EintH[J][numZ+1]= EintH[J][numZ];
        Th[J][numZ+1] = Th[J][numZ];
        ph[J][numZ+1] = ph[J][numZ];
        Potential[J][numZ+1]= Potential[J][numZ];
        VSq[J][numZ+1] = VSq[J][numZ];
        BSq[J][numZ+1] = 0.0;
        p[J][numZ+1] = p[J][numZ];
        P[J][numZ+1] = p[J][numZ+1];
        Fz[0][J][numZ+1] = Rho[J][numZ+1]*Vz[J][numZ+1];
        Fz[1][J][numZ+1] = Rho[J][numZ+1]*Vr[J][numZ+1]
        *Vz[J][numZ+1];
        Fz[2][J][numZ+1] = (Rho[J][numZ+1]*Vz[J][numZ+1]
        *Vz[J][numZ+1])+P[J][numZ+1];
        Fz[3][J][numZ+1] = Bt[J][numZ+1]*Vz[J][numZ+1];
        Fz[4][J][numZ+1] = Vz[J][numZ+1]*(E[J][numZ+1]
        +P[J][numZ+1]);
    }

    //Centerline.
    for(K=numZ; K>=Zc+1; K--)
    {
        Hr[0][0][K] = 0.0;
        Hr[1][0][K] = P[1][K];
        Hr[2][0][K] = 0.0;
        Hr[3][0][K] = 0.0;
        Hr[4][0][K] = 0.0;
        EnFluxR[0][K] = 0.0;
        SourceConv[0][0][K] = 0.0;
        SourceConv[1][0][K] = 0.0;
        SourceConv[2][0][K] = 0.0;
        SourceConv[3][0][K] = 0.0;
        SourceConv[4][0][K] = 0.0;
    }

    //Cathode tip.
    for(J=1; J<=Rc; J++)
    {
        Hz[0][J][Zc] = 0.0;
        Hz[1][J][Zc] = 0.0;
        Hz[2][J][Zc] = P[J][Zc+1];
        Hz[3][J][Zc] = 0.0;
        Hz[4][J][Zc] = 0.0;
    }

    //Cathode outer.
    for(K=Zc; K>=1; K--)
    {
        Hr[0][Rc][K] = 0.0;
        Hr[1][Rc][K] = Hr[1][Rc+1][K];
        // Hr[1][Rc][K] = P[Rc+1][K];
        Hr[2][Rc][K] = 0.0;
        Hr[3][Rc][K] = 0.0;
        Hr[4][Rc][K] = 0.0;
        EnFluxR[Rc][K] = 0.0;
        SourceConv[0][Rc][K] = 0.0;
        SourceConv[1][Rc][K] = 0.0;
        SourceConv[2][Rc][K] = 0.0;
        SourceConv[3][Rc][K] = 0.0;
        SourceConv[4][Rc][K] = 0.0;
    }
}//End of function

//-------------------------------------
void DissBound()
{
    //Backplate
    if(fabs(Jback) < Jmax)
    Vappl += 0.0005;
    else
    Vappl -= 0.0005;

    for(J=Rc+1; J<=Ra; J++)
    {
        DissipZ[0][J][0] = 0.0;
        DissipZ[1][J][0] = 0.0;
        DissipZ[2][J][0] = 0.0;

        //Total (inductive+resistive) voltage drop at the backplate.
        DissipZ[3][J][0] = Vappl/(log(Ranode/Rcathode)
        *((J-0.5)*deltaR));
        DissipZ[4][J][0] = -((-Vappl/(log(Ranode/Rcathode)
        *((J-0.5)*deltaR)) )
        -(Vz[J][0]*Bt[J][1]))*Bt[J][1]/Mu;

        //Boundary conditions for the species energy equations.
        Te[J][0] = 15000.0;
        Th[J][0] = 15000.0;
        Vz[J][0] = sqrt(gamConst*kBoltz*(Te[J][0]+Th[J][0])/mAr);
        VSq[J][0] = Vz[J][0]*Vz[J][0];
        Rho[J][0] = MassFlowRate/(Vz[J][0]*PI*(pow(Ranode,2.0)
        -pow(Rcathode,2.0)));
        n[J][0] = Rho[J][0]/mAr;
        ne[J][0] = 0.99*n[J][0];
        pe[J][0] = ne[J][0]*kBoltz*Te[J][0];
        ph[J][0] = n[J][0]*kBoltz*Th[J][0];
        p[J][0] = pe[J][0]+ph[J][0];

        //The field diffuses instantaneously inside an insulator.
        Bt[J][0] = Bt[J][1];
        BSq[J][0] = Bt[J][0]*Bt[J][0];
        EintEl[J][0] = pe[J][0]/(gamConst-1);
        EnFluxZ[J][0] = Vz[J][0]*(EintEl[J][0]+pe[J][0]);
    }

    //Anode inner.
    for(K=0; K<=Za; K++)
    {
        DissipR[0][Ra][K] = 0.0;
        DissipR[1][Ra][K] = 0.0;
        DissipR[2][Ra][K] = 0.0;

        //Bt[Ra+1][K]=0.
        DissipR[3][Ra][K] = ResTungst*((-Bt[Ra][K]/deltaR)
        +(Bt[Ra][K]/((Ra-0.5)*deltaR)))/Mu;
        ElecTemp = 2500.0;
        ElCondR[Ra][K] = kTherm[Ra][K]*(ElecTemp-Te[Ra][K])/deltaR;
        IonCondR[Ra][K] = kIon[Ra][K]*(ElecTemp-Th[Ra][K])/deltaR;
        ThermCondR[Ra][K] = ElCondR[Ra][K] + IonCondR[Ra][K];

        //Bt[Ra+1][K]=0.
        DissipR[4][Ra][K] = (DissipR[3][Ra][K]*0.5*Bt[Ra][K]/Mu)
        + ThermCondR[Ra][K];
        EnFluxR[Ra][K] = 0.0;
    }

    //Anode tip.
    for(J=Ra+1; J<=numR; J++)
    {
        DissipZ[0][J][Za]= 0.0;
        DissipZ[1][J][Za]= 0.0;
        DissipZ[2][J][Za]= 0.0;

        //Bt[J][Za] = 0.
        DissipZ[3][J][Za]= ResTungst*(Bt[J][Za+1]/deltaZ)/Mu;
        ElecTemp = 2500.0;
        ElCondZ[J][Za] = kTherm[J][Za+1]*(Te[J][Za+1]-ElecTemp)
        /deltaZ;
        IonCondZ[J][Za] = kIon[J][Za+1]*(Th[J][Za+1]-ElecTemp)
        /deltaZ;
        ThermCondZ[J][Za]= ElCondZ[J][Za] + IonCondZ[J][Za];

        //Bt[J][Za] = 0.
        DissipZ[4][J][Za]= (DissipZ[3][J][Za]*0.5*Bt[J][Za+1]/Mu)
        +ThermCondZ[J][Za];
        EnFluxZ[J][Za] = 0.0;
    }

    //Upper freestream.
    for(K=Za+1; K<=numZ; K++)
    {
        Rho[numR+1][K] = Rho[numR][K];
        n[numR+1][K] = n[numR][K];
        ne[numR+1][K] = ne[numR][K];
        Vr[numR+1][K] = Vr[numR][K];
        Vz[numR+1][K] = Vz[numR][K];
        Bt[numR+1][K] = 0.0;
        E[numR+1][K] = E[numR][K];
        VSq[numR+1][K] = (Vr[numR+1][K]*Vr[numR+1][K])
        +(Vz[numR+1][K]*Vz[numR+1][K]);
        BSq[numR+1][K] = 0.0;
        p[numR+1][K] = p[numR][K];
        EintEl[numR+1][K]= EintEl[numR][K];
        pe[numR+1][K] = pe[numR][K];
        Te[numR+1][K] = Te[numR][K];
        ph[numR+1][K] = ph[numR][K];
        Th[numR+1][K] = Th[numR][K];
        kTherm[numR+1][K]= kTherm[numR][K];
        kIon[numR+1][K]= kIon[numR][K];
        Res[numR+1][K] = Res[numR][K];
        EintEl[numR+1][K]= EintEl[numR][K];
        FrEn[numR+1][K]= Vr[numR+1][K]*(EintEl[numR+1][K]
        +pe[numR+1][K]);
    }

    //Exit.
    for(J=numR; J>=1; J--)
    {
        Rho[J][numZ+1] = Rho[J][numZ];
        n[J][numZ+1] = n[J][numZ];
        ne[J][numZ+1] = ne[J][numZ];
        Vr[J][numZ+1] = Vr[J][numZ];
        Vz[J][numZ+1] = Vz[J][numZ];
        Bt[J][numZ+1] = 0.0;
        //This is for the sake of the output file.
        Jencl[J][numZ+1]= 0.0;
        E[J][numZ+1] = E[J][numZ];
        gam[J][numZ+1] = gam[J][numZ];
        EintEl[J][numZ+1]=EintEl[J][numZ];
        pe[J][numZ+1] = pe[J][numZ];
        Te[J][numZ+1] = Te[J][numZ];
        EintH[J][numZ+1]= EintH[J][numZ];
        Th[J][numZ+1] = Th[J][numZ];
        ph[J][numZ+1] = ph[J][numZ];
        Potential[J][numZ+1]= Potential[J][numZ];
        VSq[J][numZ+1] = VSq[J][numZ];
        BSq[J][numZ+1] = 0.0;
        p[J][numZ+1] = p[J][numZ];
        P[J][numZ+1] = p[J][numZ+1];
        E[J][numZ+1] = E[J][numZ];
        gam[J][numZ+1] = gam[J][numZ];
        EintEl[J][numZ+1]=EintEl[J][numZ];
        pe[J][numZ+1] = pe[J][numZ];
        Te[J][numZ+1] = Te[J][numZ];
        EintH[J][numZ+1]= EintH[J][numZ];
        Th[J][numZ+1] = Th[J][numZ];
        ph[J][numZ+1] = ph[J][numZ];
        Potential[J][numZ+1]= Potential[J-1][numZ]
        +((0.5*(DissipZ[3][J][numZ]
        +DissipZ[3][J][K-1]))-
        (Vz[J][numZ]*Bt[J][numZ]));
        EIcollFreq[J][numZ+1]= EIcollFreq[J][numZ];
        kTherm[J][numZ+1]= kTherm[J][numZ];
        kIon[J][numZ+1] = kIon[J][numZ];
        Res[J][numZ+1] = Res[J][numZ];
        FzEn[J][numZ+1] = Vz[J][numZ+1]*(EintEl[J][numZ+1]
        +pe[J][numZ+1]);
    }

    //Centerline.
    for(K=numZ; K>=Zc+1; K--)
    {
        DissipR[0][0][K] = 0.0;
        DissipR[1][0][K] = 0.0;
        DissipR[2][0][K] = 0.0;

        //Symmetry implies Bt[1][K] = -Bt[-1][K].
        if(fabs(Bt[1][K]) > eps)
        DissipR[3][0][K]= Res[1][K]*(4*Bt[1][K]/deltaR)/Mu;
        else
        DissipR[3][0][K]= 0.0;
        DissipR[4][0][K] = 0.0;
        SourceDiss[0][0][K] = 0.0;
        SourceDiss[1][0][K] = 0.0;
        SourceDiss[2][0][K] = 0.0;
        SourceDiss[3][0][K] = 0.0;
        SourceDiss[4][0][K] = 0.0;
        EnSourceR[0][K] = 0.0;
    }

    //Cathode tip.
    for(J=1; J<=Rc; J++)
    {
        DissipZ[0][J][Zc] = 0.0;
        DissipZ[1][J][Zc] = 0.0;
        DissipZ[2][J][Zc] = 0.0;
        //Bt[J][Zc] = 0.
        DissipZ[3][J][Zc] = ResTungst*(Bt[J][Zc+1]/deltaZ)/Mu;
        ElecTemp = 2500.0;
        ElCondZ[J][Zc] = kTherm[J][Zc+1]*(Te[J][Zc+1]-ElecTemp)
        /deltaZ;
        IonCondZ[J][Zc] = kIon[J][Zc+1]*(Th[J][Zc+1]-2500)
        /deltaZ;
        ThermCondZ[J][Zc] = ElCondZ[J][Zc] + IonCondZ[J][Zc];
        //Bt[J][Zc] = 0.
        DissipZ[4][J][Zc] = (DissipZ[3][J][Zc]*0.5*Bt[J][Zc+1]/Mu)
        +ThermCondZ[J][Zc];
        EnFluxZ[J][Zc] = 0.0;
    }
    
    //Cathode outer.
    for(K=Zc; K>=0; K--)
    {
        DissipR[0][Rc][K] = 0.0;
        DissipR[1][Rc][K] = 0.0;
        DissipR[2][Rc][K] = 0.0;

        //Bt[Rc][K] = 0.
        DissipR[3][Rc][K] = ResTungst*(Bt[Rc+1][K]/deltaR)/Mu;
        ElecTemp = 2500.0;
        ElCondR[Rc][K] = kTherm[Rc+1][K]*(Te[Rc+1][K]-ElecTemp)/deltaR;
        IonCondR[Rc][K]= kIon[Rc+1][K]*(Th[Rc+1][K]-5000)/deltaR;
        ThermCondR[Rc][K] = ElCondR[Rc][K] + IonCondR[Rc][K];
        
        //Bt[Rc][K] = 0.
        DissipR[4][Rc][K] = (DissipR[3][Rc][K]*0.5*Bt[Rc+1][K]/Mu)
        +ThermCondR[Rc][K];
        SourceDiss[0][Rc][K] = 0.0;
        SourceDiss[1][Rc][K] = 0.0;
        SourceDiss[2][Rc][K] = 0.0;
        SourceDiss[3][Rc][K] = 0.0;
        SourceDiss[4][Rc][K] = DissipR[4][Rc][K]/((Rc+1)*deltaR);
        EnSourceR[Rc][K] = 0.0;
    }  
}
//-------------------------------------