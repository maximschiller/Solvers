/* INITIALISING THE ALL PARAMETERS AND BOUNDARIES */

#include "Vars.h"


void SetInitial()
{
    ConvBound();
    EvalParam();
    DissBound();
    for(J=Rc+1; J<=Ra; J++)
    {
        aSq[J][0] = gamConst*p[J][0]/Rho[J][0];
        CmSq[J][0] = aSq[J][0]+(Bt[J][0]*Bt[J][0]
        /(Mu*Rho[J][0]));
        Cf[J][0] = sqrt(CmSq[J][0]);
        //Now, Cf = Cfz. So Cfz is no longer computed
        //Largest eigenvalue
        SpecRadZ[J][0] = fabs(Vz[J][0]+Cf[J][0]);
    }
    
    RhoIn= Rho[Rc+25][1];
    VzIn = Vz[Rc+25][1];
    peIn = pe[Rc+25][1];
    TeIn = Te[Rc+25][1];
    phIn = ph[Rc+25][1];
    ThIn = Th[Rc+25][1];
}
