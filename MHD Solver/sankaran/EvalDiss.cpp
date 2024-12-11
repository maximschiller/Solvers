//---------------------------------------------------------
//EvalDiss.CPP; Written by Kameshwaran Sankaran
//---------------------------------------------------------
void EvalDiss()
{
    //Physical dissipation
    //CHAMBER
    for(Zone = 1; Zone <= 3; Zone++)
    {
        for(j=J=Rmin[Zone]+1; j<=Rmax[Zone]; j++,J++)
        {
            for (K=Zmin[Zone]+1; K<=Zmax[Zone]; K++)
            {
                //There is no need to waste time computing in places without magnetic field.
                if(BSq[J][K] > eps)
                {
                    jz[j][K] = (((J+0.5)*Bt[J+1][K])
                    -((J-0.5)*Bt[J][K]))
                    /(Mu*(J+0.5)*deltaR);
                    Ez[j][K] = 0.5*(Res[J][K]+Res[J+1][K])*jz[j][K];
                    EzH[j][K]= (0.5*(jr[J][K]+jr[J][K-1]))
                    *(0.5*((Bt[J][K]/n[J][K])
                    +(Bt[J+1][K]/n[J+1][K]))/q);
                    
                    if((Zone==1)&&(K==1))
                    EzH[j][K] = 0.0;

                    if(Zone==2)
                    {
                        if(Zc <= Za)
                        {
                            if((j<=Rc)&&(K==Zc+1))
                            EzH[j][K] = 0.0;
                        }
                        else
                        {
                            if( ((j>Ra)&&(K==Za+1)) || (j==numR) )
                            EzH[j][K] = 0.0;
                        }
                    }

                    if(Zone==3)
                    {
                        if(Zc <= Za)
                        {
                            if( ((j>Ra)&&(K==Za+1)) || (j==numR) )
                            EzH[j][K] = 0.0;
                        }
                        else
                        {
                            if( ((j<=Rc)&&(K==Zc+1)) || (j==numR) )
                            EzH[j][K] = 0.0;
                        }
                    }
                    
                    DissipR[3][j][K]= Ez[j][K] + EzH[j][K];
                    ElCondR[j][K]= 0.5*(kTherm[J][K]+kTherm[J+1][K])
                    *(Te[J+1][K]-Te[J][K])/deltaR;
                    IonCondR[j][K]= 0.5*(kIon[J][K]+kIon[J+1][K])
                    *(Th[J+1][K]-Th[J][K])/deltaR;
                    ThermCondR[j][K]= ElCondR[j][K] + IonCondR[j][K];
                    DissipR[4][j][K]= (DissipR[3][j][K]*
                    0.5*(Bt[J][K]+Bt[J+1][K])/Mu)
                    + ThermCondR[j][K];
                    SourceDiss[4][j][K] = DissipR[4][j][K]
                    /((j+1)*deltaR);
                }//end of ’if’.
            }//end of ’K’.
        }//end of ’j,J’.
        for(J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
        {
            for (k=K=Zmin[Zone]; k<=Zmax[Zone]-1; k++,K++)
            {
                if(BSq[J][K] > eps)
                {
                    jr[J][k] = -(Bt[J][K+1]-Bt[J][K])/(Mu*deltaZ);
                    Er[J][k] = -0.5*(Res[J][K]+Res[J][K+1])*jr[J][k];
                    ErH[J][k]= -(0.5*(jz[J][K]+jz[J-1][K]))
                    *(0.5*((Bt[J][K+1]/n[J][K+1])
                    +(Bt[J][K]/n[J][K]))/q);

                    if((Zone==1)&&((J==Rc+1)||(J==Ra)))
                    ErH[J][k] = 0.0;

                    if(Zone==2)
                    {
                        if(Zc<=Za)
                        {
                            if((J==1)||(J==Ra))
                            ErH[J][k] = 0.0;
                        }
                        else
                        {
                            if((J==Rc+1)||(J==numR))
                            ErH[J][k] = 0.0;
                        }
                    }

                    if(Zone==3)
                    {
                        //No Hall effect on the boundary cells.
                        if((J==1) || (J==numR) || (k == numZ))
                        ErH[J][k] = 0.0;
                    }

                    DissipZ[3][J][k]= Er[J][k] + ErH[J][k];
                    ElCondZ[J][k] = 0.5*(kTherm[J][K]+kTherm[J][K+1])
                    *(Te[J][K+1]-Te[J][K])/deltaZ;
                    IonCondZ[J][k] = 0.5*(kIon[J][K]+kIon[J][K+1])
                    *(Th[J][K+1]-Th[J][K])/deltaZ;
                    ThermCondZ[J][k] = ElCondZ[J][k] + IonCondZ[J][k];
                    DissipZ[4][J][k]= (DissipZ[3][J][k]
                    * 0.5*(Bt[J][K]+Bt[J][K+1])/Mu)
                    + ThermCondZ[J][k];
                }
            }//end of ’k,K’.
        }//end of ’J’.
    }//end for all zones.

    DissBound();//Computes the dissipative terms at the boundaries.
}
//-------------------------------------
