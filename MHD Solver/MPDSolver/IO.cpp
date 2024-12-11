/* GETTING THE INPUTS */

#include <fstream>
#include "Vars.h"

double MassFlowRate, Jmax, RAD, Ranode, Rcathode, Z, AnodeLength, CathodeLength, numR, numZ, TotTime;

void GetInput()
{
    //Accepting input from file...
    std::ifstream in("INPUT.DAT");
    in >> MassFlowRate >> Jmax >> RAD >> Ranode >> Rcathode
       >> Z >> AnodeLength >> CathodeLength >> numR >> numZ
       >> TotTime ;
    Grid();
}

void WriteFile()
{
    //Temporary output
    if(numSteps%1000 == 0)
    {
        std::ofstream OutTemp("OUT-Temp.DAT");
        OutTemp << " \"z\", \"r\", \"ne\", \"n1\", \"n2\", \"n3\", "
        << " \"nA\", \"Vr\", \"Vz\", \"Bt\", \"Jencl\", "
        << " \"p\", \"Te\",\"Th\", \"gam\", \"Pot\" "
        << "\n";

        for(Zone = 1; Zone <= 3; Zone++)
        {
            for (J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
            {
                for(K=Zmin[Zone]+1; K<=Zmax[Zone]; K++)
                {
                    OutTemp<<(K-1)*deltaZ<<"\t"<<std::setw(7)
                        <<(J-1)*deltaR<<"\t"<<std::setw(7)
                        <<ne[J][K]<<"\t"<<std::setw(7)<<ni[J][K]
                        <<"\t"<<std::setw(7)<<nii[J][K]<<"\t"<<std::setw(7)
                        <<niii[J][K]<<"\t"<<std::setw(7)<<nA[J][K]
                        <<"\t"<<std::setw(7)<<Vr[J][K]<<"\t"<<std::setw(7)
                        <<Vz[J][K]<<"\t"<<std::setw(7)<<Bt[J][K]
                        <<"\t"<<std::setw(7)<<Jencl[J][K]<<"\t"<<std::setw(7)
                        <<p[J][K]<<"\t"<<std::setw(7)<<Te[J][K]<<"\t"
                        <<std::setw(7)<<Th[J][K]<<"\t"<<std::setw(7)
                        <<gam[J][K]<<"\t"<<std::setw(7)
                        <<Potential[J][K]<<std::endl;
                }
            }
        }
    }
    
    //Writing the output.
    //The output file is formatted to be an input file for TECPLOT.
    //Note that even though the variable are stored at cell centers,
    //the output values for (z,r) correspond to the lowest indexed
    //corner of the cell. So, variables at point (J,K), which is
    //(j-0.5, k-0.5), are identified by (j,k).
    //For more information, read TECPLOT User’s Manual.
    if (t>=TotTime)
    {
        std::cerr<<"Completed "<<numSteps<<" time steps in t = "
            <<t<<" seconds; deltaT = "<<deltaT<<std::endl;
        std::ofstream OutFinal("OUTPUT.DAT");

        for(Zone = 1; Zone <= 3; Zone++)
        {
        OutFinal<<"VARIABLES = "<<" \"z\", \"r\", \"ne\", \"n1\","
                <<" \"n2\", \"n3\", \"nA\", \"Vr\","
                <<" \"Vz\", \"Bt\", \"Jencl\","
                <<" \"p\", \"Te\", \"Th\", \"gam\","
                <<" \"Pot\" " <<"\n"
                <<"ZONE I="<<Zmax[Zone]-Zmin[Zone]+1<<", J="
                <<Rmax[Zone]-Rmin[Zone]<<", F=POINT)"<<"\n";

            for (J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
            {
                for(K=Zmin[Zone]+1; K<=Zmax[Zone]+1; K++)
                {
                    OutFinal<<(K-1)*deltaZ<<"\t"<<std::setw(7)
                    <<(J-1)*deltaR<<"\t"<<std::setw(7)<<ne[J][K]
                    <<"\t"<<std::setw(7)<<ni[J][K]<<"\t"<<std::setw(7)
                    <<nii[J][K]<<"\t"<<std::setw(7)<<niii[J][K]
                    <<"\t"<<std::setw(7)<<nA[J][K]<<"\t"<<std::setw(7)
                    <<Vr[J][K]<<"\t"<<std::setw(7)<<Vz[J][K]
                    <<"\t"<<std::setw(7)<<Bt[J][K]<<"\t"<<std::setw(7)
                    <<Jencl[J][K]<<"\t"<<std::setw(7)<<p[J][K]
                    <<"\t"<<std::setw(7)<<Te[J][K]<<"\t"<<std::setw(7)
                    <<Th[J][K]<<"\t"<<std::setw(7)<<gam[J][K]<<"\t"
                    <<std::setw(7)<<Potential[J][K]<<std::endl;
                }
            }
        }

        TerminalChars();
        std::ofstream OutTerm("THRUST.DAT");
        OutTerm<<"Isp = "<<Isp<<" s"<<std::endl<<"Tupstream = "<<Tupstream
        <<" N"<<std::endl<<"Texit = "<<Texit<<" N"<<std::endl
        <<"Tinlet = "<<-Tinlet<<" N"<<std::endl<<"Thrust = "
        <<Thrust<<" N"<<std::endl<<"Maecker Thrust = "<<MaeckerT
        <<" N "<<std::endl;
        //All the variables are written to a data file, so that they can
        //be reloaded.
        std::ofstream OutNe("NeOut.DAT");
        std::ofstream OutNo("NoOut.DAT");
        std::ofstream OutTe("TeOut.DAT");
        std::ofstream Outh("PhOut.DAT");
        std::ofstream OutPh("PhOut.DAT");
        std::ofstream OutTh("ThOut.DAT");
        std::ofstream OutGam("GamOut.DAT");
        std::ofstream OutVr("VrOut.DAT");
        std::ofstream OutVz("VzOut.DAT");
        std::ofstream OutBt("BtOut.DAT");

        //The variables are stored at cell centers. So, ’J’ refers to ’J-0.5’
        //and ’K’ refers to ’K-0.5’.
        for(Zone = 1; Zone <= 3; Zone++)
        {
            for(J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
            {
                for(K=Zmin[Zone]+1; K<=Zmax[Zone]; K++)
                {
                    OutNe<<ne[J][K]<<std::endl;
                    OutNo<<n[J][K]<<std::endl;
                    OutTe<<Te[J][K]<<std::endl;
                    OutPh<<ph[J][K]<<std::endl;
                    OutTh<<Th[J][K]<<std::endl;
                    OutGam<<gam[J][K]<<std::endl;
                    OutVr<<Vr[J][K]<<std::endl;
                    OutVz<<Vz[J][K]<<std::endl;
                    OutBt<<Bt[J][K]<<std::endl;
                }
            }
        }
        std::ofstream OutVolt("Voltage.DAT");
        OutVolt<<Vappl;
    }//end of "TotTime" if...the else follows.
    else
    {
        if ( (t>=9*TotTime/10) &&
        (t<=((9*TotTime/10)+(MaxStep*deltaT))) )
        {
            std::cerr<<"Completed "<<numSteps<<" time steps in t = "<<t
            <<" seconds; deltaT = "<<deltaT<<std::endl;
            std::ofstream Out90("90PERCENT.DAT");
            for(Zone = 1; Zone <= 3; Zone++)
            {
                Out90<<"VARIABLES = "<<" \"z\", \"r\", \"ne\", \"n1\","
                " \"n2\", \"n3\", \"nA\", \"Vr\", \"Vz\", \"Bt\","
                " \"Jencl\", \"p\", \"Te\", \"Th\" "<<"\n"
                <<"ZONE I="<<Zmax[Zone]-Zmin[Zone]+1<<", J="
                <<Rmax[Zone]-Rmin[Zone]<<", F=POINT"<<"\n";

                for (J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
                {
                    for(K=Zmin[Zone]+1; K<=Zmax[Zone]+1; K++)
                    {
                        Out90<<(K-1)*deltaZ<<"\t"<<std::setw(7)
                        <<(J-1)*deltaR<<"\t"<<std::setw(7)<<ne[J][K]
                        <<"\t"<<std::setw(7)<<ni[J][K]<<"\t"<<std::setw(7)
                        <<nii[J][K]<<"\t"<<std::setw(7)<<niii[J][K]
                        <<"\t"<<std::setw(7)<<nA[J][K]<<"\t"<<std::setw(7)
                        <<Vr[J][K]<<"\t"<<std::setw(7)<<Vz[J][K]<<"\t"
                        <<std::setw(7)<<Bt[J][K]<<"\t"<<std::setw(7)
                        <<Jencl[J][K]<<"\t"<<std::setw(7)<<p[J][K]
                        <<"\t"<<std::setw(7)<<Te[J][K]<<"\t"<<std::setw(7)
                        <<Th[J][K]<<std::endl;
                    }
                }
            }
        }//end of "90PERCENT" if...the else follows.
        else
        {
            if ( (t>=8*TotTime/10) &&
            (t<=((8*TotTime/10)+(MaxStep*deltaT))) )
            {
                std::cerr<<"Completed "<<numSteps<<" time steps in t = "
                <<t<<" seconds; deltaT = "<<deltaT<<std::endl;
                std::ofstream Out80("80PERCENT.DAT");
                for(Zone = 1; Zone <= 3; Zone++)
                {

                    Out80<<"VARIABLES = "<<" \"z\", \"r\", \"ne\", \"n1\","
                    " \"n2\", \"n3\", \"nA\", \"Vr\", \"Vz\", \"Bt\","
                    " \"Jencl\", \"p\", \"Te\", \"Th\" "<<"\n"
                    <<"ZONE I="<<Zmax[Zone]-Zmin[Zone]+1<<", J="
                    <<Rmax[Zone]-Rmin[Zone]<<", F=POINT"<<"\n";

                    for (J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
                    {
                        for(K=Zmin[Zone]+1; K<=Zmax[Zone]+1; K++)
                        {
                            Out80<<(K-1)*deltaZ<<"\t"<<std::setw(7)
                            <<(J-1)*deltaR<<"\t"<<std::setw(7)
                            <<ne[J][K]<<"\t"<<std::setw(7)<<ni[J][K]
                            <<"\t"<<std::setw(7)<<nii[J][K]<<"\t"
                            <<std::setw(7)<<niii[J][K]<<"\t"<<std::setw(7)
                            <<nA[J][K]<<"\t"<<std::setw(7)<<Vr[J][K]
                            <<"\t"<<std::setw(7)<<Vz[J][K]<<"\t"
                            <<std::setw(7)<<Bt[J][K]<<"\t"<<std::setw(7)
                            <<Jencl[J][K]<<"\t"<<std::setw(7)
                            <<p[J][K]<<"\t"<<std::setw(7)<<Te[J][K]
                            <<"\t"<<std::setw(7)<<Th[J][K]<<std::endl;
                        }
                    }
                }
            }//end of "80PERCENT" if...the else follows.
            else
            {
                if ( (t>=7*TotTime/10) &&
                (t<=((7*TotTime/10)+(MaxStep*deltaT))) )
                {
                    std::cerr<<"Completed "<<numSteps<<" time steps in t = "
                    <<t<<" seconds; deltaT = "<<deltaT<<std::endl;
                    std::ofstream Out70("70PERCENT.DAT");
                    for(Zone = 1; Zone <= 3; Zone++)
                    {

                        Out70<<"VARIABLES = "<<" \"z\", \"r\", \"ne\", \"n1\","
                        " \"n2\", \"n3\", \"nA\", \"Vr\", \"Vz\", \"Bt\","
                        " \"Jencl\", \"p\", \"Te\", \"Th\" "<<"\n"
                        <<"ZONE I="<<Zmax[Zone]-Zmin[Zone]+1<<", J="
                        <<Rmax[Zone]-Rmin[Zone]<<", F=POINT"<<"\n";

                        for (J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
                        {
                            for(K=Zmin[Zone]+1; K<=Zmax[Zone]+1; K++)
                            {
                                Out70<<(K-1)*deltaZ<<"\t"<<std::setw(7)
                                <<(J-1)*deltaR<<"\t"<<std::setw(7)
                                <<ne[J][K]<<"\t"<<std::setw(7)
                                <<ni[J][K]<<"\t"<<std::setw(7)
                                <<nii[J][K]<<"\t"<<std::setw(7)
                                <<niii[J][K]<<"\t"<<std::setw(7)
                                <<nA[J][K]<<"\t"<<std::setw(7)
                                <<Vr[J][K]<<"\t"<<std::setw(7)
                                <<Vz[J][K]<<"\t"<<std::setw(7)
                                <<Bt[J][K]<<"\t"<<std::setw(7)
                                <<Jencl[J][K]<<"\t"<<std::setw(7)
                                <<p[J][K]<<"\t"<<std::setw(7)
                                <<Te[J][K]<<"\t"<<std::setw(7)
                                <<Th[J][K]<<std::endl;
                            }
                        }
                    }
                }//end of "70PERCENT" if...the else follows.
                //The tab allignments are intentionally moved off.
                else
                {
                    if( (t>=6*TotTime/10) &&
                    (t<=((6*TotTime/10)+(MaxStep*deltaT))) )
                    {
                        std::cerr<<"Completed "<<numSteps<<" time steps in t = "
                        <<t<<" seconds; deltaT = "<<deltaT<<std::endl;
                        std::ofstream Out60("60PERCENT.DAT");
                        for(Zone = 1; Zone <= 3; Zone++)
                        {

                            Out60<<"VARIABLES = "<<" \"z\", \"r\", \"ne\", \"n1\","
                            " \"n2\", \"n3\", \"nA\", \"Vr\", \"Vz\", \"Bt\","
                            " \"Jencl\", \"p\", \"Te\", \"Th\" "<<"\n"
                            <<"ZONE I="<<Zmax[Zone]-Zmin[Zone]+1<<", J="
                            <<Rmax[Zone]-Rmin[Zone]<<", F=POINT"<<"\n";
                    
                            for (J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
                            {
                                for(K=Zmin[Zone]+1; K<=Zmax[Zone]+1; K++)
                                {
                                    Out60<<(K-1)*deltaZ<<"\t"<<std::setw(7)
                                    <<(J-1)*deltaR<<"\t"<<std::setw(7)<<ne[J][K]
                                    <<"\t"<<std::setw(7)<<ni[J][K]<<"\t"<<std::setw(7)
                                    <<nii[J][K]<<"\t"<<std::setw(7)<<niii[J][K]
                                    <<"\t"<<std::setw(7)<<nA[J][K]<<"\t"<<std::setw(7)
                                    <<Vr[J][K]<<"\t"<<std::setw(7)<<Vz[J][K]
                                    <<"\t"<<std::setw(7)<<Bt[J][K]<<"\t"<<std::setw(7)
                                    <<Jencl[J][K]<<"\t"<<std::setw(7)<<p[J][K]
                                    <<"\t"<<std::setw(7)<<Te[J][K]<<"\t"<<std::setw(7)
                                    <<Th[J][K]<<std::endl;
                                }
                            }
                        }
                    }//end of "60PERCENT" if...the else follows.
                    else
                    {
                        if( (t>=5*TotTime/10) &&
                        (t<=((5*TotTime/10)+(MaxStep*deltaT))) )
                        {
                            std::cerr<<"Completed "<<numSteps<<" time steps in t = "
                            <<t<<" seconds; deltaT = "<<deltaT<<std::endl;
                            std::ofstream Out50("50PERCENT.DAT");
                            for(Zone = 1; Zone <= 3; Zone++)
                            {

                                Out50<<"VARIABLES = "<<" \"z\", \"r\", \"ne\", \"n1\","
                                " \"n2\", \"n3\", \"nA\", \"Vr\", \"Vz\", \"Bt\","
                                " \"Jencl\", \"p\", \"Te\", \"Th\" "<<"\n"
                                <<"ZONE I="<<Zmax[Zone]-Zmin[Zone]+1<<", J="
                                <<Rmax[Zone]-Rmin[Zone]<<", F=POINT"<<"\n";
                                
                                for (J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
                                {
                                    for(K=Zmin[Zone]+1; K<=Zmax[Zone]+1; K++)
                                    {
                                        Out50<<(K-1)*deltaZ<<"\t"<<std::setw(7)
                                        <<(J-1)*deltaR<<"\t"<<std::setw(7)
                                        <<ne[J][K]<<"\t"<<std::setw(7)<<ni[J][K]
                                        <<"\t"<<std::setw(7)<<nii[J][K]<<"\t"
                                        <<std::setw(7)<<niii[J][K]<<"\t"<<std::setw(7)
                                        <<nA[J][K]<<"\t"<<std::setw(7)<<Vr[J][K]
                                        <<"\t"<<std::setw(7)<<Vz[J][K]<<"\t"
                                        <<std::setw(7)<<Bt[J][K]<<"\t"<<std::setw(7)
                                        <<Jencl[J][K]<<"\t"<<std::setw(7)
                                        <<p[J][K]<<"\t"<<std::setw(7)<<Te[J][K]
                                        <<"\t"<<std::setw(7)<<Th[J][K]<<std::endl;
                                    }
                                }
                            }
                        }//end of "50PERCENT" if...the else follows.
                        else
                        {
                            if( (t>=4*TotTime/10) &&
                            (t<=((4*TotTime/10)+(MaxStep*deltaT))) )
                            {
                                std::cerr<<"Completed "<<numSteps<<" time steps in t = "
                                <<t<<" seconds; deltaT = "<<deltaT<<std::endl;
                                std::ofstream Out40("40PERCENT.DAT");
                                
                                for(Zone = 1; Zone <= 3; Zone++)
                                {

                                    Out40<<"VARIABLES = "<<" \"z\", \"r\", \"ne\", \"n1\","
                                    " \"n2\", \"n3\", \"nA\", \"Vr\", \"Vz\", \"Bt\","
                                    " \"Jencl\", \"p\", \"Te\", \"Th\" "<<"\n"
                                    <<"ZONE I="<<Zmax[Zone]-Zmin[Zone]+1<<", J="
                                    <<Rmax[Zone]-Rmin[Zone]<<", F=POINT"<<"\n";

                                    for (J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
                                    {
                                        for(K=Zmin[Zone]+1;K<=Zmax[Zone]+1;K++)
                                        {
                                            Out40<<(K-1)*deltaZ<<"\t"<<std::setw(7)
                                            <<(J-1)*deltaR<<"\t"<<std::setw(7)
                                            <<ne[J][K]<<"\t"<<std::setw(7)
                                            <<ni[J][K]<<"\t"<<std::setw(7)
                                            <<nii[J][K]<<"\t"<<std::setw(7)
                                            <<niii[J][K]<<"\t"<<std::setw(7)
                                            <<nA[J][K]<<"\t"<<std::setw(7)
                                            <<Vr[J][K]<<"\t"<<std::setw(7)
                                            <<Vz[J][K]<<"\t"<<std::setw(7)
                                            <<Bt[J][K]<<"\t"<<std::setw(7)
                                            <<Jencl[J][K]<<"\t"<<std::setw(7)
                                            <<p[J][K]<<"\t"<<std::setw(7)
                                            <<Te[J][K]<<"\t"<<std::setw(7)
                                            <<Th[J][K]<<std::endl;
                                        }
                                    }
                                }
                            }//end of "40PERCENT" if...the else follows.
                            //The tab allignments are intentionally offset.
                            else
                            {
                                if( (t>=3*TotTime/10) &&
                                (t<=((3*TotTime/10)+(MaxStep*deltaT))) )
                                {
                                    std::cerr<<"Completed "<<numSteps<<" time steps in t = "
                                    <<t<<" seconds; deltaT = "<<deltaT<<std::endl;
                                    std::ofstream Out30("30PERCENT.DAT");
                                    for(Zone = 1; Zone <= 3; Zone++)
                                    {

                                        Out30<<"VARIABLES = "<<" \"z\", \"r\", \"ne\", \"n1\","
                                        " \"n2\", \"n3\", \"nA\", \"Vr\", \"Vz\", \"Bt\","
                                        " \"Jencl\", \"p\", \"Te\", \"Th\" "<<"\n"
                                        <<"ZONE I="<<Zmax[Zone]-Zmin[Zone]+1<<", J="
                                        <<Rmax[Zone]-Rmin[Zone]<<", F=POINT"<<"\n";

                                        for (J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
                                        {
                                            for(K=Zmin[Zone]+1; K<=Zmax[Zone]+1; K++)
                                            {
                                                Out30<<(K-1)*deltaZ<<"\t"<<std::setw(7)
                                                <<(J-1)*deltaR<<"\t"<<std::setw(7)<<ne[J][K]
                                                <<"\t"<<std::setw(7)<<ni[J][K]<<"\t"<<std::setw(7)
                                                <<nii[J][K]<<"\t"<<std::setw(7)<<niii[J][K]
                                                <<"\t"<<std::setw(7)<<nA[J][K]<<"\t"<<std::setw(7)
                                                <<Vr[J][K]<<"\t"<<std::setw(7)<<Vz[J][K]<<"\t"
                                                <<std::setw(7)<<Bt[J][K]<<"\t"<<std::setw(7)
                                                <<Jencl[J][K]<<"\t"<<std::setw(7)<<p[J][K]
                                                <<"\t"<<std::setw(7)<<Te[J][K]<<"\t"<<std::setw(7)
                                                <<Th[J][K]<<std::endl;
                                            }
                                        }
                                    }
                                }//end of "30PERCENT" if...the else follows.
                                else
                                {
                                    if( (t>=2*TotTime/10) &&
                                    (t<=((2*TotTime/10)+(MaxStep*deltaT))) )
                                    {
                                        std::cerr<<"Completed "<<numSteps<<" time steps in t = "
                                        <<t<<" seconds; deltaT = "<<deltaT<<std::endl;
                                        std::ofstream Out20("20PERCENT.DAT");
                                        for(Zone = 1; Zone <= 3; Zone++)
                                        {

                                            Out20<<"VARIABLES = "<<" \"z\", \"r\", \"ne\", \"n1\","
                                            " \"n2\", \"n3\", \"nA\", \"Vr\", \"Vz\", \"Bt\","
                                            " \"Jencl\", \"p\", \"Te\", \"Th\" "<<"\n"
                                            <<"ZONE I="<<Zmax[Zone]-Zmin[Zone]+1<<", J="
                                            <<Rmax[Zone]-Rmin[Zone]<<", F=POINT"<<"\n";

                                            for (J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
                                            {
                                                for(K=Zmin[Zone]+1; K<=Zmax[Zone]+1; K++)
                                                {
                                                Out20<<(K-1)*deltaZ<<"\t"<<std::setw(7)
                                                    <<(J-1)*deltaR<<"\t"<<std::setw(7)
                                                    <<ne[J][K]<<"\t"<<std::setw(7)<<ni[J][K]
                                                    <<"\t"<<std::setw(7)<<nii[J][K]<<"\t"
                                                    <<std::setw(7)<<niii[J][K]<<"\t"<<std::setw(7)
                                                    <<nA[J][K]<<"\t"<<std::setw(7)<<Vr[J][K]
                                                    <<"\t"<<std::setw(7)<<Vz[J][K]<<"\t"
                                                    <<std::setw(7)<<Bt[J][K]<<"\t"<<std::setw(7)
                                                    <<Jencl[J][K]<<"\t"<<std::setw(7)<<p[J][K]
                                                    <<"\t"<<std::setw(7)<<Te[J][K]<<"\t"
                                                    <<std::setw(7)<<Th[J][K]<<std::endl;
                                                }
                                            }
                                        }
                                    }//end of "20PERCENT" if...the else follows.
                                    //The tab allignments are intentionally offset.
                                    else
                                    {
                                        if( (t>=TotTime/10) &&
                                        (t<=((TotTime/10)+(MaxStep*deltaT))) )
                                        {
                                            std::cerr<<"Completed "<<numSteps<<" time steps in t = "<<t
                                            <<" seconds; deltaT = "<<deltaT<<std::endl;
                                            std::ofstream Out10("10PERCENT.DAT");
                                            for(Zone = 1; Zone <= 3; Zone++)
                                            {

                                                Out10<<"VARIABLES = "<<" \"z\", \"r\", \"ne\", \"n1\","
                                                " \"n2\", \"n3\", \"nA\", \"Vr\", \"Vz\", \"Bt\","
                                                " \"Jencl\", \"p\", \"Te\", \"Th\" "<<"\n"
                                                <<"ZONE I="<<Zmax[Zone]-Zmin[Zone]+1<<", J="
                                                <<Rmax[Zone]-Rmin[Zone]<<", F=POINT"<<"\n";

                                                for (J=Rmin[Zone]+1; J<=Rmax[Zone]; J++)
                                                {
                                                    for(K=Zmin[Zone]+1; K<=Zmax[Zone]+1; K++)
                                                    {
                                                        Out10<<(K-1)*deltaZ<<"\t"<<std::setw(7)
                                                        <<(J-1)*deltaR<<"\t"<<std::setw(7)<<ne[J][K]
                                                        <<"\t"<<std::setw(7)<<ni[J][K]<<"\t"<<std::setw(7)
                                                        <<nii[J][K]<<"\t"<<std::setw(7)<<niii[J][K]
                                                        <<"\t"<<std::setw(7)<<nA[J][K]<<"\t"<<std::setw(7)
                                                        <<Vr[J][K]<<"\t"<<std::setw(7)<<Vz[J][K]<<"\t"
                                                        <<std::setw(7)<<Bt[J][K]<<"\t"<<std::setw(7)
                                                        <<Jencl[J][K]<<"\t"<<std::setw(7)<<p[J][K]
                                                        <<"\t"<<std::setw(7)<<Te[J][K]<<"\t"<<std::setw(7)
                                                        <<Th[J][K]<<std::endl;
                                                    }
                                                }
                                            }
                                        }//end of 10 percent loop.
                                    }//end of 20 percent loop.
                                }//end of 30 percent loop.
                            }//end of 40 percent loop.
                        }//end of 50 percent loop.
                    }//end of 60 percent loop.
                }//end of 70 percent loop.
            }//end of 80 percent loop.
        }//end of 90 percent loop.
    }//end of total-time loop.
}



