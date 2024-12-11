/* SCRIPT TAKES THE GRID COORDINATES AND EVALUEATES THE CELL VOLUME, AREA AND SKEWNESS */
#include "Variables.h"
#include <iostream>
#include <cmath>
#include <fstream>

// Declarations
double xave, yave, s1, s2, xVol, yVol, ang1, ang2; // variables used to compute centroid position
double dzeta, deta; // cell centre to cell centre length and face area length
int Ra, Rc, Zc, Za; // cells on which the anode and cathod are in the computational domain
int Rmin[4], Rmax[4], Zmin[4], Zmax[4]; //Vertices of the zone.


// Note that order of storage for components is along rows i and columns j - i is the lengthwise direction and j is in the radial direction
// As a result, the values are printed out with the left curve on top and bottom curve on left
// function to compute the cell centres of each cell for flux calculations
// Calculates the centroid of each cell
void cellCentres(){

    for(j=0;j<jmax;j++){
        for(i=0;i<imax;i++){
        
            // claculating the average x and y positions in each cell given that they are four-sided
            xave = (x[j][i] + x[j+1][i] + x[j][i+1] + x[j+1][i+1])/4;
            yave = (y[j][i] + y[j+1][i] + y[j][i+1] + y[j+1][i+1])/4;

            // Calculating the area in each cell by splitting into two triangles
            // Triangle 1 is top left with angle ang1
            // Triangle 2 is bottom right with angle ang2
            // Angles are calculated by using the law of cosines
            ang1 = acos(pow((x[j+1][i+1]-x[j][i]),2) - pow((x[j+1][i]-x[j][i]),2) - pow((x[j+1][i+1]-x[j+1][i]),2)
                    + 2*(x[j+1][i]-x[j][i])*(x[j+1][i+1]-x[j+1][i]));
            ang2 = acos(pow((x[j+1][i+1]-x[j][i]),2) - pow((x[j][i+1]-x[j][i]),2) - pow((x[j+1][i+1]-x[j][i+1]),2)
                    + 2*(x[j][i+1]-x[j][i])*(x[j+1][i+1]-x[j][i+1]));
            xVol = sin(ang1)*(x[j+1][i]-x[j][i])*(x[j+1][i+1]-x[j+1][i])*0.5
                        + sin(ang2)*(x[j][i+1]-x[j][i])*(x[j+1][i+1]-x[j][i+1])*0.5;

            ang1 = acos(pow((y[j+1][i+1]-y[j][i]),2) - pow((y[j+1][i]-y[j][i]),2) - pow((y[j+1][i+1]-y[j+1][i]),2)
                    + 2*(y[j+1][i]-y[j][i])*(y[j+1][i+1]-y[j+1][i]));
            ang2 = acos(pow((y[j+1][i+1]-y[j][i]),2) - pow((y[j][i+1]-y[j][i]),2) - pow((y[j+1][i+1]-y[j][i+1]),2)
                    + 2*(y[j][i+1]-y[j][i])*(y[j+1][i+1]-y[j][i+1]));
            yVol = sin(ang1)*(y[j+1][i]-y[j][i])*(y[j+1][i+1]-y[j+1][i])*0.5
                        + sin(ang2)*(y[j][i+1]-y[j][i])*(y[j+1][i+1]-y[j][i+1])*0.5;

            xc[j][i] = (xave*xVol)/xVol; 
            yc[j][i] = (yave*yVol)/yVol;
        }
    }

    // std::cout << "after computing grid centers" << std::endl;
    // for(j=0; j<jmax; j++) {
    //     for(i=0; i<imax; i++) {
    //         std::cout << xc[j][i] << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "--------------------------------" << std::endl;
    // for(j=0; j<jmax; j++) {
    //     for(i=0; i<imax; i++) {
    //         std::cout << yc[j][i] << ' ';
    //     }
    //     std::cout << std::endl;
    // } 

     /*--------------------------------------------------------------*/
    // Outputing data so that the mesh can be plotted in python

    std::ofstream outfile_x("centrexvals.txt");
    // Output array values to file in a format that can be easily imported into Python
    for(j=0; j<jmax; j++) {
        for(i=0; i<imax; i++) {
            outfile_x << xc[j][i] << " ";
        }
        outfile_x << std::endl;
    }
    outfile_x.close();

    std::ofstream outfile_y("centreyvals.txt");
    // Output array values to file in a format that can be easily imported into Python
    for(j=0; j<jmax; j++) {
        for(i=0; i<imax; i++) {
            outfile_y << yc[j][i] << " ";
        }
        outfile_y << std::endl;
    }
    outfile_y.close();

}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

// function to compute the volume, area and normal vectors
// furthermore finds tensor products and checks for skewness of mesh cells
void cellGeom(){

    // Using the location of the cell centres to compute the area of the faces between each cell
    // method given by: An Introduction to Computational Fluid Dynamics - THE FINITE VOLUME METHOD Second Edition (H K Versteeg and W Malalasekera)
    // Note that the face areas are for each of the cells - i.e. fAx[0][0] is cell 1, fAx[1][0] is cell 2
    
    // Current compass directions - north upstream, east towards anode, south downstream and west towards cathode
    for(j=0; j<jmax; j++) {
        for(i=0; i<imax; i++) {
            // volume of each cell - check if needed
            // computed using Gauss's formula
            totVol[j][i] = 0.5*( (x[j][i+1]-x[j+1][i])*(y[j+1][i+1]-y[j][i])
                                + (x[j][i]-x[j+1][i+1])*(y[j][i+1]-y[j+1][i]) );



            // face areas north, east, south and west
            fAn[j][i] = sqrt(pow(x[j+1][i]-x[j][i],2) + pow(y[j+1][i]-y[j][i],2));
            fAe[j][i] = sqrt(pow(x[j+1][i+1]-x[j+1][i],2) + pow(y[j+1][i+1]-y[j+1][i],2));
            fAs[j][i] = sqrt(pow(x[j+1][i+1]-x[j][i+1],2) + pow(y[j+1][i+1]-y[j][i+1],2));
            fAw[j][i] = sqrt(pow(x[j][i+1]-x[j][i],2) + pow(y[j][i+1]-y[j][i],2));



            // normal vectors given by face area
            // mormal vectors in i and j direction - north, east, south and west
            nvni[j][i] = (y[j+1][i]-y[j][i])/fAn[j][i];
            nvnj[j][i] = (x[j+1][i]-x[j][i])/fAn[j][i];
            nvei[j][i] = (y[j+1][i+1]-y[j+1][i])/fAe[j][i];
            nvej[j][i] = (x[j+1][i+1]-x[j+1][i])/fAe[j][i]; 
            nvsi[j][i] = (y[j+1][i+1]-y[j][i+1])/fAs[j][i];
            nvsj[j][i] = (x[j+1][i+1]-x[j][i+1])/fAs[j][i];
            nvwi[j][i] = (y[j][i+1]-y[j][i])/fAw[j][i];
            nvwj[j][i] = (x[j][i+1]-x[j][i])/fAw[j][i];



            // NOTE - for ezeta and eta need to check if cells are on the boundaries
            // Recal that location of cell centres are actually one above the i,j counter
            // ezeta (vector parallel to line joining cell centers) north, east, south and west
            // dzeta is distance between cell centres
            // In order to coincide ezeta vecotrs with cell centres need to account for row that is being removed
            if(i>0){ // cells with north faces on boundary cannot have ezetan
                dzeta = sqrt(pow(xc[j][i-1]-xc[j][i],2) + pow(yc[j][i-1]-yc[j][i],2));
                ezetani[j][i] = (xc[j][i-1]-xc[j][i])/dzeta;
                ezetanj[j][i] = (yc[j][i-1]-yc[j][i])/dzeta;
            }
            else{
                ezetani[j][i] = 0;
                ezetanj[j][i] = 0;
            }

            if(j<jmax-1){
                dzeta = sqrt(pow(xc[j+1][i]-xc[j][i],2) + pow(yc[j+1][i]-yc[j][i],2));
                ezetaei[j][i] = (xc[j+1][i]-xc[j][i])/dzeta;
                ezetaej[j][i] = (yc[j+1][i]-yc[j][i])/dzeta;
            }
            else{
                ezetaei[j][i] = 0;
                ezetaej[j][i] = 0;
            }

            if(i<imax-1){
                dzeta = sqrt(pow(xc[j][i+1]-xc[j][i],2) + pow(yc[j][i+1]-yc[j][i],2));
                ezetasi[j][i] = (xc[j][i+1]-xc[j][i])/dzeta;
                ezetasj[j][i] = (yc[j][i+1]-yc[j][i])/dzeta;
            }
            else{
                ezetasi[j][i] = 0;
                ezetasj[j][i] = 0;
            }

            if(j>0){ // cells with west faces on boundary cannot have ezetaw
                dzeta = sqrt(pow(xc[j-1][i]-xc[j][i],2) + pow(yc[j-1][i]-yc[j][i],2));
                ezetawi[j][i] = (xc[j-1][i] - xc[j][i])/dzeta;
                ezetawj[j][i] = (yc[j-1][i] - yc[j][i])/dzeta;          
            }
            else{
                ezetawi[j][i] = 0;
                ezetawj[j][i] = 0;
            }



            // eeta (vector parallel to cell boundary) north, east, south and west
            // deta is distance of face areas
            deta = sqrt(pow(x[j+1][i]-x[j][i],2) + pow(y[j+1][i]-y[j][i],2));
            eetani[j][i] = (x[j+1][i]-x[j][i])/deta;
            eetanj[j][i] = (y[j+1][i]-y[j][i])/deta;

            deta = sqrt(pow(x[j+1][i+1]-x[j+1][i],2) + pow(y[j+1][i+1]-y[j+1][i],2));
            eetaei[j][i] = (x[j+1][i+1]-x[j+1][i])/deta;
            eetaej[j][i] = (y[j+1][i+1]-y[j+1][i])/deta;

            deta = sqrt(pow(x[j+1][i+1]-x[j][i+1],2) + pow(y[j+1][i+1]-y[j][i+1],2));
            eetasi[j][i] = (x[j+1][i+1]-x[j][i+1])/deta;
            eetasj[j][i] = (y[j+1][i+1]-y[j][i+1])/deta;

            deta = sqrt(pow(x[j][i+1]-x[j][i],2) + pow(y[j][i+1]-y[j][i],2));
            eetawi[j][i] = (x[j][i+1]-x[j][i])/deta;
            eetawj[j][i] = (y[j][i+1]-y[j][i])/deta;

            // need to calculate the angle between the normal and face vector for flux calculation
            thetan[j][i] = atan( ((ezetani[j][i]*eetani[j][i]) + (ezetanj[j][i]*eetanj[j][i]))/
                            ((nvni[j][i]*ezetani[j][i]) + (nvnj[j][i]*ezetanj[j][i])) );
            thetae[j][i] = atan( ((ezetaei[j][i]*eetaei[j][i]) + (ezetaej[j][i]*eetaej[j][i]))/
                            ((nvei[j][i]*ezetaei[j][i]) + (nvej[j][i]*ezetaej[j][i])) );
            thetas[j][i] = atan( ((ezetasi[j][i]*eetasi[j][i]) + (ezetasj[j][i]*eetasj[j][i]))/
                            ((nvsi[j][i]*ezetasi[j][i]) + (nvsj[j][i]*ezetasj[j][i])) );
            thetaw[j][i] = atan( ((ezetawi[j][i]*eetawi[j][i]) + (ezetawj[j][i]*eetawj[j][i]))/
                            ((nvwi[j][i]*ezetawi[j][i]) + (nvwj[j][i]*ezetawj[j][i])) );

        }
    }

    // std::cout << "--------------------------------" << std::endl;
    // std::cout << "checking vectors" << std::endl;
    // for(j=0; j<jmax; j++) {
    //     for(i=0; i<imax; i++) {
    //         std::cout << eetani[j][i] << ' ';
    //     }
    //     std::cout << std::endl;
    // } 

}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

// Computes the skewness of the cells in order to determine mesh quality
// may not potentiall need since mesh is fairly structured
void gridError(){

    // check the angle of each cell - using fluid mechanics youtube video
    // require that cell angle is less than 70 deg
    // for(j=0; j<jmax; j++) {
    //     for(i=0; i<imax; i++) {
            // thetac[j][i] = 
    //     }
    // }

}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

// function to compute the dimensions of the thruster based on the mesh
void dimensions(){

    // computing the bounds of cell centre counters
    Jmax = sizeof(xc)/sizeof(xc[0]);
    Imax = sizeof(xc[0])/sizeof(xc[0][0]); 


    // computing the values for the radius of the anode and cathode
    // This run a check that if the value on the either extremes in the radial direction
    // is above or below some value then this is the anode/cathode
    int LenRc, LenRa;
    for(i=0;i<=imax;i++){
        if(y[0][i] >= 0){ // currently set to =0 but will update later since the cathode will not be at 0
            Rcathode[i] = y[0][i];
            LenRc = i;
        }
        if(y[jmax][i] < jmax){ // set the bound to jmax as will eventually have the thrster domain exit equal to jmax
            Ranode[i] = y[jmax][i];
            LenRa = i;
        }
    }

    std::cout << "--------------------------------" << std::endl;
    std::cout << "printing positions Rcathode at centre cells" << std::endl;
    for(i=0;i<LenRc;i++){
        Rcathode[i] = (Rcathode[i+1]+Rcathode[i])/2;
        std::cout << Rcathode[i] << ' ';
    }
    std::cout << std::endl;

    std::cout << "printing positions Ranode at centre cells" << std::endl;
    for(i=0;i<LenRa;i++){
        Ranode[i] = (Ranode[i+1]+Ranode[i])/2;
        std::cout << Ranode[i] << ' ';
    }
    std::cout << std::endl;


    // Cell position of anode and cathode
    // Will need to update when use actual geometry
    Rc = 0;
    Ra = numR;
    Zc = LenRc; 
    Za = LenRa;


    // Need to set up different zones as the bounds for calculation change along the thruster-length
    // first zone will be between the anode and cathode
    // second zone will be between anode and boundary
    // third zone will  be between top and bottom boundary
    // need to identify the corners of each zone - Rmin, Rmax, Zmin, Zmax

    //Corners of ZONE1
    Rmin[1] = Rc;
    Rmax[1] = Ra;
    Zmin[1] = 0;
    if(Zc <= Za){
        Zmax[1] = Zc;
    }     
    else{
        Zmax[1] = Za;
    }
        
    //Corners of ZONE2
    if(Zc <= Za){
        Rmin[2] = 0;
    }
    else{
        Rmin[2] = Rc;
    }
    if(Zc <= Za){
        Rmax[2] = Ra;
    }
    else{
        Rmax[2] = numR;
    }
    if(Zc <= Za){
        Zmin[2] = Zc;
    }
    else{
        Zmin[2] = Za;
    }
    if(Zc <= Za){
        Zmax[2] = Za;
    }
    else{
        Zmax[2] = Zc;
    }

    //Corners of ZONE3
    Rmin[3] = 0;
    Rmax[3] = numR;
    if(Zc <= Za){
        Zmin[3] = Za;
    }
    else{
        Zmin[3] = Zc;
        Zmax[3] = numZ;
    }


    // Printing extra stuff out

    // int LenRa = sizeof(Rcathode)/sizeof(Rcathode[0]);
    // std::cout << LenRc << std::endl;
    // int RcathodeSize = sizeof(Rcathode)/sizeof(Rcathode[0]);
    // int RanodeSize = sizeof(Ranode)/sizeof(Ranode[0]);
    // std::cout << RcathodeSize << std::endl;

    // std::cout << "--------------------------------" << std::endl;
    // std::cout << "printing positions Rcathode" << std::endl;
    // for(i=0;i<=imax;i++){
    //     if(y[0][i] >= 0){
    //         std::cout << Rcathode[i] << ' ';
    //     }
    // }
    // std::cout << std::endl;
    
    // std::cout << "printing positions Ranode" << std::endl;
    // for(i=0;i<=imax;i++){
    //     if(y[jmax][i] < jmax){
    //         std::cout << Ranode[i] << ' ';
    //     }
    // }
    // std::cout << std::endl;


}