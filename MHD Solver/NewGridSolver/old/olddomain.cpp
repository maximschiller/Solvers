/* SCRIPT WHICH INITIALISES THE DOMAIN AND THEN DISCRETISES IT NON-ORTHOGONALLY USING TFI AND LAPLACE */
#include "Variables.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>

// Declarations
FILE *pt1, *pt2;
// const int xmax=300, ymax = 100;
char filename1[] = "trnsFinGrid.dat";
char filename2[] = "ellipGrid.dat";


/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

// Take input file from software such as GMSH and only use the boundary points so that discretisation can be done in-house
void DomainInput(){


}


/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

// Potentially take input from a software which builds the outline of a surface or somehow streamline the geometry building process
// Could also be a manual way of doing geometry building
void discretiseGeom(){



}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

// Note that in development of the domain the origin has assumed to be top left corner
// Function which generates a grid from boundary
void gridgen(){

    int iter; // iteration counter
    double disxi, diseta, dx, dy, height, t, s; // domain variables + variables for TFI
    double omega, resid, alf, beta, gamdomain, xtemp, ytemp, a, aa, c, cc, xeta, yeta, xxi, yxi,
            J, PP = 0., QQ = 0.; // variables for laplace grid generation

    // Definitions
    // Variables used for grid dimensions in the eta and xi space
    // imax = 6; // number of discretisations in computational domain
    // jmax = 6; // number of discretisations in computational domain
    dx = 1./imax;
    dy = 1./jmax;
    disxi = 1./(imax);
    diseta = 1./(jmax);
    height = 0.5; // height parabola starts at

    /*--------------------------------------------------------------*/

    // Opening write files for TFI and Laplace generated grids
    pt1 = fopen(filename1,"w");
    pt2 = fopen (filename2, "w");


    // Compute boundary points of domain - note that i is the columns (in the xi direction) and j is the rows (in the eta direction)
    // At the moment just using a basic geometry - hopefully can input a data file with discretised points on a geometry
    for (i=0;i<=imax;i++){
        // Bottom Curve - xi = (0,1), eta = 0
        x[i][0] = dx*i;
        y[i][0] = 0;
        // Top Curve
        x[i][jmax] = dx*i;
        y[i][jmax] = 0.5*(pow(dx*i,2)) + 0.5;
    }
    for (j=0;j<=jmax;j++){
        // Left curve
        x[0][j] = 0;
        y[0][j] = (dy/2.)*j;
        // Right curve
        x[imax][j] = 1;
        y[imax][j] = dy*j;
    }

    // for (int i = 0; i <=imax; i++){
    //     for (int j = 0; j <= jmax; j++){
    //         std::cout << x[i][j] << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "--------------------------------" << std::endl;
    // for (int i = 0; i <=imax; i++){
    //     for (int j = 0; j <= jmax; j++){
    //         std::cout << y[i][j] << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "--------------------------------" << std::endl;

    /*--------------------------------------------------------------*/

    // Trans-Finite interpolation for inner points
    for(j=1;j<jmax;j++){
        t = (j-1)*dy;
        for(i=1;i<imax;i++){
            s = (i-1)*dx;
            // ii = double(i);
            // jj = double(j);
            // x[i][j] = (ii/imax)*x[imax][j] + ((imax-ii)/imax)*x[0][j] + (jj/jmax)*x[i][jmax] + ((jmax-jj)/jmax)*x[i][0] - ((ii/imax)*(jj/jmax))*x[imax][jmax]
            //             -((ii/imax)*((jmax-jj)/jmax))*x[imax][0]-(((imax-ii)/imax)*(jj/jmax))*x[0][jmax] - (((imax-ii)/imax)*((jmax-jj)/jmax))*x[0][0];
            // y[i][j] = (ii/imax)*y[imax][j] + ((imax-ii)/imax)*y[0][j] + (jj/jmax)*y[i][jmax] + ((jmax-jj)/jmax)*y[i][0] - ((ii/imax)*(jj/jmax))*y[imax][jmax]
            //             -((ii/imax)*((jmax-jj)/jmax))*y[imax][0] - (((imax-ii)/imax)*(jj/jmax))*y[0][jmax] - (((imax-ii)/imax)*((jmax-jj)/jmax))*y[0][0];

            x[i][j] = (1.-s)*x[0][j] + s*x[imax][j] + (1.-t)*x[i][0] + t*x[i][jmax]
                        -(1.-s)*(1.-t)*x[0][0] - (1.-s)*t*x[0][jmax]
                        -s*(1.-t)*x[imax][0] - s*t*x[imax][jmax];
            y[i][j] = (1.-s)*y[0][j] + s*y[imax][j] + (1.-t)*y[i][0] + t*y[i][jmax]
                        -(1.-s)*(1.-t)*y[0][0] - (1.-s)*t*y[0][jmax]
                        -s*(1.-t)*y[imax][0] - s*t*y[imax][jmax];
        }
    }

    // std::cout << "after TFI" << std::endl;
    // for (int i = 0; i <=imax; i++){
    //     for (int j = 0; j <= jmax; j++){
    //         std::cout << x[i][j] << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "--------------------------------" << std::endl;
    // for (int i = 0; i <=imax; i++){
    //     for (int j = 0; j <= jmax; j++){
    //         std::cout << y[i][j] << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "--------------------------------" << std::endl;

    // Writing trans-finite data file in TECPLOT format
    fprintf(pt1,"VARIABLES = \"x\",\"y\"\n");
    fprintf(pt1,"ZONE T = \"0\" I = %d J = %d\n",imax,jmax);
    for(j=0;j<=jmax;j++){
        for(i=0;i<=imax;i++){
            fprintf(pt1,"%f %f\n",x[i][j],y[i][j]); 
        }
    }

    /*--------------------------------------------------------------*/

    // Solving Laplace Equation using Gauss-Seidel Algorithm
    omega = 0.3;
    resid = 1.;
    alf = 0.;
    beta = 0.;
    gamdomain = 0.;
    xtemp = 0.;
    ytemp = 0.;
    iter = 0;
    a = 0.1;    aa = 225.;
    c = 3.6;    cc = 8.5;

    while(resid>1e-6){ // Convergence requirement
        resid = 0.;
        iter++;
        for(i=1;i<imax;i++){
            for(j=1;j<jmax;j++){
            
                xeta = (x[i][j+1]-x[i][j-1])/(2.*diseta);
                yeta = (y[i][j+1]-y[i][j-1])/(2.*diseta);
                xxi = (x[i+1][j]-x[i-1][j])/(2.*disxi);
                yxi = (y[i+1][j]-y[i-1][j])/(2.*disxi);
                J = xxi*yeta - xeta*yxi;

                alf = (xeta*xeta+yeta*yeta);
                beta = (xxi*xeta+yxi*yeta);
                gamdomain = (xxi*xxi+yxi*yxi); 

                if(fabs((double)i/(imax)-0.5) == 0.){ 
                    PP = 0.;
                }
                else{
                     PP = -a*((double)i/(imax)-0.5)/(fabs((double)i/
                            (imax)-0.5))*exp(-c*fabs((double)i/(imax)-
                            0.5));
                }
                if(fabs((double)j/(jmax)-0.0) == 0.){
                    QQ = 0.;
                }
                else{
                    QQ = -aa*((double)j/(jmax)-0.0)/(fabs((double)j/
                            (jmax)-0.0))*exp(-cc*fabs((double)j/(jmax)-
                            0.0));
                }

                xtemp = pow(disxi*diseta,2.)/(2.*(alf*diseta*diseta+gamdomain*disxi*disxi))*
                            (alf/(disxi*disxi)*(x[i+1][j]+x[i-1][j]) + gamdomain/(diseta*
                            diseta)*(x[i][j+1]+x[i][j-1]) - beta/(2.*disxi*diseta)*
                            (x[i+1][j+1]+x[i-1][j-1]-x[i-1][j+1]-x[i+1][j-1])
                            + (J*J)*(xxi*PP+xeta*QQ));
                ytemp = pow(disxi*diseta,2.)/(2.*(alf*diseta*diseta+gamdomain*disxi*disxi))*
                            (alf/(disxi*disxi)*(y[i+1][j]+y[i-1][j]) + gamdomain/(diseta*
                            diseta)*(y[i][j+1]+y[i][j-1]) - beta/(2.*disxi*diseta)*
                            (y[i+1][j+1]+y[i-1][j-1]-y[i-1][j+1]-y[i+1][j-1])
                            + (J*J)*(yxi*PP+yeta*QQ));

                resid = resid + pow((x[i][j]-xtemp),2.) + pow((y[i][j]-ytemp),2.);

                xtemp = omega*xtemp + (1.-omega)*x[i][j];
                ytemp = omega*ytemp + (1.-omega)*y[i][j];
                x[i][j] = xtemp;
                y[i][j] = ytemp;
            }
        }
        resid = sqrt(resid);
    } 
    
    // Writing data file in TECPLOT format
    fprintf(pt2,"VARIABLES = \"x\",\"y\"\n");
    fprintf(pt2,"ZONE T = \"0\" I = %d J = %d\n",imax,jmax);
    for(j=0;j<=jmax;j++){
        for(i=0;i<=imax;i++){
            fprintf(pt2,"%f %f\n",x[i][j],y[i][j]); 
        }
    }

    std::cout << "after laplace" << std::endl;
    for (int i = 0; i <=imax; i++){
        for (int j = 0; j <= jmax; j++){
            std::cout << x[i][j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << "--------------------------------" << std::endl;
    for (int i = 0; i <=imax; i++){
        for (int j = 0; j <= jmax; j++){
            std::cout << y[i][j] << ' ';
        }
        std::cout << std::endl;
    }

    std::cout << y[3][5] << std::endl;
    std::cout << y[4][5] << std::endl;

    fclose(pt1);
    fclose(pt2);

    /*--------------------------------------------------------------*/
    // Outputing data so that the mesh can be plotted in python

    std::ofstream outfile_x("meshxvals.txt");
    // Output array values to file in a format that can be easily imported into Python
    for(i=0; i<=imax; i++) {
        for(j=0; j<=jmax; j++) {
            outfile_x << x[i][j] << " ";
        }
        outfile_x << std::endl;
    }
    outfile_x.close();

    std::ofstream outfile_y("meshyvals.txt");
    // Output array values to file in a format that can be easily imported into Python
    for(i=0; i<=imax; i++) {
        for(j=0; j<=jmax; j++) {
            outfile_y << y[i][j] << " ";
        }
        outfile_y << std::endl;
    }
    outfile_y.close();

}





