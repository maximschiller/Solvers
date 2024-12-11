/* SCRIPT WHICH INITIALISES THE DOMAIN AND THEN DISCRETISES IT NON-ORTHOGONALLY USING TFI AND LAPLACE */
#include "Variables.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <vector>

// Declarations
FILE *pt1, *pt2;
// const int xmax=300, ymax = 100;
char filename1[] = "trnsFinGrid.dat";
char filename2[] = "ellipGrid.dat";

// Declarations for domain inputting and discretisation
std::string data = "condi"; // variable name of domain geometry text file
double col1, col2;
std::vector<double> xtop;
std::vector<double> ytop;
std::vector<double> xleft;
std::vector<double> yleft;
std::vector<double> xbot;
std::vector<double> ybot;
std::vector<double> xright;
std::vector<double> yright;


/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

// Take input file from software such as GMSH and only use the boundary points so that discretisation can be done in-house
// given a certain spacing, i.e. the number of data points in each boundary can read in a single text file which outputs the corresponding boundary
// i.e. text file could have boundaries in this order, top, left, bottom, right and write function which takes these and splits them into corresponding boundaries

// currently have it so that it reads in seperate text files for each boundary - need to fix, finalise and test this
void DomainInput(){

    // xtop and ytop data input
    std::string filenametop = "geometry/topbound.txt"; // file name and its location
    std::ifstream file(filenametop);
    while (file >> col1 >> col2) {
        xtop.push_back(col1);
        ytop.push_back(col2);
    }

    file.close();

    // xleft and yleft data input
    std::string filenameleft = "geometry/leftbound.txt"; // file name and its location
    std::ifstream file2(filenameleft);
    while (file2 >> col1 >> col2) {
        xleft.push_back(col1);
        yleft.push_back(col2);
    }

    file2.close();

        // xbot and ybot data input
    std::string filenamebot = "geometry/botbound.txt"; // file name and its location
    std::ifstream file3(filenamebot);
    while (file3 >> col1 >> col2) {
        xbot.push_back(col1);
        ybot.push_back(col2);
    }

    file3.close();

    // xright and yright data input
    std::string filenameright = "geometry/rightbound.txt"; // file name and its location
    std::ifstream file4(filenameright);
    while (file4 >> col1 >> col2) {
        xright.push_back(col1);
        yright.push_back(col2);
    }

    file4.close();



    // checking if inputting works
    for (auto& x : xtop) {
    std::cout << x << " ";
    }
    std::cout << std::endl;

    std::cout << "printed values for boundary from python" << std::endl;  


}


/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

// Potentially take input from a software which builds the outline of a surface or somehow streamline the geometry building process
// Could also be a manual way of doing geometry building
// given two endpoints of a boundary/shape of the boundary use this function to discretise along the boundary
// could potentially update so that imax = spacing, would make it a lot easier
void discretiseGeom(){

    // Compute boundary points of domain - note that i is the columns (in the xi direction) and j is the rows (in the eta direction)
    // At the moment just using a basic geometry - hopefully can input a data file with discretised points on a geometry
    // read in data from domain input

    // Basic geom generation
    // for (i=0;i<=imax;i++){
    //     // Bottom Curve - xi = (0,1), eta = 0
    //     x[0][i] = dx*i;
    //     y[0][i] = 0;
    //     // Top Curve
    //     x[jmax][i] = dx*i;
    //     y[jmax][i] = 0.5*(pow(dx*i,2)) + 0.5;
    // }
    // for (j=0;j<=jmax;j++){
    //     // Left curve
    //     x[j][0] = 0;
    //     y[j][0] = (dy/2.)*j;
    //     // Right curve
    //     x[j][imax] = 1;
    //     y[j][imax] = dy*j;
    // }


    /*--------------------------------------------------------------*/


    // Python script geom generation
    // getting the values of the data out up to imax
    // xtop, ytop, xleft, yleft, xbot, ybot, xright, yright
    std::vector<double> spaced_values_xtop;
    std::vector<double> spaced_values_ytop;
    std::vector<double> spaced_values_xleft;
    std::vector<double> spaced_values_yleft;
    std::vector<double> spaced_values_xbot;
    std::vector<double> spaced_values_ybot;
    std::vector<double> spaced_values_xright;
    std::vector<double> spaced_values_yright;



    // note that all the boundaries do not have the same number of points from the python script
    // need to have imax and jmax seperate

    // top
    int total_length_top = xtop.size();
    for (int i = 0; i <= imax-1; i++) {
        float step = (float)(total_length_top - 1) / (float)(imax);
        int index = (int)(i * step + 0.5);

        spaced_values_xtop.push_back(xtop[index]);
        spaced_values_ytop.push_back(ytop[index]);
    }

    spaced_values_xtop.push_back(xtop[xtop.size()-1]); // add the last value to the end
    spaced_values_ytop.push_back(ytop[ytop.size()-1]); // add the last value to the end


    // bottom
    int total_length_bot = xbot.size();
    for (int i = 0; i <= imax-1; i++) {
        float step = (float)(total_length_bot - 1) / (float)(imax);
        int index = (int)(i * step + 0.5);

        spaced_values_xbot.push_back(xbot[index]);
        spaced_values_ybot.push_back(ybot[index]);
    }

    spaced_values_xbot.push_back(xbot[xbot.size()-1]); // add the last value to the end
    spaced_values_ybot.push_back(ybot[ybot.size()-1]); // add the last value to the end

    // left
    int total_length_left = xleft.size();
    for (int j = 0; j <= jmax-1; j++) {
        float step = (float)(total_length_left - 1) / (float)(jmax);
        int index = (int)(j * step + 0.5);

        spaced_values_xleft.push_back(xleft[index]);
        spaced_values_yleft.push_back(yleft[index]);

    }

    spaced_values_xleft.push_back(xleft[xleft.size()-1]); // add the last value to the end
    spaced_values_yleft.push_back(yleft[yleft.size()-1]); // add the last value to the end


    // right
    int total_length_right = xright.size();
    for (int j = 0; j <= jmax-1; j++) {
        float step = (float)(total_length_right - 1) / (float)(jmax);
        int index = (int)(j * step + 0.5);

        spaced_values_xright.push_back(xright[index]);
        spaced_values_yright.push_back(yright[index]);
    }

    spaced_values_xright.push_back(xright[xright.size()-1]); // add the last value to the end
    spaced_values_yright.push_back(yright[yright.size()-1]); // add the last value to the end

    // print out the result
    for (int i = 0; i < spaced_values_yleft.size(); i++) {
        std::cout << spaced_values_yleft[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "printed values for discretised array" << std::endl;


    /*--------------------------------------------------------------*/

    
    // assigning the value to x and y
    // x vals
    // bottom 
    for (int i = 0; i <= spaced_values_xbot.size(); i++) {
        x[0][i] = spaced_values_xbot[i];
        y[0][i] = spaced_values_ybot[i];
    }

    // top
    for (int i = 0; i <= spaced_values_xtop.size(); i++) {
        x[jmax][i] = spaced_values_xtop[i];
        y[jmax][i] = spaced_values_ytop[i];
    }


    // left
    for (int j = 1; j <= spaced_values_xleft.size(); j++) {
        x[j][0] = spaced_values_xleft[j];
        y[j][0] = spaced_values_yleft[j];
    }

    // right
    for (int j = 1; j <= spaced_values_xright.size(); j++) {
        x[j][imax] = spaced_values_xright[j];
        y[j][imax] = spaced_values_yright[j];
    }
    


    /*--------------------------------------------------------------*/

    // printing the values to the screen
    for (int j = 0; j <= jmax; j++){    
        for (int i = 0; i <=imax; i++){
            std::cout << x[j][i] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << "--------------------------------" << std::endl;
    for (int j = 0; j <= jmax; j++){    
        for (int i = 0; i <=imax; i++){
            std::cout << y[j][i] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << "--------------------------------" << std::endl;



    // normalising the x and y matrices for grid generation
    // can either normalise the arrays or just set Z and RAD to be equivalent to the domain size
    // i think normalising is easier and gives an independent study of domain size
    // x matrix
    double min_val = 0, max_val = 0;
    double range = 0;
    for (int j = 0; j <= jmax; j++){    
        for (int i = 0; i <=imax; i++){
            if (x[j][i] < min_val) {
                min_val = x[j][i];
            }
            if (x[j][i] > max_val) {
                max_val = x[j][i];
            }
        }
    }

    range = max_val - min_val;
    for (int j = 0; j <= jmax; j++){    
        for (int i = 0; i <=imax; i++){
            x[j][i] = (x[j][i] - min_val) / range;
        }
    }

    // y matrix
    for (int j = 0; j <= jmax; j++){    
        for (int i = 0; i <=imax; i++){
            if (y[j][i] < min_val) {
                min_val = y[j][i];
            }
            if (x[j][i] > max_val) {
                max_val = y[j][i];
            }
        }
    }

    range = max_val - min_val;
    for (int j = 0; j <= jmax; j++){    
        for (int i = 0; i <=imax; i++){
            y[j][i] = (y[j][i] - min_val) / range;
        }
    }

}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

// Note that in development of the domain the origin has assumed to be top left corner
// Function which generates a grid from boundary
void gridgen(){

    int iter; // iteration counter
    double disxi, diseta, height, t, s; // domain variables + variables for TFI
    double omega, resid, alf, beta, gamdomain, xtemp, ytemp, a, aa, c, cc, xeta, yeta, xxi, yxi,
            J, PP = 0., QQ = 0.; // variables for laplace grid generation

    // Definitions
    // Variables used for grid dimensions in the eta and xi space
    // imax = 6; // number of discretisations in computational domain
    // jmax = 6; // number of discretisations in computational domain
    // dx = deltaZ;
    // dy = 1./jmax;
    // disxi = 1./(imax);
    // diseta = 1./(jmax);

    disxi = deltaZ;
    diseta = deltaR;
    height = 0.5; // height parabola starts at

    /*--------------------------------------------------------------*/

    // Opening write files for TFI and Laplace generated grids
    pt1 = fopen(filename1,"w");
    pt2 = fopen (filename2, "w");

    /*--------------------------------------------------------------*/

    // Trans-Finite interpolation for inner points
    for(i=1;i<imax;i++){
        s = (i-1)*dx;
        for(j=1;j<jmax;j++){
            t = (j-1)*dy;
    
            x[j][i] = (1.-s)*x[j][0] + s*x[j][imax] + (1.-t)*x[0][i] + t*x[jmax][i]
                        -(1.-s)*(1.-t)*x[0][0] - (1.-s)*t*x[jmax][0]
                        -s*(1.-t)*x[0][imax] - s*t*x[jmax][imax];
            y[j][i] = (1.-s)*y[j][0] + s*y[j][imax] + (1.-t)*y[0][i] + t*y[jmax][i]
                        -(1.-s)*(1.-t)*y[0][0] - (1.-s)*t*y[jmax][0]
                        -s*(1.-t)*y[0][imax] - s*t*y[jmax][imax];
        }
    }

    std::cout << "after TFI" << std::endl;
    for (int j = 0; j <= jmax; j++){    
        for (int i = 0; i <=imax; i++){
            std::cout << x[j][i] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << "--------------------------------" << std::endl;
    for (int j = 0; j <= jmax; j++){    
        for (int i = 0; i <=imax; i++){
            std::cout << y[j][i] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << "--------------------------------" << std::endl;

    // Writing trans-finite data file in TECPLOT format
    fprintf(pt1,"VARIABLES = \"x\",\"y\"\n");
    fprintf(pt1,"ZONE T = \"0\" I = %d J = %d\n",imax,jmax);
    for(i=0;i<=imax;i++){
        for(j=0;j<=jmax;j++){  
            fprintf(pt1,"%f %f\n",x[j][i],y[j][i]); 
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
        for(j=1;j<jmax;j++){
            for(i=1;i<imax;i++){
            
                xeta = (x[j+1][i]-x[j-1][i])/(2.*diseta);
                yeta = (y[j+1][i]-y[j-1][i])/(2.*diseta);
                xxi = (x[j][i+1]-x[j][i-1])/(2.*disxi);
                yxi = (y[j][i+1]-y[j][i-1])/(2.*disxi);
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
                            (alf/(disxi*disxi)*(x[j][i+1]+x[j][i-1]) + gamdomain/(diseta*
                            diseta)*(x[j+1][i]+x[j-1][i]) - beta/(2.*disxi*diseta)*
                            (x[j+1][i+1]+x[j-1][i-1]-x[j+1][i-1]-x[j-1][i+1])
                            + (J*J)*(xxi*PP+xeta*QQ));
                ytemp = pow(disxi*diseta,2.)/(2.*(alf*diseta*diseta+gamdomain*disxi*disxi))*
                            (alf/(disxi*disxi)*(y[j][i+1]+y[j][i-1]) + gamdomain/(diseta*
                            diseta)*(y[j+1][i]+y[j-1][i]) - beta/(2.*disxi*diseta)*
                            (y[j+1][i+1]+y[j-1][i-1]-y[j+1][i-1]-y[j-1][i+1])
                            + (J*J)*(yxi*PP+yeta*QQ));

                resid = resid + pow((x[j][i]-xtemp),2.) + pow((y[j][i]-ytemp),2.);

                xtemp = omega*xtemp + (1.-omega)*x[j][i];
                ytemp = omega*ytemp + (1.-omega)*y[j][i];
                x[j][i] = xtemp;
                y[j][i] = ytemp;
            }
        }
        resid = sqrt(resid);
    } 
    
    // Writing data file in TECPLOT format
    fprintf(pt2,"VARIABLES = \"x\",\"y\"\n");
    fprintf(pt2,"ZONE T = \"0\" I = %d J = %d\n",imax,jmax);
    for(i=0;i<=imax;i++){
        for(j=0;j<=jmax;j++){
            fprintf(pt2,"%f %f\n",x[j][i],y[j][i]); 
        }
    }

    std::cout << "after laplace" << std::endl;
    for (int j = 0; j <= jmax; j++){    
        for (int i = 0; i <=imax; i++){
            std::cout << x[j][i] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << "--------------------------------" << std::endl;
    for (int j = 0; j <= jmax; j++){    
        for (int i = 0; i <=imax; i++){
            std::cout << y[j][i] << ' ';
        }
        std::cout << std::endl;
    }

    // checking to see if numbers are in the correct order in the array
    // std::cout << y[2][2] << std::endl;
    // std::cout << y[3][2] << std::endl;

    fclose(pt1);
    fclose(pt2);

    /*--------------------------------------------------------------*/
    // Outputing data so that the mesh can be plotted in python

    std::ofstream outfile_x("meshxvals.txt");
    // Output array values to file in a format that can be easily imported into Python
    for(i=0; i<=imax; i++) {
        for(j=0; j<=jmax; j++) {
            outfile_x << x[j][i] << " ";
        }
        outfile_x << std::endl;
    }
    outfile_x.close();

    std::ofstream outfile_y("meshyvals.txt");
    // Output array values to file in a format that can be easily imported into Python
    for(i=0; i<=imax; i++) {
        for(j=0; j<=jmax; j++) {
            outfile_y << y[j][i] << " ";
        }
        outfile_y << std::endl;
    }
    outfile_y.close();

}





