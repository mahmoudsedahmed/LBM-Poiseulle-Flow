#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>

using namespace std;

int i;
int j;
int k;
int kj;
int xnew;
int ynew;
int kopp;

//////////////////////////////////////////////////////////////
// Unit Conversion
double l_p=0.01 ;
const int N = 50;
int dx_lbm = 1 ;
double dx_p = l_p * dx_lbm / N ;
double u_p = 0.01;
double nio_p = 1*pow(10,-6);
double Re=l_p*u_p/nio_p;
int dt_lbm = 1 ;
double taw = 0.75;
double omega = 1 / taw;
double c_s = 1 / sqrt(3) ;
double nio_lbm = (pow(c_s,2) * (taw-0.5)) * dt_lbm ;
double dt_p = nio_lbm * (pow(dx_p,2)) / nio_p ;
double u_lbm = u_p * dt_p / dx_p ;
//////////////////////////////////////////////////////////////

const int a=9;
double w[9] ={4./9. , 1./9. , 1./9. , 1./9. , 1./9. , 1./36. , 1./36. , 1./36. , 1./36. };
int cx[9] = {0 , 1 , 0 , -1 , 0 , 1 , -1 , -1 , 1};
int cy[9] = {0 , 0 , 1 , 0 , -1 , 1 , 1 , -1 , -1};
int opp[9] = {0 , 3 , 4 , 1 , 2 , 7 , 8 , 5 , 6};
const int im = (N+1);
const int jm = (7*N);
int solid[im][jm];
double u_x[im][jm];
double u_y[im][jm];
double er[im][jm];
double erabs[im][jm];
double uer[im][jm];
double rho[im][jm];
double rho_o[im][jm];
double cu[im][jm];
double temp;
double fin[im][jm][a];
double fout[im][jm][a];
double feq[im][jm][a];

int main()
{
// Print Grid Dimensions
ofstream imsv("im.txt", ios::out);
imsv << im;
imsv.close();
ofstream jmsv("jm.txt", ios::out);
jmsv << jm;
jmsv.close();

// Initialization
for (i=0; i<im; ++i){
    for (j=0; j<jm; ++j){
        solid[i][j]=0;
        rho[i][j]=0;
        rho_o[i][j]=0;
        u_x[i][j]=0;
        u_y[i][j]=0;
        uer[i][j]=0;
        er[i][j]=0;
        erabs[i][j]=0;
        cu[i][j]=0;
    }
}
// Initialization - fin,fout,feq
for (i=0; i<im; ++i){
     for( j=0; j<jm; ++j){
         for ( k=0; k<9; ++k){
                fin[i][j][k]=0;
                fout[i][j][k]=0;
                feq[i][j][k]=0;
            }
        }
    }
// Geometry Details
for (j=0; j<jm; ++j){
    solid[0][j]=1;
    solid[im-1][j]=1;
}
// Initial Conditions
for (i=0; i<im; ++i){
    for (int j=0; j<jm; ++j){
            if (!solid[i][j]){
                rho[i][j]=1;
                rho_o[i][j]=1;
                for (k=0; k<9; ++k){
                        fin[i][j][k]= w[k] * rho [i][j];
                }
            }
    }
}
//Macroscopic Parameters
for (i=0; i<im; ++i){
    for (j=0; j<jm; ++j){
            rho[i][j]=0;
            u_x[i][j]=0;
            u_y[i][j]=0;
            if (!solid[i][j]){
            for (k=0; k<9; ++k){
                rho[i][j] += fin[i][j][k];
                u_x[i][j] += cx[k] * fin[i][j][k];
                u_y[i][j] += cy[k] * fin[i][j][k];
                }
                u_x[i][j] = u_x[i][j] / rho_o[i][j];
                u_y[i][j] = u_y[i][j] / rho_o[i][j];
            }
    }
}

// Main Calculation Loop
int t;
int tm = 50000;
double conv = pow(10,-8);
int tp = 1000;

for (t=0; t<=tm; ++t){
////////////////////////////////

//////////////////////////////////////////////////////////////
// Boundary Conditions//

for (i=1; i<im-1; ++i){
        // Macroscopic Parameters
        // West Boundary - Inlet Velocity
        u_x[i][0] = u_lbm;
        rho[i][0] = (rho_o[i][0] * u_x[i][0]) + (fin[i][0][0]+fin[i][0][2]+fin[i][0][4]) + (2*(fin[i][0][3]+fin[i][0][6]+fin[i][0][7]));
        u_y[i][0] = 0;
        // East Boundary - Extrapolation & [Zho-He]
        u_x[i][jm-1] = u_x[i][jm-2];
        u_y[i][jm-1] = u_y[i][jm-2];
        rho[i][jm-1] = (-rho_o[i][jm-1] * u_x[i][jm-1]) + (fin[i][jm-1][0]+fin[i][jm-1][2]+fin[i][jm-1][4]) + (2*(fin[i][jm-1][1]+fin[i][jm-1][5]+fin[i][jm-1][8])) ;
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Distribution Functions
        //West Side - Velocity Boundary [Zho - He]
        fin[i][0][1] = fin[i][0][3] + ((2./3.) * rho_o[i][0] * u_x[i][0]);
        fin[i][0][5] = fin[i][0][7] - (0.5 * (fin[i][0][2] - fin[i][0][4])) + ((1./6.) * rho_o[i][0] * u_x[i][0]) + (0.5 * rho_o[i][0] * u_y[i][0]);
        fin[i][0][8] = fin[i][0][6] + (0.5 * (fin[i][0][2] - fin[i][0][4])) + ((1./6.) * rho_o[i][0] * u_x[i][0]) - (0.5 * rho_o[i][0] * u_y[i][0]);
        //East Side - Velocity Boundary [Zho - He]
        fin[i][jm-1][3] = fin[i][jm-1][1] - ((2./3.) * rho_o[i][jm-1] * u_x[i][jm-1]);
        fin[i][jm-1][7] = fin[i][jm-1][5] + (0.5 * (fin[i][jm-1][2] - fin[i][jm-1][4])) - ((1./6.) * rho_o[i][jm-1] * u_x[i][jm-1]) - (0.5 * rho_o[i][jm-1] * u_y[i][jm-1]);
        fin[i][jm-1][6] = fin[i][jm-1][8] - (0.5 * (fin[i][jm-1][2] - fin[i][jm-1][4])) - ((1./6.) * rho_o[i][jm-1] * u_x[i][jm-1]) + (0.5 * rho_o[i][jm-1] * u_y[i][jm-1]);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Macroscopic Parameters
for (i=0; i<im; ++i){
    for (j=0; j<jm; ++j){
            rho[i][j]=0;
            u_x[i][j]=0;
            u_y[i][j]=0;
            if (!solid[i][j]){
            for (k=0; k<9; ++k){
                rho[i][j] += fin[i][j][k];
                u_x[i][j] += cx[k] * fin[i][j][k];
                u_y[i][j] += cy[k] * fin[i][j][k];
                }
                u_x[i][j] = u_x[i][j] / rho_o[i][j];
                u_y[i][j] = u_y[i][j] / rho_o[i][j];
            }
    }
}
// Results Print
if (t%tp==0 && t>1){
// Display Run Status
//cout <<"Iteration = " << t << " " << "Max. Error = " << temp << endl;
// Print Run Status
if (t==tp){ // First Print
ofstream Con("Console.txt", ios::out);
Con << "Iteration = " << t << " " << "Max. Error = " << temp << endl;
Con.close();
}
if (t>tp){  // Next Prints
ofstream Con("Console.txt", ios::app);
Con << "Iteration = " << t << " " << "Max. Error = " << temp << endl;
Con.close();
}

// Save Results
ofstream rhot("rhot.txt", ios::out);
for (int j=0; j<jm; ++j){
    for (int i=0; i<im; ++i){
      rhot << rho[i][j] << endl ;
        }
}
rhot.close();
ofstream u_xt("u_xt.txt", ios::out);
for (int j=0; j<jm; ++j){
    for (int i=0; i<im; ++i){
      u_xt << u_x[i][j] << endl ;
        }
}
u_xt.close();
}

if (t>1){ //Convergence Criteria
temp = 0;
for (i=0; i<im; ++i){
    for (j=0; j<jm; ++j){
        er[i][j] = (u_x[i][j] - uer[i][j]) / u_x[i][j];
        erabs[i][j] = abs(er[i][j]);
            if(erabs[i][j]>temp){
            temp = erabs[i][j];
        }
    }
}
if(temp<conv){  // Solution Converged
        //cout << "Solution Converged" << endl;
        //cout <<"Iteration = " << t << " " << "Max. Error = " << temp << endl;
        ofstream Con("Console.txt", ios::app);
        Con <<"Solution Converged" << endl;
        Con << "Iteration = " << t << " " << "Max. Error = " << temp << endl;
        Con.close();
    break;
}
}

// uer
for (i=0; i<im; ++i){
    for(j=0; j<jm; ++j){
        uer[i][j] = u_x[i][j];
            }
}

//Collision Step
for (i=0; i<im; ++i){
    for (j=0; j<jm; ++j){
     if(!solid[i][j]){
      for (k=0; k<9; ++k) {
        cu[i][j] = 3*(cx[k] * u_x[i][j] + cy[k] * u_y[i][j]);
        feq[i][j][k] = w[k] * (rho[i][j] + (rho_o[i][j] * (cu[i][j] + 0.5 * (cu[i][j] * cu[i][j]) - 1.5 * ((u_x[i][j] * u_x[i][j]) +(u_y[i][j] * u_y[i][j])))));
        fout[i][j][k] = fin[i][j][k] - ( omega * (fin[i][j][k] - feq[i][j][k])) ;
      }
      }
      //On-Grid Bounceback
        if(solid[i][j]){
         for(k=0; k<9; ++k){
         kj = opp[k];
         fout[i][j][k] = fin[i][j][kj];
         }
        }
    }
}

//Streaming Step
for (i=0; i<im; ++i){
    for(j=0; j<jm; ++j){
            xnew = 0;
            ynew = 0;
            for(k=0; k<=8; ++k){
                xnew = j + cx[k];
                ynew = i - cy[k];
        if (xnew >= 0 && xnew <= (jm-1) && ynew >= 0 && ynew <= (im-1) ) {
                                fin[ynew][xnew][k] = fout[i][j][k];
                }
        }
    }
}

////////////////////////////////////////////////////

} // main loop

/////////////////////////////////////////////////////////////////////////////
// Final Results
////////////////////////////////////////////////////////////////////////////
ofstream rhot("rhot.txt", ios::out);
for (int j=0; j<jm; ++j){
    for (int i=0; i<im; ++i){
      rhot << rho[i][j] << endl ;
        }
}
rhot.close();
ofstream u_xt("u_xt.txt", ios::out);
for (int j=0; j<jm; ++j){
    for (int i=0; i<im; ++i){
      u_xt << u_x[i][j] << endl ;
        }
}
u_xt.close();

} //main
