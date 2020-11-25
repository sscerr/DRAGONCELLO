#ifndef _DRAGONCELLO_H
#define _DRAGONCELLO_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "constants.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>

using namespace std;
using namespace DRAGON;

struct diagEntry {
  double upper;
  double central;
  double lower;
};

struct norms {
    double max;
    double mean;
    double norm2;
};

class DRAGONCELLO {
  
public:
  
  DRAGONCELLO(int,int);
  ~DRAGONCELLO() { cout << "I am the destructor " << endl; }
  
  void buildXGrid();
  void buildYGrid();
  void buildZGrid();
  void buildPGrid();
  void initEmptyGrids();
  void initGrids();
  void initSource();
  void initBfield();
  void initAnisoSpatialDiffusion();
  void initDriftLikeVelocities();
  void initEloss();
  
  void dumpSpectra();
  void dumpSpectra_2(int counter);
  void dumpProfiles();  
  void dumpProfiles2D();

  void compute_epsilon(int counter);  

  void computeResidualLosses(int);
  void computeResidualDiffAniso(int);

  void computeSourceEnergyIntegral(int, int, int);
  void computeSolutionEnergyIntegral(int, int, int);
   
  void initCN2Daniso();
  
  void Run2D();
    
  void propagateInZaniso(double dt, double dtbar, int counter);
  void propagateInRaniso(double dt, double dtbar, int counter);
  void propagateMixedDerivatives(double dt, double dtbar, int counter);
  void propagateEloss(double dt, double dtbar);
  
  void printLocalFlux();
  void dumpSolution(int counter);

  
  inline long int index(int ix, int iy, int iz, int ip) {
    return  (((ix * ny + iy) * nz + iz) * np + ip);
  }
  
  inline long int indexSpatial(int ix, int iy, int iz) {
    return  (((ix * ny + iy) * nz + iz));
  }
  
  inline double computeVelocity(double momentum) {  
	  return c_light*sqrt( 1. - mass*mass*c_light*c_light/(momentum*momentum)  );
  }
  
  inline vector<double> getP() {
	  return p_vec;
  }
  
  void solveTridiagonalSystem(vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& r, vector<double>& u, int n);
  int gsl_linalg_solve_tridiag(const vector<double> & diag, const vector<double> & abovediag, const vector<double> & belowdiag,
                                const vector<double> & rhs, vector<double> & x);
  int solve_tridiag_nonsym(const vector<double> & diag, const vector<double> & abovediag, const vector<double> & belowdiag,
                            const vector<double> & rhs, vector<double> & x, size_t N);
   
  
private:
  
  int Z, A;
  int species_id;
  
  double deltax, deltay, deltaz, deltaLogp;
  vector<double> x_vec, y_vec, z_vec, p_vec;
  vector<double> deltaxUp, deltayUp, deltazUp;
  vector<double> deltaxDown, deltayDown, deltazDown;
  vector<double> deltaxCentral, deltayCentral, deltazCentral;
  
  unsigned int ixsun, iysun, izsun, ipLocal;
  
  vector<double> N;          //(r,z,p) or (x,y,z,p)
  vector<double> N_previous; //(r,z,p) or (x,y,z,p)
  
  double sourceEnergyIntegral;
  double totalInjectedEnergy;
  double solutionEnergyIntegral;
  double reaccelerationPower;  
  double totalReaccelerationEnergy;
  
  vector<double> analyticalSolution;    
  
  vector<double> sourceTerm;

  vector<double> b_r;  // (r,z) or (r,phi,z); b_r = B_r/|B| magnetic field unit vector projected along r   
  vector<double> b_phi;  // (r,z) or (r,phi,z);
  vector<double> b_x;  // b_x = B_x/|B| magnetic field unit vector projected along x   
  vector<double> b_y;  // b_y = B_y/|B| 
  vector<double> b_z;  // same as before, but along z

  vector<double> diffusionCoefficient;      //(r,z,p) or (x,y,z,p)
  vector<double> diffusionCoefficient_xx;      //(r,z,p) or (x,y,z,p)
  vector<double> diffusionCoefficient_zz;      //(r,z,p) or (x,y,z,p)  

  vector<double> Dpara;  //(r,z,p) or (x,y,z,p)
  vector<double> Dperp;  //(r,z,p) or (x,y,z,p)  

  vector<double> D_rr;  //(r,z,p) or (r,phi,z,p)
  vector<double> D_rz;  //(r,z,p) or (r,phi,z,p)  
  vector<double> D_xx;  //(x,z,p) or (x,y,z,p)
  vector<double> D_xz;  //(x,z,p) or (x,y,z,p)  
  vector<double> D_zz;  //(r,z,p) or (x,z,p) or ...  
  vector<double> u_r;  //(r,z,p) or (x,y,z,p)  
  vector<double> u_x;  //(x,z,p) or (x,y,z,p)  
  vector<double> u_z;  //(r,z,p) or (x,z,p) or ...  
  
  vector<double> reaccelerationCoefficient;      //(r,z,p) or (x,y,z,p)
  vector<double> energyLossTerm;      //(r,z,p) or (x,y,z,p)

  vector<double> upperDiagonalZ;
  vector<double> lowerDiagonalZ;
  vector<double> centralDiagonalZ;
  
  vector<double> upperDiagonalX;
  vector<double> lowerDiagonalX;
  vector<double> centralDiagonalX;

  vector<double> upperDiagonalR;
  vector<double> lowerDiagonalR;
  vector<double> centralDiagonalR;
  
  vector<double> upperDiagonalDpp;
  vector<double> lowerDiagonalDpp;
  vector<double> centralDiagonalDpp;
  
  vector<double> upperDiagonalEloss;
  vector<double> lowerDiagonalEloss;
  vector<double> centralDiagonalEloss;
};

#endif

