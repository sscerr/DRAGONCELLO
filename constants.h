#ifndef _CONSTANTS_H
#define _CONSTANTS_H

#include <cmath>
#include <string>

#include "cgs.h"



// *** SOURCE SELECTION ***

//#define GREEN_FUNCTION  // point-like source term (Gaussian profile(s))
#define X_0 3.*kpc
#define Z_0 0.
#define SIGMA_0 0.25*kpc
#define X_1 3.*kpc
#define Z_1 0.
#define SIGMA_1 0.25*kpc

//#define ANALYTIC_SOLUTION_TEST // see Kissmann, APP 55:37 (2014) 

#define REALISTIC_SOURCE   // realistic source distribution at and around the disk



// *** B-field CONFIGURATION SELECTION ***

//#define USER_DEFINED_BFIELD
#define FARRAR_SIMPLE_2D

#define PARALLEL

using namespace DRAGON;

#define pow2(A) ((A)*(A))
#define pow3(A) ((A)*(A)*(A))
#define pow4(A) ((A)*(A)*(A)*(A))

// CONSTANTS

const int NUM_THREADS = 4;

const double r_Sun   = 8.3 * kpc;  // value of r at Sun's position
const double x_Sun   = 8.3 * kpc;  // value of x at Sun's position
const double y_Sun   = 0.0 * kpc;  // value of y at Sun's position
const double z_Sun   = 0.0 * kpc;  // value of z at Sun's position

const double aTinyNumber       = 0.01 * kpc;

const double p       = 3.0;     // Number of differential operators to be updated

//PARTICLE
const int Zparticle = -1; //atomic number
const int Aparticle = 0;  //mass number
const double mass = mass_electron;
 
// ALGORITHM
const double dtmin = 0.0001 * Myr;
const double dtmax = 0.005 * Myr; 
const double dtfactor = 1.; 
const int Nrept = 100;
const int interruptAfter = 10;

// GRID
const int nx =  151;
const int ny =  1;
const int nz =  101;
const int np =  15;
const double xmin = 0. * kpc;
const double xmax = 15. * kpc;
const double ymin = 0. * kpc;
const double ymax = 0. * kpc;
const double zmin = -5. * kpc;
const double zmax = 5. * kpc;
const double pmin =    10. * GeV;
const double pmax =  100.0 * GeV;

const double pFactor = exp(log(pmax/pmin)/(np-1));

//ORDERED MAGNETIC FIELD
const double B_0 = 10. * microgauss;
const double B_inf = 0.1 * nanogauss;
const double Bscale_r = 3. * kpc;

//SOURCE
const double z_disk = 0.5 * kpc;
const double epsilon = 0.1;
const double SNR = 1. / 40. / year;
const double ESN = 1e51 * erg;
const double pminSource = 0.1*GeV;
const double pmaxSource = pmax;
const double p_ref = 1.* GeV;
const double injSlope = -2.3;
const double timestepSourceIsActiveUntil = 1.e10;

// DIFFUSION
const double D0 = 1.e28 * cm * cm / s;
const double Dzz = D0;
const double delta = 0.5;
const double reference_rigidity = 1. * GeV; //GV
const double D0para = 1.e28 * cm * cm / s;
const double D0perp = 0.1*D0para; 
const double delta_para = 0.3; 
const double delta_perp = 0.5; 

//REACCELERATION
const double vA = 50. * km / s; 
const double Dpp0 = 4./(3.*delta*(4.-delta*delta)*(4.-delta)) * (vA*vA / D0) * p_ref * p_ref;

//ENERGY LOSSES
const double gasDensity        = 1. / cm3;
const double magneticField     = 1. * muG;
const double hadronicElossNorm =  3.e-7   * GeV/s; //erg/s
const double leptonicElossNorm =  1.e-15  * GeV/s; //erg/s
const double leptonicElossConstant = 1.e-17 * GeV/s;
const double z_losses = 1.0 * kpc;
    
//NORMALIZATION
//const double sourceNorm = 1.06e16
const double sourceNorm =  1.e-17; //cm^-2 s^-1 erg^-1
// L^2/t/L^2 [N] = [Q] 
// [Q] = Q0 1/L {Q = Q0 delta(z)}
// [Q0] = 1/L^2/t/E
// [N] = 1/L^3/E

// N = Q H^2/D

const double localFlux = 5.5e-6; //cm^-2 s^-1 sr^-1 GV^-1
const double localFluxMomentum = 100. * GeV;

#endif
