#include "dragoncello.h"

#include <string>
#include <cmath>

using namespace std;
using namespace DRAGON;

DRAGONCELLO::DRAGONCELLO(int Z_, int A_) : Z(Z_), A(A_), species_id(1000 * Z + A) {
  cout << "I am the constructor " << endl;
}

void DRAGONCELLO::initGrids(void) {

  cout << "Creating grids..." << endl;
  
  buildXGrid();
  buildYGrid(); //code is 3D-ready for next release (currently: ny = 1)
  buildZGrid();
  buildPGrid();
  initEmptyGrids();

}


void DRAGONCELLO::initBfield(void) {
  cout << "Initializing ordered B field..." << endl;

#ifdef USER_DEFINED_BFIELD
  //
  // The user can implement any desired B-field here (i.e., different from the one used in Cerri et al., JCAP 10:019 (2017) //
  //
  cout << " *** user-defined B-field ***" << endl;

  for (int ix = 0; ix < nx; ++ix)
          for (int iy = 0; iy < ny; ++iy)
                  for (int iz = 0; iz < nz; ++iz) {

                           long int ind_sp = indexSpatial(ix,iy,iz);

                           // define B field (example)
                           b_r[ind_sp] = B_0; 
                           b_phi[ind_sp] = 0.;  
                           b_z[ind_sp] = B_0; 

                           // compute local unit vector b = B/|B|
                           double modB = sqrt( b_r[ind_sp]*b_r[ind_sp] + b_phi[ind_sp]*b_phi[ind_sp] + b_z[ind_sp]*b_z[ind_sp] );
                           b_r[ind_sp] /= modB;
                           b_phi[ind_sp] /= modB;
                           b_z[ind_sp] /= modB;
                  }
#endif

#ifdef FARRAR_SIMPLE_2D
  //
  // X-shaped B-field model: Jansson & Farrar, ApJ 757:14 (2012)
  // B-field ring model: Pshirkov et al., ApJ 738:192 (2011) 
  //
  cout << " *** Simplified Farrar-like B-field (axisymmetric: X-shaped halo + rings) ***" << endl;

  const double rcX = 4.8*kpc;
  const double rcD = 5.*kpc;
  const double rX = 3.*kpc;
  const double R0 = 10.*kpc;
  const double R0H = 8.*kpc;
  const double Z0 = 1.*kpc;
  const double Z0H = 1.3*kpc;
  const double ThetaX0 = 49.*M_PI/180.;
  const double BX0 = 4.6; 
  const double BD0 = 4.; 
  const double BH0 = 1.; 

  double Bxx, Bdisk, Bhalo;
  double ThetaX, rp;

  for (int ix = 0; ix < nx; ++ix)
	  for (int iy = 0; iy < ny; ++iy)
		  for (int iz = 0; iz < nz; ++iz) {
			  
			  long int ind_sp = indexSpatial(ix,iy,iz);
			  
			  double r = x_vec[ix];
			  double z = z_vec[iz];
			  
			  //Bxx
			  if (r > rcX) {
				  rp = r-fabs(z)/tan(ThetaX0);
				  ThetaX = ThetaX0;
				  Bxx = BX0*exp(-rp/rX)*(rp/max(aTinyNumber,r));
			  }
			  else {
				  rp = r*rcX/(rcX+fabs(z)/tan(ThetaX0));
				  ThetaX = atan( (fabs(z)/max(0.1*aTinyNumber,r)) + (rcX/max(0.1*aTinyNumber,r))*tan(ThetaX0) );
				  Bxx = BX0*exp(-rp/rX)*pow(rp/max(aTinyNumber,r),2);
			  }
			  
			  //Bdisk
			  Bdisk = BD0*exp(-fabs(z)/Z0);
			  if (r > rcD) 
				  Bdisk *= exp(-(r-r_Sun)/R0);
			  
			  //Bhalo
			  if (fabs(z) < Z0H) { 
				  double Z1H = 0.2*kpc;
				  Bhalo = (BH0/(1.+pow((fabs(z)-Z0H)/Z1H,2)))*(r/R0H)*exp(1.-r/R0H);
			  }
			  else {
				  double Z1H = 0.4*kpc;
				  Bhalo = (BH0/(1.+pow((fabs(z)-Z0H)/Z1H,2)))*(r/R0H)*exp(1.-r/R0H);
			  }
			  
			  // B_r should change sign at z < 0 in order to have "X shape"
			  if (z<0)
			  {
			  	  b_r[ind_sp] =  - Bxx*cos(ThetaX);				  
			  }
		  	  else
			  {
			  	  b_r[ind_sp] =  Bxx*cos(ThetaX);				  
			  }
			  b_z[ind_sp] = Bxx*sin(ThetaX);
			  b_phi[ind_sp] = Bhalo + Bdisk;

			  b_r[ind_sp]   *= microgauss;
			  b_phi[ind_sp] *= microgauss;
			  b_z[ind_sp]   *= microgauss;

			  // this ensure that |B| doesn't go to zero (setting lower limit to 0.01*B_inf ~ 0.001*nanogauss)
			  double modB = sqrt( b_r[ind_sp]*b_r[ind_sp] + b_phi[ind_sp]*b_phi[ind_sp] + b_z[ind_sp]*b_z[ind_sp] ) + 0.01*B_inf;
			  
			  if (z_vec[iz] == 0) {
				  cout << x_vec[ix]/kpc << "\t" << y_vec[iy]/kpc << "\t" << z_vec[iz]/kpc << "\t" << " br =  " << b_r[ind_sp] << " bz =  " << b_z[ind_sp] << endl;
				  cout << rp << "\t" << ThetaX << "\t" << Bxx << endl;
			  }
			  
			  if (modB == 0)
				  cout << "WARNING" << endl;
			  
			  b_r[ind_sp] /= modB;
			  b_phi[ind_sp] /= modB;
			  b_z[ind_sp] /= modB;
		  }
#endif
}


void DRAGONCELLO::initSource(void) {
   cout << "Initializing source term..." << endl;
  
   double sourceProfile = 0.;
   double injSpectrum = 0.;
  
      
   for (int ix = 0; ix < nx; ++ix)
      for (int iy = 0; iy < ny; ++iy)
	 for (int iz = 0; iz < nz; ++iz) {
	    for (int ip = 0; ip < np; ++ip) {

	       sourceProfile = 1.;
	       double sigma = 0.3*kpc;
				  
#ifdef REALISTIC_SOURCE 
               //
               // realistic source distribution: Lorimer et al., MNRAS 372:777 (2006)
               //
	       sourceProfile = exp( -(z_vec[iz]*z_vec[iz]) / (sigma*sigma) );
               sourceProfile *= pow( x_vec[ix]/r_Sun , 1.9 )*exp( -5.00*(x_vec[ix]-r_Sun)/r_Sun - fabs(z_vec[iz])/sigma );
#endif				  


#ifdef GREEN_FUNCTION
               //
               // Green-funtion test: point-like source located at (R,z) = (X_0,Z_0) and width ~ SIGMA_0
               // [ see Appendix B.1 in Cerri et al., JCAP 10:019 (2017) ]
               //
               double alpha0 = M_PI/4.;
               sourceProfile = ( 1./(SIGMA_0*SIGMA_0*2.* M_PI) ) * exp( -( pow( (x_vec[ix] - X_0), 2.) + pow( (z_vec[iz] - Z_0), 2.) ) / ( 2.*SIGMA_0*SIGMA_0 )  );
#endif				  

#ifdef ANALYTIC_SOLUTION_TEST 
               //
               // analytical-solution test described in Kissmann, APP 55:37 (2014)
               // [ see Appendix B.2 in Cerri et al., JCAP 10:019 (2017) ]
               //
               long int ind = index(ix,iy,iz,ip); 

               double kR = 0.5*M_PI/(xmax-xmin);
               double kH = M_PI/(zmax-zmin);
               double psi = cos(kR*x_vec[ix])*cos(kH*z_vec[iz]);
               double dpsidr = - kR*sin(kR*x_vec[ix])*cos(kH*z_vec[iz]);
               double dpsidz = - kH*cos(kR*x_vec[ix])*sin(kH*z_vec[iz]); 
               double d2psidrdz = kR*kH*sin(kR*x_vec[ix])*sin(kH*z_vec[iz]); 

	       sourceProfile = ( kR*kR*D_rr[ind] + kH*kH*D_zz[ind] )*psi - 2.*D_rz[ind]*d2psidrdz  
                               - u_r[ind]*dpsidr - u_z[ind]*dpsidz ;
#endif		
				  
	       injSpectrum = pow( p_vec[ip]/p_ref, injSlope);
	       sourceTerm[index(ix,iy,iz,ip)] = injSpectrum * sourceProfile;
				  
	    }
	}  
}



void DRAGONCELLO::initAnisoSpatialDiffusion(void) {
  cout << "Initializing *ANISOTROPIC* diffusion term..." << endl;
  
  double slope_para = 0;
  double slope_perp = 0;
  
  for (int ix = 0; ix < nx; ++ix)
    for (int iy = 0; iy < ny; ++iy)
      for (int iz = 0; iz < nz; ++iz) {

        long int ind_sp = indexSpatial(ix,iy,iz); 

        for (int ip = 0; ip < np; ++ip) {

          long int ind_tot = index(ix,iy,iz,ip); 

          slope_para = pow((p_vec[ip] / reference_rigidity), delta_para);
          slope_perp = pow((p_vec[ip] / reference_rigidity), delta_perp);

          Dpara[ind_tot] = D0para*slope_para; 
          Dperp[ind_tot] = D0perp*slope_perp; 

          D_rr[ind_tot] = Dperp[ind_tot] + (Dpara[ind_tot]-Dperp[ind_tot])*b_r[ind_sp]*b_r[ind_sp]; 
          D_rz[ind_tot] = (Dpara[ind_tot]-Dperp[ind_tot])*b_r[ind_sp]*b_z[ind_sp]; 
          D_zz[ind_tot] = Dperp[ind_tot] + (Dpara[ind_tot]-Dperp[ind_tot])*b_z[ind_sp]*b_z[ind_sp]; 

          // check: diffusion timescales should be larger than time step!
          if (ix==nx/2 && iz==nz/2 && ip%1==0) {
        	  cout << "*PARALLEL* diffusion timescale at p = " <<  p_vec[ip]/GeV << ": " << zmax*zmax/Dpara[index(ix,iy,iz,ip)]/Myr << " Myr " << endl; 
        	  cout << "*PERPENDICULAR* diffusion timescale at p = " <<  p_vec[ip]/GeV << ": " << xmax*xmax/Dperp[index(ix,iy,iz,ip)]/Myr << " Myr " << endl;
        	  cout << "*PARALLEL* diffusion timescale across a bin at p = " <<  p_vec[ip]/GeV << ": " << deltaz*deltaz/Dpara[index(ix,iy,iz,ip)]/Myr << " Myr " << endl; 
        	  cout << "*PERPENDICULAR* diffusion timescale across a bin at p = " <<  p_vec[ip]/GeV << ": " << deltax*deltax/Dperp[index(ix,iy,iz,ip)]/Myr << " Myr " << endl;
          }
        }
      }
}


void DRAGONCELLO::initDriftLikeVelocities(void) {
  //
  // see Appendix A in Cerri et al., JCAP 10:019 (2017)
  // 
  cout << "Initializing drift-like velocities..." << endl;
    
  for (int ix = 0; ix < nx; ++ix)
    for (int iy = 0; iy < ny; ++iy)
      for (int iz = 0; iz < nz; ++iz)
        for (int ip = 0; ip < np; ++ip) {
          long int ind = index(ix,iy,iz,ip); 
          long int ind_rup, ind_rdown, ind_zup, ind_zdown; 

          ind_rup = (ix < nx-1) ? index(ix+1,iy,iz,ip) : ind;
          ind_zup = (iz < nz-1) ? index(ix,iy,iz+1,ip) : ind;
          ind_rdown = (ix > 0) ? index(ix-1,iy,iz,ip) : ind;
          ind_zdown = (iz > 0) ? index(ix,iy,iz-1,ip) : ind;

          u_r[ind] = D_rr[ind]/max(aTinyNumber,x_vec[ix]) 
                     + (D_rr[ind_rup]-D_rr[ind_rdown])/(2.*deltaxCentral[ix])
                     + (D_rz[ind_zup]-D_rz[ind_zdown])/(2.*deltazCentral[iz]);            

          u_z[ind] = D_rz[ind]/max(aTinyNumber,x_vec[ix]) 
                     + (D_rz[ind_rup]-D_rz[ind_rdown])/(2.*deltaxCentral[ix])
                     + (D_zz[ind_zup]-D_zz[ind_zdown])/(2.*deltazCentral[iz]);            

        }
}



void DRAGONCELLO::initEloss(void) {
  if (A > 0)
    cout << "Calculating Hadronic energy losses..." << endl;
  else
    cout << "Calculating Leptonic energy losses..." << endl;
  
  double value = 0.;
  
  for (int ix = 0; ix < nx; ++ix)
    for (int iy = 0; iy < ny; ++iy)
      for (int iz = 0; iz < nz; ++iz)
        for (int ip = 0; ip < np; ++ip) {
        	
          double profile = 1.;
          
#ifdef TEST_REALISTIC
          profile = 1. * exp(-(z_vec[iz]*z_vec[iz]/(2.*z_losses*z_losses)));          
#endif          
          
          if (A>0) {
        	  value = hadronicElossNorm * gasDensity;
        	  energyLossTerm[index(ix,iy,iz,ip)] = 0.;
          } else {
              value = leptonicElossConstant*profile + leptonicElossNorm * pow2(p_vec[ip]/GeV) * pow2(magneticField*profile/muG); //erg/s
              energyLossTerm[index(ix,iy,iz,ip)] = value; 
          }
          
          if (ix==nx/2 && iz==nz/2 && ip%10==0)
        	  cout << "eloss timescale at p = " <<  p_vec[ip]/GeV << ": " << p_vec[ip]/energyLossTerm[index(ix,iy,iz,ip)]/Myr << " Myr " << endl; 
          
        }
}


void DRAGONCELLO::dumpSpectra(void) {
  cout << "Dumping spectra ..." << endl;
  
  ofstream outfile;
  outfile.open("output/galaxy_spectra.txt");
  outfile << "#p [GeV]" << "   " << "N       " << "   " << "Source  " << "   " << "Dxx     " << "   " << "Dpp     " << endl;
  outfile << scientific << setprecision(2);
  for (int ip = 0; ip < np; ++ip) {
    outfile << p_vec[ip] / GeV << "   ";
    outfile << N[index(ixsun,0,izsun,ip)] << "   ";
    outfile << sourceTerm[index(ixsun,0,izsun,ip)] << "   ";
    outfile << diffusionCoefficient[index(ixsun,0,izsun,ip)] << "   ";
    outfile << reaccelerationCoefficient[index(ixsun,0,izsun,ip)] << "   ";
    outfile << energyLossTerm[index(ixsun,0,izsun,ip)] << endl;
  }
  outfile.close();
}


void DRAGONCELLO::dumpProfiles(void) {
  cout << "Dumping profiles ..." << endl;
  
  ofstream outfile;
  outfile.open("output/galaxy_profiles.txt");
  outfile << "#z [kpc] " << "   " << "N       " << "   " << "Source  " << "   " << "Dxx     " << "   " << "Dpp     " << endl;
  outfile << scientific << setprecision(2);
  for (int iz = 0; iz < nz; ++iz) {
    outfile << z_vec[iz] / kpc << "   ";
    outfile << N[index(ixsun,0,iz,10)] << "   ";
    outfile << sourceTerm[index(ixsun,0,iz,10)] << "   ";
    outfile << diffusionCoefficient[index(ixsun,0,iz,10)] << "   ";
    outfile << reaccelerationCoefficient[index(ixsun,0,iz,10)] << "   ";
    outfile << energyLossTerm[index(ixsun,0,iz,10)] << endl;
  }
  outfile.close();
}


void DRAGONCELLO::dumpProfiles2D(void) {
  cout << "Dumping r-z maps ..." << endl;
  
  ofstream outfile;
  outfile.open("output/2DAnisoMaps.txt");
  outfile << "#r" << "\t" << "z" << "\t" << "N" << "\n"; // "   " << "Source  " << "   " << "Dxx     " << "   " << "Dpp     " << endl;
  outfile << scientific << setprecision(2);
  for (int ix = 0; ix < nx; ++ix) 
	  for (int iz = 0; iz < nz; ++iz) {
		outfile << x_vec[ix] / kpc << "   ";
		outfile << z_vec[iz] / kpc << "   ";
		outfile << N[index(ix,0,iz,0)] * cm3 * GeV << "   ";
		outfile << sourceTerm[index(ix,0,iz,0)] << "   ";
		outfile << diffusionCoefficient[index(ix,0,iz,0)] << "   ";
		outfile << reaccelerationCoefficient[index(ix,0,iz,0)] << "   ";
		outfile << energyLossTerm[index(ix,0,iz,0)] << endl;
	  }
  outfile.close();
}


void DRAGONCELLO::initCN2Daniso(void) {
  //
  // see Appendix A in Cerri et al., JCAP 10:019 (2017)
  // 
  cout << "Calculating the CN coefficients for *ANISOTROPIC* case..." << endl;

  int iy = 0;

  for (int ix = 0; ix < nx; ++ix)
    for (int iz = 0; iz < nz; ++iz)
      for (int ip = 0; ip < np; ++ip) {

        long int ind = index(ix,iy,iz,ip);

        upperDiagonalR.push_back( D_rr[ind] / (deltaxUp[ix]*deltaxCentral[ix])
                                  + u_r[ind]/(2.*deltaxCentral[ix]) );

        centralDiagonalR.push_back( 2. * D_rr[ind] / (deltaxUp[ix]*deltaxDown[ix]) );

        lowerDiagonalR.push_back( D_rr[ind] / (deltaxDown[ix]*deltaxCentral[ix]) 
                                  - u_r[ind]/(2.*deltaxCentral[ix]) );
						
        upperDiagonalZ.push_back( D_zz[ind] / (deltazUp[iz]*deltazCentral[iz])
                                  + u_z[ind]/(2.*deltazCentral[iz]) );

        centralDiagonalZ.push_back( 2. * D_zz[ind] / (deltazUp[iz]*deltazDown[iz]) );

        lowerDiagonalZ.push_back( D_zz[ind] / (deltazDown[iz]*deltazCentral[iz]) 
                                  - u_z[ind]/(2.*deltazCentral[iz]) );
	                
	}
}






