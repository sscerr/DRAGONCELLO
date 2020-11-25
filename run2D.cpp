#include "dragoncello.h"
#include "omp.h"


void DRAGONCELLO::propagateInZaniso(double dt, double dtbar, int counter) {
	
   int iz, ix, ip;
   long int ind;
   double value;

#ifdef PARALLEL 
#pragma omp parallel default(shared) private(ix,iz,ip,ind,value) num_threads(NUM_THREADS) 
#endif         
  {

   vector<double> knownTerm(nz - 2, 0.0);
   vector<double> diagonal(nz - 2, 0.0);
   vector<double> upperDiagonal(nz - 3, 0.0);
   vector<double> lowerDiagonal(nz - 3, 0.0);
   vector<double> solution(nz - 2, 0.0);

#ifdef PARALLEL
#pragma omp for schedule(dynamic)
#endif    	    			
   for (ip = 0; ip < np; ip++) {
      for (ix = 0; ix < nx - 1; ix++) { //notice the boundary condition
         for (iz = 1; iz < nz - 1; iz++) {

            ind = index(ix,0,iz,ip);

            diagonal.at(iz - 1) = 1. + dt/2. * centralDiagonalZ.at(ind);

            if (iz != nz - 2) {
               upperDiagonal.at(iz - 1) = -dt/2. * upperDiagonalZ.at(ind);
	    }
  	    if (iz != 1) {
	       lowerDiagonal.at(iz - 2) = -dt/2. * lowerDiagonalZ.at(ind);
	    }
 	    knownTerm.at(iz - 1) = N.at(ind) * (2. - diagonal.at(iz -1));
	    if (counter <= timestepSourceIsActiveUntil) {
	       knownTerm.at(iz - 1) +=  dtbar * sourceTerm.at(ind);
	    }
	    if (iz != nz - 2)
	       knownTerm.at(iz - 1) -= N[index(ix,0,iz+1,ip)] * upperDiagonal.at(iz - 1);
	    if (iz != 1)
	       knownTerm.at(iz - 1) -= N[index(ix,0,iz-1,ip)] * lowerDiagonal.at(iz - 2);
	 }
				
	 gsl_linalg_solve_tridiag(diagonal, upperDiagonal, lowerDiagonal, knownTerm, solution);
			
	 for (iz = 1; iz < nz - 1; ++iz) {
	    value = solution.at(iz - 1);
				
	    N[index(ix,0,iz,ip)] = (value > 0) ? value : 0;
					
	 }
      }
    }
  }
}



void DRAGONCELLO::propagateMixedDerivatives(double dt, double dtbar, int counter) {
    
    for (int ip = 0; ip < np; ip++) {
        for (int ix = 1; ix < nx - 1; ix++) { //notice the boundary condition
            for (int iz = 1; iz < nz - 1; iz++) {
                
                long int ind = index(ix,0,iz,ip);
                
                double S_ij = 2.*D_rz.at(ind)/(4.*deltaxCentral[ix]*deltazCentral[iz]);
                    
                double value = N.at(ind);

                value += dt * S_ij * ( N.at(index(ix+1,0,iz+1,ip)) - N.at(index(ix-1,0,iz+1,ip))
                                     - N.at(index(ix+1,0,iz-1,ip)) + N.at(index(ix-1,0,iz-1,ip)) );
                if (counter <= timestepSourceIsActiveUntil) {
                    value += dtbar * sourceTerm.at(ind);
                }  

                N.at(ind) = (value > 0) ? value : 0.0;

            }
        }
    }
}



void DRAGONCELLO::propagateEloss(double dt, double dtbar) {
  
  vector<double> knownTerm_p(np, 0.0);
  vector<double> diagonal_p(np, 0.0);
  vector<double> upperDiagonal_p(np, 0.0);
  vector<double> lowerDiagonal_p(np, 0.0);
  vector<double> solution_p(np, 0.0);
  
  double value;
  
  // propagation in momentum
  
  for (int ix = 0; ix < nx-1; ix++) {
       for (int iz = 1; iz < nz-1; iz++) {
           for (int ip = 0; ip < np-1; ++ip) {
               
               long int ind = index(ix,0,iz,ip);
               diagonal_p[ip] = 1. + dt/2.*centralDiagonalEloss[ind]; 
               upperDiagonal_p[ip] = -dt/2.*upperDiagonalEloss[ind];
               lowerDiagonal_p[ip] =  0.;
               knownTerm_p[ip] = N[ind] * (2. - diagonal_p[ip]) -  N[ind+1]*upperDiagonal_p[ip] + dtbar*sourceTerm[ind];
            }
           
            solveTridiagonalSystem(lowerDiagonal_p, diagonal_p, upperDiagonal_p, knownTerm_p, solution_p, np-1);
           
            for (int ip = 0; ip < np; ++ip) {
               value = solution_p[ip];
               N[index(ix,0,iz,ip)] = (value > 0) ? value : 0.0;
               
            }
       }
   }
}



void DRAGONCELLO::propagateInRaniso(double dt, double dtbar, int counter) {

   int iz,ix,ip;
   long int ind;
   double value;

#ifdef PARALLEL   
#pragma omp parallel default(shared) private(iz,ix,ip,ind,value) num_threads(NUM_THREADS) 
#endif          
  {

    vector<double> knownTerm(nx - 1, 0.0);
    vector<double> diagonal(nx - 1, 0.0);
    vector<double> upperDiagonal(nx - 2, 0.0);
    vector<double> lowerDiagonal(nx - 2, 0.0);
    vector<double> solution(nx - 1, 0.0);

#ifdef PARALLEL
#pragma omp for schedule(dynamic)
#endif    	

   for (ip = 0; ip < np; ip++) {
      for (iz = 1; iz < nz - 1; iz++) { //notice the boundary condition
         for (ix = 0; ix < nx - 1; ix++) {
				
            ind = index(ix,0,iz,ip);
	
            diagonal.at(ix) = 1. + 0.5*dt * centralDiagonalR.at(ind);
	    if (ix != nx - 2) {
	       upperDiagonal.at(ix) = - 0.5*dt * upperDiagonalR.at(ind);
	    }
	    if (ix == 0) {
	       upperDiagonal.at(ix) += - 0.5*dt * lowerDiagonalR.at(ind);
	    }
	    if (ix != 0) {
	       lowerDiagonal.at(ix - 1) = - 0.5*dt * lowerDiagonalR.at(ind);
	    }
	    knownTerm.at(ix) = N.at(ind) * (2. - diagonal.at(ix));
	    if (counter <= timestepSourceIsActiveUntil) {
	       knownTerm.at(ix) += dtbar * sourceTerm.at(ind);
	    }    
	    if (ix != nx - 2)
	       knownTerm.at(ix) -= N[index(ix+1,0,iz,ip)] * upperDiagonal.at(ix);
	    if (ix != 0) 
	       knownTerm.at(ix) -= N[index(ix-1,0,iz,ip)] * lowerDiagonal.at(ix - 1);
	 }
				
	 gsl_linalg_solve_tridiag(diagonal, upperDiagonal, lowerDiagonal, knownTerm, solution);
	
	 for (ix = 0; ix < nx - 1; ++ix) {
	    value = solution.at(ix);
	    N[index(ix,0,iz,ip)] = (value > 0) ? value : 0.0;
	 }
      }
   }
  }
}



void DRAGONCELLO::printLocalFlux() {

  double norm = localFlux / (N[index(ixsun,0,izsun,ipLocal)] / ( c_light/4./M_PI * GeV  ) );
 
  cout << "local Flux = " << localFlux << " 1/cm^2/s/sr/GeV " << endl;
  cout << "local Flux from the code = " << N[index(ixsun,0,izsun,ipLocal)] / ( c_light/4./M_PI * GeV  ) << " 1/cm^2/s/sr/GeV " << endl;
  cout << "Norm [NOT APPLIED] = " << norm << endl;
}



void DRAGONCELLO::dumpSolution(int counter) {
  
  stringstream sstream, sstream2;
  string filename, filename2;
  ofstream outfile, outfile2;
    
  sstream << "output/localSpectrum_" << counter << ".dat" ;
  filename = sstream.str();
  /*outfile.open(filename.c_str());
  outfile << "p" << "\t" << "N" << "\t" << "N_analytical" << endl;

  for (int ip = 0; ip < np; ip++) {
	  outfile << p_vec[ip] / GeV << "\t" << N[index(ixsun,0,izsun,ip)] * cm3 * GeV  << "\t" << analyticalSolution[index(ixsun,0,izsun,ip)] * cm3 * GeV  << endl; 
	  //outfile << p_vec[ip] << "\t" << N[index(ixsun,0,izsun,ip)] << "\t" << analyticalSolution[index(ixsun,0,izsun,ip)] << endl; //NOUNITS

	  if (ip%4==0)
		  cout << ip << "\t" << p_vec[ip] / GeV << " GeV \t" << N[index(ixsun,0,izsun,ip)] * cm3 * GeV << endl;
	  //cout << p_vec[ip]  << " GeV \t" << N[index(ixsun,0,izsun,ip)]  << endl; //NOUNITS
  }
  outfile.close();*/
    
  sstream2 << "output/2DAnisoMap_" << counter << ".dat" ;
  filename2 = sstream2.str();
  cout << "Writing partial output on this file: " << filename2 << endl;
  outfile2.open(filename2.c_str());
  outfile2 << "#r" << "   " << "z" << "   " << "N " << endl;
  outfile2 << scientific << setprecision(2);
  for (int ix = 0; ix < nx; ++ix) 
	  for (int iz = 0; iz < nz; ++iz) {
		outfile2 << x_vec[ix] / kpc << "   ";
		outfile2 << z_vec[iz] / kpc << "   ";
		outfile2 << N[index(ix,0,iz,0)] * cm3 * GeV << "   " << endl;
	  }
  outfile2.close();

}



void DRAGONCELLO::dumpSpectra_2(int counter) {
  
  stringstream sstream3;
  string filename3;
  ofstream outfile3;
    
   
  sstream3 << "output/Spectra_r-p_z=0_" << counter << ".dat" ;
  filename3 = sstream3.str();
  cout << "Writing partial output on this file: " << filename3 << endl;
  outfile3.open(filename3.c_str());
  outfile3 << "#r" << "   " << "p" << "   " << "N" << endl;
  outfile3 << scientific << setprecision(2);
  for (int ix = 0; ix < nx; ++ix) 
	  for (int ip = 0; ip < np; ++ip) {
		outfile3 << x_vec[ix] / kpc << "   ";
		outfile3 << p_vec[ip] / GeV << "   ";
		outfile3 << N[index(ix,0,(nz-1)/2,ip)] * cm3 * GeV << "   " << endl;
	  }
  outfile3.close();

}



void DRAGONCELLO::compute_epsilon(int counter) {

    stringstream sstream;
    string filename;
    ofstream outfile;
    
    for (int ip = 0; ip < np; ip++) {
        sstream << "output/epsilon_dt="<< dtmax/Myr << "_Nz_=" << nz << "_p=" << p_vec[ip]/GeV << "_" << counter << ".dat" ;
        filename = sstream.str();
        cout << "Writing partial output on this file: " << filename << endl;
        outfile.open(filename.c_str());
        outfile << "p" << "\t" << "dz" << "\t" << "eps" << endl;
        double dz = (zmax - zmin) / (nz - 1);
        double epsilon;
        epsilon = 0.;
        for (int iz = 1; iz < nz-1; iz++) {
            double deltaN = abs(N[index(ixsun,0,iz,ip)] - analyticalSolution[index(ixsun,0,iz,ip)]);
            if (deltaN > epsilon) epsilon = deltaN;
        }
        outfile << p_vec[ip] / GeV << "\t" << dz/kpc << "\t" << epsilon << endl;
        outfile.close();
        sstream.str(std::string());
    }
}



void DRAGONCELLO::computeResidualDiffAniso(int counter) {
    
    vector<norms> residual;

    int iz, ix, ip;
    long int ind;
    long npoints;     

    double meanResidual;
    double maxResidual;
    double norm2Residual;

    double Q;

    double Drr;
    double Drz;
    double Dzz;
    double Drr_up;
    double Drr_do;
    double Dzz_up;
    double Dzz_do;

    double Ni;
    double Nrup;
    double Nrdown;
    double Nzup;
    double Nzdown;
    double Nrupzup;
    double Nrupzdown;
    double Nrdownzup;
    double Nrdownzdown;

    double r;
    double rup;
    double rdown;
    double z;
    double zup;
    double zdown;

    double Lr;
    double Lz;
    double Lmix;

    double Residual;
   

    for (ip = 0; ip < np; ++ip) {
        
       meanResidual = 0.;
       maxResidual = 0.;
       norm2Residual = 0.;
        
       npoints = 0;     

       for (ix = 1; ix < nx-1; ++ix)
          for (iz = 2; iz < nz-2; ++iz) {
            
             ind = index(ix,0,iz,ip);
            
             Q = sourceTerm[ind];
            
             if (fabs(Q) > 1e-60) {
                
                Drr    = D_rr[ind];
                Drz    = D_rz[ind];
                Dzz    = D_zz[ind];
                Drr_up = D_rr[index(ix+1,0,iz,ip)];
                Drr_do = D_rr[index(ix-1,0,iz,ip)];
                Dzz_up = D_zz[index(ix,0,iz+1,ip)];
                Dzz_do = D_zz[index(ix,0,iz-1,ip)];
                
                Ni          = N[ind];
                Nrup        = N[index(ix+1,0,iz,ip)];
                Nrdown      = N[index(ix-1,0,iz,ip)];
                Nzup        = N[index(ix,0,iz+1,ip)];
                Nzdown      = N[index(ix,0,iz-1,ip)];
                Nrupzup     = N[index(ix+1,0,iz+1,ip)];
                Nrupzdown   = N[index(ix+1,0,iz-1,ip)];
                Nrdownzup   = N[index(ix-1,0,iz+1,ip)];
                Nrdownzdown = N[index(ix-1,0,iz-1,ip)];

                r     = x_vec[ix];
                rup   = x_vec[ix+1];
                rdown = x_vec[ix-1];
                z     = z_vec[iz];
                zup   = z_vec[iz+1];
                zdown = z_vec[iz-1];
                
                Lr = Drr / ((rup - rdown) / 2.) * ((Nrup - Ni) / (rup - r) - (Ni - Nrdown) / (r - rdown))
                + Drr / r * ((Nrup - Nrdown) / (rup - rdown))
                + (Drr_up - Drr_do) / (rup - rdown) * (Nrup - Nrdown) / (rup - rdown);
                
                Lz = Dzz / ((zup - zdown) / 2.) * ((Nzup - Ni) / (zup - z) - (Ni - Nzdown) / (z - zdown))
                + (Dzz_up - Dzz_do) / (zup - zdown) * (Nzup - Nzdown) / (zup - zdown);

                Lmix = Drz*(Nrupzup - Nrupzdown - Nrdownzup + Nrdownzdown)/( (rup - rdown)*(zup - zdown) / 2. );
                
                Residual =  (Lz + Lr + Lmix + Q) / Q;
                
                meanResidual += fabs(Residual);
                norm2Residual += Residual*Residual;
                if (Residual > maxResidual) maxResidual = Residual;
                
                npoints++;
             }
         }   

        norms n = { maxResidual, meanResidual / (double)npoints, sqrt(norm2Residual) / (double)npoints };
        
        if (ip%10==0)

        cout  << p_vec[ip] / GeV << " GeV -> "  << n.max << "\t" << n.mean << "\t" << n.norm2 << endl;
        
        residual.push_back(n);
        cout << ip << "\t" << n.max << "\t" << residual[0].max << endl;
        // }
        
        stringstream sstream;
        string filename;
        ofstream outfile;
        
        sstream << "output/NormResidual_dt="<< dtmax/Myr << "_Nz=" << nz << "_p=" << p_vec[ip]/GeV << ".dat" ;
        filename = sstream.str();
        cout << "Writing residuals on this file: " << filename << endl;
        
        outfile.open(filename.c_str(),std::ios_base::app);
        if (counter == 0 ) outfile << "p" << "\t" << "Max" << "\t" << "Mean" << "\t" << "Norm2" << endl;
        
        outfile << p_vec[ip] / GeV << "\t" << residual[ip].max << "\t" << residual[ip].mean << "\t" << residual[ip].norm2 << endl;
        
        outfile.close();
        sstream.str(std::string());
    }

}



void DRAGONCELLO::computeResidualLosses(int counter) {
    
    //vector<norms> residual(nx*nz);
    
    stringstream sstream;
    string filename;
    ofstream outfile;
    
    sstream << "output/normResidual_" << counter << ".dat" ;
    filename = sstream.str();
    cout << "Writing residuals on this file: " << filename << endl;

    outfile.open(filename.c_str());
    outfile << "ip" << "\t" << "p" << "\t" << "Max" << "\t" << "Mean" << "\t" << "Norm2" << endl;
    
    outfile << scientific << setprecision(5);
    
    for (int ix = 1; ix < nx-1; ++ix)
        for (int iz = 1; iz < nz-1; ++iz) {
            
            if (!((ix==nx/2) && (iz==nz/2)))
                continue;
        
            double meanResidual  = 0.;
            //double maxResidual   = 0.;
            double norm2Residual = 0.;

            long npoints = 0;
    
            for (int ip = 1; ip < np-1; ++ip) {
        
                long int ind = index(ix,0,iz,ip);
                
                double Q = sourceTerm[ind];
                
                if (fabs(Q) > 1e-60) {
                    
                    double Ni     = N[ind];
                    double Npup   = N[index(ix,0,iz,ip+1)];
                    //double Npdown = N[index(ix,0,iz,ip-1)];
                    
                    double Pdot   = energyLossTerm[index(ix,0,iz,ip)];
                    double Pdotup   = energyLossTerm[index(ix,0,iz,ip+1)];
                    //double Pdotdown   = energyLossTerm[index(ix,0,iz,ip-1)];
                    
                    double p      = p_vec[ip];
                    double pup    = p_vec[ip+1];
                    //double pdown  = p_vec[ip-1];
       
                    double Lp = ( ( Pdotup*Npup - Pdot*Ni )/( pup - p ) );
                    
                    double Residual =  (Lp + Q) / Q;
                    
                    meanResidual = fabs(Residual);
                    norm2Residual = Residual*Residual;
                    
                    npoints = 1;
                    
                    norms n = { Residual, meanResidual / (double)npoints, sqrt(norm2Residual) / (double)npoints };
        
                    if ((ix==nx/2) && (iz==nz/2)) {
                    	if (ip%8==0)
                          cout << ip << "\t" << p_vec[ip]/GeV  << "\t" << n.max << "\t" << n.mean << "\t" << n.norm2 << endl;
                        outfile << ip << "\t" << p_vec[ip]/GeV  << "\t" << n.max << "\t" << n.mean << "\t" << n.norm2 << endl;                
                    }
                    
                    //npoints++;
                }
            }           
        
    }
    
    outfile.close();
}

  

void DRAGONCELLO::computeSolutionEnergyIntegral(int ix_, int iy_, int iz_){
	
	solutionEnergyIntegral = 0.;
	for (int ip = 0; ip < np; ++ip) 
		solutionEnergyIntegral += (N[index(ix_,iy_,iz_,ip)]*p_vec[ip]*p_vec[ip]*pFactor);
				  
}



void DRAGONCELLO::computeSourceEnergyIntegral(int ix_, int iy_, int iz_){
	
	sourceEnergyIntegral = 0.;
	for (int ip = 0; ip < np; ++ip) 
		sourceEnergyIntegral += (sourceTerm[index(ix_,iy_,iz_,ip)]*p_vec[ip]*p_vec[ip]*pFactor);
				  
}







void DRAGONCELLO::Run2D() {
  
  cout << endl << "Welcome to the 2D solver for Particle: A = " << A << "\t" << "Z = " << Z << endl << endl;
  
  double dt = dtmax;
  double currentTime = 0.;
  
  int counter = -1;    
  int totalCounter = -1;
    
  while ((dt > dtmin) && (counter < interruptAfter)) {
    
    counter++;
    
    double dtbar = dt / p;

    cout << endl << "   *** dt = " << dt / Myr << " Myr & Nrept = " << Nrept << endl;
    cout << "   *** Current time = " << currentTime/Myr << " Myr" << endl;
    
    dumpSolution(counter);
    dumpSpectra_2(counter);

    time_t time_s, time_e;
    time(&time_s); // fix initial time


    cout << "-> Propagation in (Z,R): begins..." << endl;
    
    for (int Niter = 0; Niter < Nrept; ++Niter) {
    	
      totalCounter++;
      currentTime += dt;
      
      propagateInZaniso(dt, dtbar, totalCounter);
      propagateMixedDerivatives(dt, dtbar, totalCounter);
      propagateInRaniso(dt, dtbar, totalCounter);

      //propagateEloss(dt, dtbar); //IC-losses for leptons only (un-comment if needed)

    }
    cout << "   ...END" << endl;
      
    time(&time_e); // fix final time

    computeResidualDiffAniso(counter);      

    cout << Nrept << "iterations done in " << (double)(time_e-time_s) << " s." << endl;

    printLocalFlux();
      
    dt *= dtfactor;
    
  } 
  
  
  return ;
}
