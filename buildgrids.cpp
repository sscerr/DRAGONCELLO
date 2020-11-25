#include "dragoncello.h"

void DRAGONCELLO::buildXGrid() {
	cout << "x in " << xmin / kpc << " ... " << xmax / kpc << " kpc" << endl;

	deltax = (xmax - xmin)/(nx-1);
	for (int ix=0; ix<nx; ix++) {
		x_vec.push_back(xmin + ix * deltax);
	}

	ixsun = lower_bound(x_vec.begin(), x_vec.end(), x_Sun) - x_vec.begin();

	for (int ix = 0; ix < nx; ++ix) {
		if (ix < nx-1)
			deltaxUp.push_back(x_vec[ix+1] - x_vec[ix]);
		else
			deltaxUp.push_back(x_vec[ix] - x_vec[ix-1]);

		if (ix > 0)
			deltaxDown.push_back(x_vec[ix] - x_vec[ix-1]);
		else
			deltaxDown.push_back(x_vec[1] - x_vec[0]);

		deltaxCentral.push_back(0.5 * (deltaxDown.back() + deltaxUp.back()));
	}
}



void DRAGONCELLO::buildYGrid() {
        cout << "y in " << ymin / kpc << " ... " << ymax / kpc << " kpc" << endl;

        if (ny > 1) {
                cout << "   !!!   W A R N I N G   !!!  " << endl;
                cout << "  - - - - - - - - - - - - - - " << endl;
                cout << "   3D version NOT AVAILABLE   " << endl;
                exit(0);
                //deltay = (ymax - ymin) / (ny - 1);
                //for (int iy=0; iy<ny; iy++)
                //        y_vec.push_back(ymin + iy*deltay);
        }
        else {
                deltay = 0;
                y_vec.push_back(0);
        }

        iysun = lower_bound(y_vec.begin(), y_vec.end(), y_Sun) - y_vec.begin();

        for (int iy=0; iy<ny; ++iy) {
                if (iy < ny-1)
                        deltayUp.push_back(y_vec[iy+1] - y_vec[iy]);
                else
                        deltayUp.push_back(y_vec[iy] - y_vec[iy-1]);

                if (iy > 0)
                        deltayDown.push_back(y_vec[iy] - y_vec[iy-1]);
                else
                        deltayDown.push_back(y_vec[1] - y_vec[0]);

                deltayCentral.push_back(0.5 * (deltayDown.back() + deltayUp.back()));
        }
}



void DRAGONCELLO::buildZGrid() {
	cout << "z in " << zmin / kpc << " ... " << zmax / kpc << " kpc" << endl;

	deltaz = (zmax - zmin) / (nz - 1);
	for (int iz=0; iz<nz; ++iz) {
		z_vec.push_back(zmin + iz * deltaz);
		//cout << z_vec[iz]/kpc << endl;
	}

	izsun = lower_bound(z_vec.begin(), z_vec.end(), z_Sun) - z_vec.begin();

	for (int iz = 0; iz < nz; ++iz) {
		if (iz < nz-1)
			deltazUp.push_back(z_vec[iz+1] - z_vec[iz]);
		else
			deltazUp.push_back(z_vec[iz] - z_vec[iz-1]); //CHECK!!!

		if (iz > 0)
			deltazDown.push_back(z_vec[iz] - z_vec[iz-1]);
		else
			deltazDown.push_back(z_vec[1] - z_vec[0]); //CHECK!!!

		deltazCentral.push_back(0.5 * (deltazDown.back() + deltazUp.back()));
	}
}

void DRAGONCELLO::buildPGrid() {
	cout << "p in " << pmin / GeV << " ... " << pmax / GeV << " GeV" << endl;

	deltaLogp = exp(log(pmax / pmin) / (np - 1));
	for (int ip = 0; ip < np; ip++) {
		p_vec.push_back(exp(log(pmin) + ip * log(deltaLogp)));
		cout << p_vec[ip]/GeV << endl;
	}

	ipLocal = lower_bound(p_vec.begin(), p_vec.end(), localFluxMomentum) - p_vec.begin();
}

void DRAGONCELLO::initEmptyGrids() {
	int nelements = nx * ny * nz * np;

	N.resize(nelements,0.);
	N_previous.resize(nelements,0);
	sourceTerm.resize(nelements,0);
        b_r.resize(nx*ny*nz,0);
        b_phi.resize(nx*ny*nz,0);
        b_x.resize(nx*ny*nz,0);
        b_y.resize(nx*ny*nz,0);
        b_z.resize(nx*ny*nz,0);
	Dpara.resize(nelements,0);
	Dperp.resize(nelements,0);
	D_rr.resize(nelements,0);
	D_rz.resize(nelements,0);
	D_xx.resize(nelements,0);
	D_xz.resize(nelements,0);
	D_zz.resize(nelements,0);
	u_r.resize(nelements,0);
	u_x.resize(nelements,0);
	u_z.resize(nelements,0);
	diffusionCoefficient.resize(nelements,0);
	diffusionCoefficient_xx.resize(nelements,0);
	diffusionCoefficient_zz.resize(nelements,0);
	reaccelerationCoefficient.resize(nelements,0);
	energyLossTerm.resize(nelements,0);
	analyticalSolution.resize(nelements,0.);

}

