/*
 * FileIO.c
 *
 *  Created on: Oct 24, 2015
 *      Author: bazow
 */
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <cstdio>

#include "edu/osu/rhic/harness/io/FileIO.h"
#include "edu/osu/rhic/harness/lattice/LatticeParameters.h"
#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

namespace rhic
{
	void file_io::output(const precision_t var, double t, std::string &pathToOutDir, std::string &name, lattice_parameters &lattice_params) {
		// Clear fname's state
		fname.str(std::string());
		// New filename
		fname << pathToOutDir << "/" << name << std::fixed
					<< setprecision(3) << t << ".dat";

		file.open(fname); 

		int nx = lattice_params.numLatticePointsX;
		int ny = lattice_params.numLatticePointsY;
		int nz = lattice_params.numLatticePointsRapidity;
		double dx = lattice_params.latticeSpacingX;
		double dy = lattice_params.latticeSpacingY;
		double dz = lattice_params.latticeSpacingRapidity;

		double x,y,z;

		int i,j,k;
		int s;
	/*
		for(i = 2; i < nx+2; ++i) {
			x = (i-2 - (nx-1)/2.)*dx;
			for(j = 2; j < ny+2; ++j) {
				y = (j-2 - (ny-1)/2.)*dy;
				offset = (nz+4) * (j + (ny+4) * i);
				for(k = 2; k < nz+2; ++k) {
					s = k + offset;
					z = (k-2 - (nz-1)/2.)*dz;
					fprintf(fp, "%.3f\t%.3f\t%.3f\t%.8f\n",x,y,z,var[s]);
				}
			}
		}
	*/
		for(k = 2; k < nz+2; ++k) {
			z = (k-2 - (nz-1)/2.)*dz;
			for(j = 2; j < ny+2; ++j) {
				y = (j-2 - (ny-1)/2.)*dy;
				for(i = 2; i < nx+2; ++i) {
					x = (i-2 - (nx-1)/2.)*dx;
					s = col_maj_linear_index(i, j, k, nx+4, ny+4);
					file << std::fixed << setprecision(3) << x << "\t" << y << "\t" << z
							 << "\t" << var[s] << "\n";
				}
			}
		}

		file.close();
	}
}
