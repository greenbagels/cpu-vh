/*
 * LatticeParameters.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef LATTICEPARAMETERS_H_
#define LATTICEPARAMETERS_H_

#include <string>
#include <libconfig.h++>
#include "parameters.h"

#define N_GHOST_CELLS_M 2
#define N_GHOST_CELLS_P 2
#define N_GHOST_CELLS 4

namespace rhic
{
	struct lattice_parameters : parameters
	{
		lattice_parameters(libconfig::Config &cfg, std::string config_dir);
		void load_params(libconfig::Config &cfg, const char *config_dir) override;
	
		int numLatticePointsX;
		int numLatticePointsY;
		int numLatticePointsRapidity;
		int numComputationalLatticePointsX;
		int numComputationalLatticePointsY;
		int numComputationalLatticePointsRapidity;
		int numProperTimePoints;

		double latticeSpacingX;
		double latticeSpacingY;
		double latticeSpacingRapidity;
		double latticeSpacingProperTime;
	};
}
#endif /* LATTICEPARAMETERS_H_ */
