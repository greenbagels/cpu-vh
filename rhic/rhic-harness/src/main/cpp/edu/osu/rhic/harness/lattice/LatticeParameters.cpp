/*
 * LatticeParameters.c
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#include "edu/osu/rhic/harness/lattice/LatticeParameters.h"
#include "edu/osu/rhic/harness/util/Properties.h"

namespace rhic
{
	void lattice_parameters::lattice_parameters(libconfig::Config &cfg, std::string config_dir)
	{
		load_params(cfg, config_dir);
	}

	void lattice_parameters::load_params(libconfig::Config &cfg, std::string config_dir)
	{
		// Read the file
		std::string fname = config_dir + std::string("/lattice.properties");

		try
		{
			cfg.readFile(fname.c_str());
		}
		catch (const std::exception &e)
		{
			std::cerr << "No configuration file " << fname << " found for lattice parameters - "
			  << config_error_text(cfg) << "\nUsing default lattice configuration parameters.\n";
		}

		get_prop(cfg, "numLatticePointsX", &numLatticePointsX, 128);
		get_prop(cfg, "numLatticePointsY", &numLatticePointsY, 128);
		get_prop(cfg, "numLatticePointsRapidity", &numLatticePointsRapidity, 64);
		get_prop(cfg, "numProperTimePoints", &numProperTimePoints, 10);
	
		get_prop(cfg, "latticeSpacingX", &latticeSpacingX, 0.08);
		get_prop(cfg, "latticeSpacingY", &latticeSpacingY, 0.08);
		get_prop(cfg, "latticeSpacingRapidity", &latticeSpacingRapidity, 0.3);
		get_prop(cfg, "latticeSpacingProperTime", &latticeSpacingProperTime, 0.01);

		numComputationalLatticePointsX = numLatticePointsX+N_GHOST_CELLS;
		numComputationalLatticePointsY = numLatticePointsY+N_GHOST_CELLS;
		numComputationalLatticePointsRapidity = numLatticePointsRapidity+N_GHOST_CELLS;
	}
}

