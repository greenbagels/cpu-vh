/*
 * InitialConditionParameters.c
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#include "edu/osu/rhic/harness/ic/InitialConditionParameters.h"
#include "edu/osu/rhic/harness/util/Properties.h"

namespace rhic
{
	void ic_parameters::ic_parameters(libconfig::Config &cfg, std::string &config_dir)
	{
		load_params(cfg, config_dir);
	}
	
	void ic_parameters::load_params(libconfig::Config &cfg, std::string &config_dir) {
	// Read the file
		std::string fname = config_dir + std::string("/ic.properties");
		
		try
		{
			cfg.readFile(fname.c_str());
		}
		catch
		{
			std::cerr << "No configuration file " << fname << " found for initial condition parameters - "
			  << config_error_text(cfg) << "\nUsing default initial condition configuration parameters.\n";
		}

		get_prop(cfg, "initialConditionType", &initialConditionType, 2);	
		get_prop(cfg, "numberOfNucleonsPerNuclei", &numberOfNucleonsPerNuclei, 208);

		get_prop(cfg, "initialEnergyDensity", &initialEnergyDensity, 1.0);
		get_prop(cfg, "scatteringCrossSectionNN", &scatteringCrossSectionNN, 62);
		get_prop(cfg, "impactParameter", &impactParameter, 7);
		get_prop(cfg, "fractionOfBinaryCollisions", &fractionOfBinaryCollisions, 0.5);
		get_prop(cfg, "rapidityVariance", &rapidityVariance, 0.5);
		get_prop(cfg, "rapidityMean", &rapidityMean, 0.5);
	}
}
