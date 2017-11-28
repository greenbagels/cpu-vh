/*
 * HydroParameters.c
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#include "edu/osu/rhic/harness/hydro/HydroParameters.h"
#include "edu/osu/rhic/harness/util/Properties.h"

namespace rhic
{
	void hydro_parameters::hydro_parameters(libconfig::Config &cfg, std::string &config_dir)
	{
		load_params(cfg, config_dir);
	}

	void hydro_parameters::load_params(libconfig::Config &cfg, std::string &config_dir) {
		// Read the file
		std::string fname = config_dir + std::string("/hydro.properties");
		
		try
		{
			cfg.readFile(fname.c_str());
		}
		catch (const std::exception &e)
		{
			std::cerr << "No configuration file " << fname << " found for hydrodynamic parameters - "
			  << config_error_text(cfg) << "\nUsing default hydrodynamic configuration parameters.\n";
		}
	
		get_prop(cfg, "initialProperTimePoint", &initialProperTimePoint, 0.1);
		get_prop(cfg, "shearViscosityToEntropyDensity", &shearViscosityToEntropyDensity, 0.0795775);
		get_prop(cfg, "freezeoutTemperatureGeV", &freezeoutTemperatureGeV, 0.155);
	
		get_prop(cfg, "initializePimunuNavierStokes", &initializePimunuNavierStokes, 1);
	}
}
