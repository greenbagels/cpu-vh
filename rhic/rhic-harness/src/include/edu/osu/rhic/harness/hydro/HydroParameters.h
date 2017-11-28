/*
 * HydroParameters.h
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#ifndef HYDROPARAMETERS_H_
#define HYDROPARAMETERS_H_

#include <libconfig.h++>
#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"
#include "parameters.h"

namespace rhic
{
	struct hydro_parameters : parameters
	{
		hydro_parameters(libconfig::Config &cfg, std::string &config_dir);
		void load_params(libconfig::Config &cfg, std::string &config_dir) override;  

		double initialProperTimePoint;
		double shearViscosityToEntropyDensity;
		double freezeoutTemperatureGeV;
		int initializePimunuNavierStokes;
	};
}

#endif /* HYDROPARAMETERS_H_ */
