/*
 * InitialConditionParameters.h
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#ifndef INITIALCONDITIONPARAMETERS_H_
#define INITIALCONDITIONPARAMETERS_H_

#include <libconfig.h++>
#include "parameters.h"

namespace rhic
{
	struct ic_parameters
	{
		ic_parameters(libconfig::Config &cfg, std::string config_dir) override;
		void load_params(libconfig::Config &cfg, std::string config_dir);
		int initialConditionType;
		int numberOfNucleonsPerNuclei;

		double initialEnergyDensity;
		double scatteringCrossSectionNN;
		double impactParameter;
		double fractionOfBinaryCollisions;
	
		// longitudinal energy density profile parameters
		double rapidityVariance; // \sigma^{2}_{\eta}
		double rapidityMean; // flat region around \ets_s = 0
	};
}

#endif /* INITIALCONDITIONPARAMETERS_H_ */
