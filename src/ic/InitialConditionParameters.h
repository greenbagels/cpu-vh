/*
 * InitialConditionParameters.h
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#ifndef INITIALCONDITIONPARAMETERS_H_
#define INITIALCONDITIONPARAMETERS_H_

#include <libconfig.h>

struct InitialConditionParameters
{
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

void loadInitialConditionParameters(config_t *cfg, const char* configDirectory, void * params);

#endif /* INITIALCONDITIONPARAMETERS_H_ */
