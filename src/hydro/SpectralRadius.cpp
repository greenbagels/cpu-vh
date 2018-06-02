/*
 * SpectralRadius.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#include <math.h>

#include "../hydro/SpectralRadius.h"
#include "../hydro/EnergyMomentumTensor.h"
#include "../hydro/DynamicalVariables.h"

PRECISION spectralRadiusX(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un) {
	return fabs(ux/ut);
}

PRECISION spectralRadiusY(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un) {
	return fabs(uy/ut);
}

PRECISION spectralRadiusZ(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un) {
	return fabs(un/ut);
}
