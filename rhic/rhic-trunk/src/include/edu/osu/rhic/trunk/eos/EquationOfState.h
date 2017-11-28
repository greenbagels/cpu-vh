/*
 * EquationOfState.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef EQUATIONOFSTATE_H_
#define EQUATIONOFSTATE_H_

#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

#define CONFORMAL_EOS

// ideal gas of massless quarks and gluons
//#define EOS_FACTOR 15.6269 // Nc=3, Nf=3
#define EOS_FACTOR 13.8997 // Nc=3, Nf=2.5

namespace rhic
{
	class eos
	{
		public:
			template<typename T>
			T equilibriumPressure(T e);

			template<typename T>
			T speedOfSoundSquared(T e);

			template<typename T>
			T effectiveTemperature(T e);

			template<typename T>
			T equilibriumEnergyDensity(T T);
	};
}

#endif /* EQUATIONOFSTATE_H_ */
