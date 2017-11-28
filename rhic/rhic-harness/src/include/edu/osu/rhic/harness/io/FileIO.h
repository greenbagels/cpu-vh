/*
 * FileIO.h
 *
 *  Created on: Oct 24, 2015
 *      Author: bazow
 */

#ifndef FILEIO_H_
#define FILEIO_H_

#include <fstream>
#include <sstream>
#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

namespace rhic
{
	class file_io
	{
		public:
			void output(const PRECISION * const var, double t, const char *pathToOutDir, const char *name, void * latticeParams);
		private:
			std::stringstream fname;
			std::ofstream file;
	}
}

#endif /* FILEIO_H_ */
