/*
 * Properties.h
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#ifndef PROPERTIES_H_
#define PROPERTIES_H_

#include <libconfig.h++>

namespace rhic
{
	template<typename T>
	void get_prop(libconfig::Config &cfg, const char *prop_name, T &prop_val, T default_val);
}

#endif /* PROPERTIES_H_ */
