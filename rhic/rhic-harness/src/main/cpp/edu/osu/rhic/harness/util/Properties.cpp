/*
 * Properties.c
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#include "edu/osu/rhic/harness/util/Properties.h"

namespace rhic
{
	template<typename T>
	void get_prop(libconfig::Config &cfg, const char *prop_name, T &prop_val, T default_val) {
		  if(cfg::lookupValue(std::string(prop_name), prop_val))
		    return;
		  else
		    prop_val = default_val;
	}
}

