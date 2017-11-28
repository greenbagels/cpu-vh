/*
 ============================================================================
 Name        : Run.c
 Author      : Dennis Bazow
 Version     :
 Copyright   :
 Description : Run viscous hydrodynamic simulation of a relativistic heavy ion collision
 ============================================================================
 */

#include <stdlib.h>
#include <stdio.h> // for printf
#include <sys/time.h> // for timing
//#include <unistd.h>		// for current working directory
#include <libconfig.h++> // not a typo, switch to the C++ API

//#include "gtest/gtest.h" // for unit testing

#include "edu/osu/rhic/harness/cli/CommandLineArguments.h"
#include "edu/osu/rhic/harness/lattice/LatticeParameters.h"
#include "edu/osu/rhic/harness/ic/InitialConditionParameters.h"
#include "edu/osu/rhic/harness/hydro/HydroParameters.h"
#include "edu/osu/rhic/harness/hydro/HydroPlugin.h"

/*
int runTest(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
*/

void runHydro(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory, const char *outputDir) {
	run(latticeParams, initCondParams, hydroParams, rootDirectory, outputDir);
}

int main(int argc, char **argv) {

	class cli_arguments cli;
	struct LatticeParameters latticeParams;
	struct InitialConditionParameters initCondParams;
	struct HydroParameters hydroParams;

	cli::load_cli_args(argc, argv);
	cli::print_arg_vals();

	//=========================================
	// Set parameters from configuration files
	//=========================================
	libconfig::Config lattice_config, ic_config, hydro_config;

	// TODO: Switch to C++ API
	// Set lattice parameters from configuration file
	load_lattice_params(lattice_config, cli::cli_args.config_dir, lattice_params);
	// Set initial condition parameters from configuration file
	load_ic_params(ic_config, cli::cli_args.config_dir, ic_params);
	// Set hydrodynamic parameters from configuration file
	load_hydro_params(hydro_config, cli::cli_args.config_dir, hydro_params);

	//=========================================
	// Run tests
	//=========================================
	/*
	if (cli.runTest) {
		int status = runTest(argc, argv);		
		printf("Done tests.\n");
	}
	*/
	//=========================================
	// Run hydro
	//=========================================
	if (cli.runHydro) {
		runHydro(&latticeParams, &initCondParams, &hydroParams, rootDirectory, cli.outputDirectory);
		printf("Done hydro.\n");
	}

	// TODO: Probably should free host memory here since the freezeout plugin will need
	// to access the energy density, pressure, and fluid velocity.

	return 0;
}
