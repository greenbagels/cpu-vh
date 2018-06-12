//This is a c++ wrapper for cpu-vh
//Written by Derek Everett 2018

//Test if this wrapper can also be used with cpu-vah and gpu-vh

#ifndef SRC_HYDROWRAPPER_
#define SRC_HYDROWRAPPER_

#include <stdlib.h>
#include <stdio.h> // for printf
#include <sys/time.h> // for timing
#include <unistd.h>		// for current working directory
#include <libconfig.h>
#include <vector>

#include "../include/CommandLineArguments.h"
#include "../include/LatticeParameters.h"
#include "../include/InitialConditionParameters.h"
#include "../include/HydroParameters.h"
#include "../include/HydroPlugin.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

const char *version = "";
const char *address = "";

class HYDRO {
 private:

 public:
    HYDRO();
    ~HYDRO();

    std::vector<double> initial_energy_density;
    std::vector<double> initial_pressure;
    std::vector<double> initial_ut;
    std::vector<double> initial_ux;
    std::vector<double> initial_uy;
    std::vector<double> initial_un;
    std::vector<double> initial_pitt;
    std::vector<double> initial_pitx;
    std::vector<double> initial_pity;
    std::vector<double> initial_pitn;
    std::vector<double> initial_pixx;
    std::vector<double> initial_pixy;
    std::vector<double> initial_pixn;
    std::vector<double> initial_piyy;
    std::vector<double> initial_piyn;
    std::vector<double> initial_pinn;
    std::vector<double> initial_Pi;

    //run the hydro until freezeout
    int run_hydro();


    //run the hydro for one time step
    //int run_hydro_time_step();

    //save the freezeout surface info
    //void save_fo_history();

    //support to initilialize the energy density from a vector - useful for JETSCAPE
    //note units of argument should be GeV / fm^3
    //then we convert to fm^(-4)
    void initialize_ed_from_vector(std::vector<double>);


    //support to initilialize all components of T^\mu\nu from vectors - useful for JETSCAPE
    void initialize_from_vectors(std::vector<double>&, //e
                            std::vector<double>&, //p
                            std::vector<double>&, //ut
                            std::vector<double>&, //ux
                            std::vector<double>&, //uy
                            std::vector<double>&, //un
                            std::vector<double>&, //pitt
                            std::vector<double>&, //pitx
                            std::vector<double>&, //pity
                            std::vector<double>&, //pitn
                            std::vector<double>&, //pixx
                            std::vector<double>&, //pixy
                            std::vector<double>&, //pixn
                            std::vector<double>&, //piyy
                            std::vector<double>&, //piyn
                            std::vector<double>&, //pinn
                            std::vector<double>&); //Pi
};

HYDRO::HYDRO() {

}

HYDRO::~HYDRO() {
}

//use this function to initialize energy density
void HYDRO::initialize_ed_from_vector(std::vector<double> energy_density_in) {
  initial_energy_density = energy_density_in;
}

//use this function to return final hydro variables as vectors
void HYDRO::initialize_from_vectors(std::vector<double> &energy_density_in,
                                        std::vector<double> &pressure_in,
                                        std::vector<double> &ut_in,
                                        std::vector<double> &ux_in,
                                        std::vector<double> &uy_in,
                                        std::vector<double> &un_in,
                                        std::vector<double> &pitt_in,
                                        std::vector<double> &pitx_in,
                                        std::vector<double> &pity_in,
                                        std::vector<double> &pitn_in,
                                        std::vector<double> &pixx_in,
                                        std::vector<double> &pixy_in,
                                        std::vector<double> &pixn_in,
                                        std::vector<double> &piyy_in,
                                        std::vector<double> &piyn_in,
                                        std::vector<double> &pinn_in,
                                        std::vector<double> &Pi_in) {
  initial_energy_density = energy_density_in;
  initial_pressure = pressure_in;
  initial_ut = ut_in;
  initial_ux = ux_in;
  initial_uy = uy_in;
  initial_un = un_in;
  initial_pitt = pitt_in;
  initial_pitx = pitx_in;
  initial_pity = pity_in;
  initial_pitn = pitn_in;
  initial_pixx = pixx_in;
  initial_pixy = pixy_in;
  initial_pixn = pixn_in;
  initial_piyy = piyy_in;
  initial_piyn = piyn_in;
  initial_pinn = pinn_in;
  initial_Pi = Pi_in;
}


//where the magic happens
//taken from main() in Run.cpp of cpu-vh
int HYDRO::run_hydro(int argc, char **argv)
{

  struct CommandLineArguments cli;
	struct LatticeParameters latticeParams;
	struct InitialConditionParameters initCondParams;
	struct HydroParameters hydroParams;

	loadCommandLineArguments(argc, argv, &cli, version, address);

	char *rootDirectory = NULL;
	size_t size;
	rootDirectory = getcwd(rootDirectory,size);

	// Print argument values
	printf("configDirectory = %s\n", cli.configDirectory);
	printf("outputDirectory = %s\n", cli.outputDirectory);
	if (cli.runHydro)
		printf("runHydro = True\n");
	else
		printf("runHydro = False\n");
	if (cli.runTest)
		printf("runTest = True\n");
	else
		printf("runTest = False\n");

	//=========================================
	// Set parameters from configuration files
	//=========================================
	config_t latticeConfig, initCondConfig, hydroConfig;

	// Set lattice parameters from configuration file
	config_init(&latticeConfig);
	loadLatticeParameters(&latticeConfig, cli.configDirectory, &latticeParams);
	config_destroy(&latticeConfig);
	// Set initial condition parameters from configuration file
	config_init(&initCondConfig);
	loadInitialConditionParameters(&initCondConfig, cli.configDirectory, &initCondParams);
	config_destroy(&initCondConfig);
	// Set hydrodynamic parameters from configuration file
	config_init(&hydroConfig);
	loadHydroParameters(&hydroConfig, cli.configDirectory, &hydroParams);
	config_destroy (&hydroConfig);

	//=========================================
	// Run hydro
	//=========================================
	if (cli.runHydro) {
    run(&latticeParams, &initCondParams, &hydroParams, rootDirectory, cli.outputDir);
		printf("Done hydro.\n");
	}

	// TODO: Probably should free host memory here since the freezeout plugin will need
	// to access the energy density, pressure, and fluid velocity.

	return 0;

}

#endif  // SRC_HYDROWRAPPER_
