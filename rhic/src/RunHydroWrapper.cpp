//this is an example of how the HYDRO class can be instantiated in a larger program
//by default this file is compiled, which can then be run in the same way that the standalone executable
//given by Run.cpp was run...

//by default the Makefile produces the executable cpu-vh

#include "HydroWrapper.h"

int main(int argc, char **argv)
{
  //Declare an instance of HYDRO class
  HYDRO vh;

  //pass the initial energy density vector
  //vh.initialize_ed_from_vector(init_e);

  //run the hydro until freezeout
  vh.run_hydro(argc, argv);

  //save the fo surface info to a struct which can be passed to CF/Sampler module

}
