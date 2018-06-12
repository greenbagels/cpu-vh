#include "HydroWrapper.cpp"

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
