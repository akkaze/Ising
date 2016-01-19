#ifndef SIM_DATASTRUCTURE_HPP
#define SIM_DATASTRUCTURE_HPP
// ----- SIMULATION DATA STRUCTURES -----

struct sim_parameters {
  // ----- INPUT STRUCT FOR SIMRUN() -----

  unsigned int N;
  bool periodic;

  char init;
  unsigned int drysweeps;

  unsigned int bins;
  unsigned int binwidth;
  unsigned int sweeps;
  unsigned int use_fsize_correction;

  double J;
  double B;
  double T;
};


struct sim_results {
  // ----- OUTPUT STRUCT FOR SIMRUN() -----

  bool success;

  double B;
};


struct bin_results {
  // ----- RESULTS OF A SINGLE BIN -----
    double q2,q4;
};

#endif
