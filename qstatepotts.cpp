#include "qstatepotts.h"
#include <chrono>

QStatePotts::QStatePotts(const string &net_file,const uint8_t& num_states,const int& rank)
{
    unsigned int rngseed;
    rngseed = chrono::system_clock::now().time_since_epoch().count();
    rngseed += rank;
    this->rng = mt19937(rngseed);

    this->num_states = num_states;
}




