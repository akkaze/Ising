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


void QStatePotts::mcstep()
{
    size_t num_vertices = vertices.size();

    uniform_int_distribution<> dis(0,num_vertices - 1);
    size_t vertex_idx = dis(rng);
    vertex_t vertex = vertices[vertex_idx];

    uint8_t prev_spin  = scale_free_net[vertex];
    uint8_t cur_spin = 0;

    dis = uniform_int_distribution<>(0,num_states - 1);
    do {
        cur_spin = dis(rng);
    }while(cur_spin != prev_spin);

    double prev_H = 0;
    double cur_H = 0;
}

