#ifndef QSTATEPOTTS_H
#define QSTATEPOTTS_H

#include <random>
#include <vector>
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_traits.hpp"

using namespace boost;
using namespace std;

class QStatePotts
{
public:
    QStatePotts(const string& net_file,const uint8_t& num_states,const int& rank);
    mt19937 rng;

    uint8_t num_states;

    double J;
    double T;

    void mcstep();
    void mcs();

    double interaction(const uint8_t& s1,const uint8_t& s2)
    {
        return J * (s1 == s2 ? 1: 0);
    }

    Graph scale_free_net;
    vector<vertex_t> vertices;

    typedef adjacency_list<vecS, vecS, undirectedS,uint8_t> Graph;
    typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef boost::graph_traits<Graph>::edge_descriptor edge_t;
};

#endif // QSTATEPOTTS_H
