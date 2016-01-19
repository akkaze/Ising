#ifndef SCHED_HPP
#define SCHED_HPP
#include <boost/mpi/communicator.hpp>

#include "sim_datastruct.hpp"


namespace mpi = boost::mpi;
void sched_master( const sim_parameters& par, const mpi::communicator& mpicomm );
void sched_slave( const sim_parameters& par, const mpi::communicator& mpicomm );
#endif // SCHED_HPP

