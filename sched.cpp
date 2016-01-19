#include "sched.hpp"
#include "simrun.hpp"
#include "msgtag.hpp"

#include <boost/mpi/collectives.hpp>
void sched_slave( const sim_parameters& par, const mpi::communicator& mpicomm )
{
    schedmsg_t schedmsg;
    do {
        mpi::broadcast( mpicomm, schedmsg, 0 );

        if ( schedmsg == SCHEDMSG_START_MCC ) {
            // check if master wants us to initialize the system in a specific pconf
            simrun_slave( par, mpicomm );
        }
    } while ( schedmsg != SCHEDMSG_EXIT );

    return;
}


void sched_master( const sim_parameters& par, const mpi::communicator& mpicomm )
{
  // everything done, tell everyone to quit!
  schedmsg_t schedmsg;
  schedmsg = SCHEDMSG_START_MCC;
  mpi::broadcast( mpicomm, schedmsg, 0 );

  const sim_results& res =
    simrun_master(
      par, mpicomm
    );

  schedmsg = SCHEDMSG_EXIT;
  mpi::broadcast( mpicomm, schedmsg, 0 );
}
