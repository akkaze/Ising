#include <iostream>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/filesystem.hpp>
#include <boost/timer/timer.hpp>
#include "sched.hpp"

using namespace std;
namespace mpi = boost::mpi;
namespace fs  = boost::filesystem;

int main(int argc,char** argv)
{
    string file = "stat.txt";

    mpi::environment  mpienv( argc, argv );
    mpi::communicator mpicomm;

    sim_parameters par;

    par.N = atoi( argv[1] );
    par.periodic = atoi( argv[2] );
    par.init = *argv[3];
    par.drysweeps = atoi( argv[4] );

    par.bins = atoi( argv[5] );
    par.binwidth = atoi( argv[6] );
    par.sweeps = atoi( argv[7] );


    par.use_fsize_correction = atoi( argv[8] );
    par.T = atof( argv[9] );
    par.J = atof( argv[10] );
    par.B = atof( argv[11] );

    if ( mpicomm.rank() == 0 ) { // only master prints welcome message
        cout << endl;
        cout << "    ==========================================" << endl;
        cout << "    | hVMC - hubbard Variational Monte Carlo |" << endl;
        cout << "    ==========================================" << endl;
        cout << endl;
    }

    if ( mpicomm.rank() == 0 ) { // this process is the master
        boost::timer::cpu_timer timer;
        sched_master( par, mpicomm );
        cout << timer.format(6) << endl;
    } else { // this process is a slave
        sched_slave( par, mpicomm );
    }

    if ( mpicomm.rank() == 0 ) { // only master prints exit message
        cout << endl;
        cout << ":: All done, exiting!" << endl;
        cout << endl;
    }

    //    if ( !res.success ) {
    //        cout << "WARNING: simulation terminated unsuccessfully!\n";
    //    }

    //    cout << res.B << endl;
    //    cout << "... done!";
    return 0;
}

