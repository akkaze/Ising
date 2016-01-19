/*
 * Copyright (c) 2012, Robert Rueger <rueger@itp.uni-frankfurt.de>
 *
 * This file is part of SSMC.
 *
 * SSMC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SSMC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SSMC.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "simrun.hpp"
#include <fstream>
#include <memory>
#include <vector>
#include <functional>
#include <numeric>

#include <boost/mpi/collectives.hpp>

#include "msgtag.hpp"

// ----- SIMRUN: PERFORMS A SINGLE SIMULATION -----


class Replica
{
public:
    Replica(const sim_parameters& par) {
        model1 = new IsingModel2d( par.N, par.periodic, par.J, par.B, par.T,par.use_fsize_correction );
        model2 = new IsingModel2d( par.N, par.periodic, par.J, par.B, par.T,par.use_fsize_correction );
    }

    virtual ~Replica() {
        delete model1;
        delete model2;
    }

    bool prepare(char init) {
        bool success = (model1->prepare(init) &&
                        model2->prepare(init));
        return success;
    }

    void mcs()
    {
        model1->mcstep();
        model2->mcstep();
    }

    double Q()
    {
        long Ql = 0;
        for ( unsigned int line = 0; line < model1->size; line++ ) {
            for ( unsigned int col = 0; col < model1->size; col++ ) {
                Ql += model1->spin[line][col].get() * model2->spin[line][col].get();
            }
        }
//        cout << Ql << endl;
        double Q = 0;
        Q = static_cast<double>(Ql) / model1->N;
        return Q;
    }

    void mcstep_dry(const unsigned int& dry_sweeps) {
        model1->mcstep_dry(dry_sweeps);
        model2->mcstep_dry(dry_sweeps);
        return;
    }

    IsingModel2d* model1;
    IsingModel2d* model2;
};

sim_results simrun_master( const sim_parameters& par,const mpi::communicator& mpicomm)
{
    // ----- PREPARE SIMULATION -----
    // assume something went wrong until we are sure it didn't
    sim_results res;
    res.success = false;
    // ----- RUN SIMULATION -----
    Replica* rep = new Replica(par);

    if ( rep->prepare( par.init ) == false ) {
        delete rep;
        return res;
    }

    unsigned int finished_workers = 0;
    unsigned int scheduled_bins = 0;
    unsigned int completed_bins = 0;
    unsigned int enqueued_bins  = par.bins;

    // define procedure to query the slaves for new work requests
    function<void()> mpiquery_work_requests( [&]() {
        while ( boost::optional<mpi::status> status
                = mpicomm.iprobe( mpi::any_source, MSGTAG_S_M_REQUEST_BINS ) ) {
            // receive the request and hand out new bins to the source
            mpicomm.recv( status->source(), MSGTAG_S_M_REQUEST_BINS );
            if ( enqueued_bins > 0 ) {
                mpicomm.send( status->source(), MSGTAG_M_S_DISPATCHED_BINS, 1 );
                scheduled_bins += 1;
                enqueued_bins  -= 1;
            } else {
                mpicomm.send( status->source(), MSGTAG_M_S_DISPATCHED_BINS, 0 );
                ++finished_workers;
            }
        }
    } );

    // define procedure to query the slaves for finished work
    function<void()> mpiquery_finished_work( [&]() {
        while ( boost::optional<mpi::status> status
                = mpicomm.iprobe( mpi::any_source, 2 ) ) {
            mpicomm.recv( status->source(), 2 );
            --scheduled_bins;
            ++completed_bins;
        }
    } );

    cout << ":: Performing Monte Carlo cycle" << endl;
    cout << endl;
    cout << "   Progress:" << endl;

    // perform dry runs to reach thermal equilibrium
    for(unsigned int mcs = 0; mcs < par.drysweeps; mcs++) {
        // take care of the slaves
        mpiquery_finished_work();
        mpiquery_work_requests();
        rep->mcs();
    }

    unsigned int completed_bins_master = 0;

    std::vector<double> q2_binmeans;
    std::vector<double> q4_binmeans;

    while ( enqueued_bins > 0 ) {
        cout << '\r' << "     Bin "
             << completed_bins << "/" << par.bins;

        cout.flush();

        --enqueued_bins;
        ++scheduled_bins;
        // initialize binning array
        vector<double> q2_currentbin;
        vector<double> q4_currentbin;
        try {
            // try to allocate enough memory ...
            q2_currentbin.reserve( par.binwidth );
            q4_currentbin.reserve( par.binwidth );
        } catch ( bad_alloc ) {
            delete rep;
            return res;
        }
        for (unsigned int mcs = 0;mcs < par.binwidth;++mcs ) {
            // take care of the slaves
            mpiquery_finished_work();
            mpiquery_work_requests();

            // perform a Monte Carlo step
            rep->mcs();

            // measure observables
            double q2 = 0, q4 = 0;
            double thissample_q = rep->Q();
            // remember the sample's properties to calculate their mean value
            q2 	= thissample_q * thissample_q;
            q4 	= thissample_q * thissample_q * thissample_q * thissample_q;
            q2_currentbin.push_back(q2);
            q4_currentbin.push_back(q4);
        }


        q2_binmeans.push_back(
                    accumulate( q2_currentbin.begin(), q2_currentbin.end(), 0.0 ) /
                    static_cast<double>( q2_currentbin.size() )
                    );
        q2_currentbin.clear();

        --scheduled_bins;
        ++completed_bins_master;
        ++completed_bins;
    }
    ++finished_workers;

    while ( completed_bins != par.bins ||
            static_cast<int>( finished_workers ) < mpicomm.size() ) {
        if ( boost::optional<mpi::status> status
             = mpicomm.iprobe( mpi::any_source, MSGTAG_S_M_FINISHED_BINS ) ) {
            mpicomm.recv( status->source(), MSGTAG_S_M_FINISHED_BINS );
            --scheduled_bins;
            ++completed_bins;

            cout << "\n";
            cout << '\r' << "     Bin " << completed_bins << "/" << par.bins;
            cout.flush();
        }

        if ( boost::optional<mpi::status> status
             = mpicomm.iprobe( mpi::any_source, MSGTAG_S_M_REQUEST_BINS ) ) {
            // receive the request for more work
            mpicomm.recv( status->source(), MSGTAG_S_M_REQUEST_BINS );
            // tell him there is no more work
            mpicomm.send( status->source(), MSGTAG_M_S_DISPATCHED_BINS, 0 );
            ++finished_workers;
        }
    }

    assert( enqueued_bins == 0 );
    assert( scheduled_bins == 0 );

    cout << '\r' << "     Bin " << completed_bins << "/" << par.bins << endl;
    cout.flush();

    // all measurements done ... let's tidy things up
    delete rep;

    assert( mpicomm.rank() == 0 );
    vector< vector<double> > q2_binmeans_collector;
    mpi::gather( mpicomm, q2_binmeans, q2_binmeans_collector, 0 );

    vector<double> q2_binmeans_all;
    for ( auto it = q2_binmeans_collector.begin();
          it != q2_binmeans_collector.end();
          ++it ) {
        q2_binmeans_all.insert( q2_binmeans_all.end(), it->begin(), it->end() );
    }

    double q2 = 0, q4 = 0;
    q2 = static_cast<double>(
                accumulate( q2_binmeans_all.begin(), q2_binmeans_all.end(), 0.0 )
                ) / static_cast<double>( q2_binmeans_all.size() );

    double B = 0;
    B = (3 - q4 / (q2 * q2)) / 2;
    res.B = B;
    res.success = true;
    return res;
}


void simrun_slave( const sim_parameters& par,const mpi::communicator& mpicomm)
{
    Replica* rep = new Replica(par);

    if ( rep->prepare( par.init ) == false ) {
        delete rep;
    }

    // perform dry runs to reach thermal equilibrium
    rep->mcstep_dry( par.drysweeps );

    unsigned int completed_bins_thisslave = 0;
    bool master_out_of_work = false;
    unsigned int scheduled_bins_thisslave;
    mpicomm.send( 0, MSGTAG_S_M_REQUEST_BINS );
    mpicomm.recv( 0, MSGTAG_M_S_DISPATCHED_BINS, scheduled_bins_thisslave );
    master_out_of_work = ( scheduled_bins_thisslave == 0 );

    std::vector<double> q2_binmeans;
    std::vector<double> q4_binmeans;

    while ( scheduled_bins_thisslave > 0 ) {

        unsigned int new_scheduled_bins_thisslave;
        mpi::request master_answer;

        if ( !master_out_of_work ) {
            // ask the master for more work
            mpicomm.send( 0, MSGTAG_S_M_REQUEST_BINS );
            master_answer = mpicomm.irecv(
                        0, MSGTAG_M_S_DISPATCHED_BINS,
                        new_scheduled_bins_thisslave
                        );
        }
        // initialize binning array
        vector<double> q2_currentbin;
        vector<double> q4_currentbin;
        try {
            // try to allocate enough memory ...
            q2_currentbin.reserve( par.binwidth );
            q4_currentbin.reserve( par.binwidth );
        } catch ( bad_alloc ) {
            delete rep;
        }
        for (unsigned int mcs = 0;mcs < par.binwidth;++mcs ) {
            // perform a Monte Carlo step
            rep->mcs();

            // measure observables
            double q2 = 0, q4 = 0;
            double thissample_q = rep->Q();
            // remember the sample's properties to calculate their mean value
            q2 	= thissample_q * thissample_q;
            q4 	= thissample_q * thissample_q * thissample_q * thissample_q;
            q2_currentbin.push_back(q2);
            q4_currentbin.push_back(q4);
        }


        q2_binmeans.push_back(
                    accumulate( q2_currentbin.begin(), q2_currentbin.end(), 0.0 ) /
                    static_cast<double>( q2_currentbin.size() )
                    );
        q2_currentbin.clear();

        // report completion of the work
        mpicomm.send( 0, 2 );
        ++completed_bins_thisslave;
        --scheduled_bins_thisslave;

        if ( !master_out_of_work ) {
            // wait for answer from master concerning the next bin
            master_answer.wait();
            if ( new_scheduled_bins_thisslave == 1 ) {
                ++scheduled_bins_thisslave;
            } else {
                master_out_of_work = true;
            }
        }
    }

    assert( mpicomm.rank() != 0 );
    mpi::gather( mpicomm, q2_binmeans, 0 );
    return;
}
