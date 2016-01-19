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

#include "model_ising2dsqrmet.hpp"


// ----- 2D ISING-MODEL ON A SQUARE LATTICE (METROPOLIS-ALGORITHM)

IsingModel2d::IsingModel2d( const unsigned int& N, const bool& periodic,
                            const double& J, const double& B, const double& T,
                            const unsigned int& fsize_correction_mode)
  : periodic( periodic ),J( J ),B( B ),N( N* N ),T( T ),size( N ),
    fsize_correction_mode( fsize_correction_mode ),
    fsize_ordered_phase( false )
{
  // set up two dimensional NxN array: spin[line][column]
  spin.resize( size );
  for ( unsigned int line = 0; line < N; line++ ) {
    spin[line].resize( size );
  }

  rng = gsl_rng_alloc( gsl_rng_mt19937 );
  gsl_rng_set( rng, rand() );
}

IsingModel2d::~IsingModel2d()
{
    // free rng's memory
    gsl_rng_free( rng );
}


unsigned int IsingModel2d::spin_count() const
{
  // returns the total number of spins in the system
  return N;
}


bool IsingModel2d::prepare( const char& mode )
{
  switch ( mode ) {

  case 'r':	// completely random state ...
    for ( unsigned int line = 0; line < size; line++ ) {
      for ( unsigned int col = 0; col < size; col++ ) {
        if ( gsl_rng_uniform_int( rng, 2 ) == 0 ) {
          spin[line][col].flip();
        }
      }
    }
    break;

  case 'u': // sets all spins up
    for ( unsigned int line = 0; line < size; line++ ) {
      for ( unsigned int col = 0; col < size; col++ ) {
        spin[line][col].set( +1 );
      }
    }
    break;

  case 'd': // sets all spins down
    for ( unsigned int line = 0; line < size; line++ ) {
      for ( unsigned int col = 0; col < size; col++ ) {
        spin[line][col].set( -1 );
      }
    }
    break;
  }
  return true;
}


void IsingModel2d::metropolis_singleflip()
{
  // find a random spin to flip
  unsigned int flip_line = gsl_rng_uniform_int( rng, size );
  unsigned int flip_col  = gsl_rng_uniform_int( rng, size );

  // flip it!
  spin[flip_line][flip_col].flip();

  // calculate energy difference
  double deltaH = - 2 * B * spin[flip_line][flip_col].get();
  if ( periodic ) {
    // top neighbour
    deltaH += -2 * J * ( spin[flip_line][flip_col]
                         * spin[( flip_line + size - 1 ) % size][flip_col] );
    // bottom neighbour
    deltaH += -2 * J * ( spin[flip_line][flip_col]
                         * spin[( flip_line + 1 ) % size][flip_col] );
    // left neighbour
    deltaH += -2 * J * ( spin[flip_line][flip_col]
                         * spin[flip_line][( flip_col + size - 1 ) % size] );
    // right neighbour
    deltaH += -2 * J * ( spin[flip_line][flip_col]
                         * spin[flip_line][( flip_col + 1 ) % size] );
  } else {
    if ( flip_line != 0 ) {
      // top neighbour
      deltaH += -2 * J * ( spin[flip_line][flip_col]
                           * spin[flip_line - 1][flip_col] );
    }
    if ( flip_line != size - 1 ) {
      // bottom neighbour
      deltaH += -2 * J * ( spin[flip_line][flip_col]
                           * spin[flip_line + 1][flip_col] );
    }
    if ( flip_col != 0 ) {
      // left neighbour
      deltaH += -2 * J * ( spin[flip_line][flip_col]
                           * spin[flip_line][flip_col - 1] );
    }
    if ( flip_col != size - 1 ) {
      // right neighbour
      deltaH += -2 * J * ( spin[flip_line][flip_col]
                           * spin[flip_line][flip_col + 1] );
    }
  }

  if ( deltaH > 0 ) {
    // read or calculate exp(- deltaH / T)
    double exp_deltaHoverT;
    if ( exp_precalc.count( deltaH ) == 1 ) {
      exp_deltaHoverT = exp_precalc[deltaH];
    } else {
      exp_deltaHoverT = exp( - deltaH / T );
      exp_precalc[deltaH] = exp_deltaHoverT;
    }
    // accept the new state?
    if ( gsl_rng_uniform( rng ) > exp_deltaHoverT ) {
      // new state rejected ... reverting!
      spin[flip_line][flip_col].flip();
      return;
    }
  }
}


void IsingModel2d::mcstep()
{
  // empty precalculated exponentials (needed for SA)
  exp_precalc.clear();

  for ( unsigned long int n = 1; n <= N; n++ ) {
    metropolis_singleflip();
  }
}


void IsingModel2d::mcstep_dry( const unsigned int& k_max )
{
  for ( unsigned int k = 0; k < k_max; k++ ) {
    mcstep();
  }
}


double IsingModel2d::H() const
{
  // measures the system's energy

  double H = 0;

  // energy in external magnetic field
  if ( B != 0 ) {
    for ( unsigned int line = 0; line < size; line++ ) {
      for ( unsigned int col = 0; col < size; col++ ) {
        H += - B * spin[line][col].get();
      }
    }
  }

  // energy due to interaction within the lattice
  for ( unsigned int line = 0; line < size - 1; line++ ) {
    for ( unsigned int col = 0; col < size - 1; col++ ) {
      H += - J * ( spin[line][col] * spin[line + 1][col] ); // below
      H += - J * ( spin[line][col] * spin[line][col + 1] ); // right
    }
  }
  for ( unsigned int col = 0; col < size - 1; col++ ) {
    // horizontal neighbours in the last line
    H += - J * ( spin[size - 1][col] * spin[size - 1][col + 1] );
  }
  for ( unsigned int line = 0; line < size - 1; line++ ) {
    // vertical neighbours in the last column
    H += - J * ( spin[line][size - 1] * spin[line + 1][size - 1] );
  }

  if ( periodic ) {
    // interaction over vertical and horizontal borders
    for ( unsigned int col = 0; col < size; col++ ) {
      H += - J * ( spin[0][col] * spin[size - 1][col] );
    }
    for ( unsigned int line = 0; line < size; line++ ) {
      H += - J * ( spin[line][0] * spin[line][size - 1] );
    }
  }
  return H;
}


double IsingModel2d::h() const
{
  // measures the system's energy per spin
  return H() / N;
}


int IsingModel2d::M() const
{
  // measures the system's magnetization
  int M = 0;
  for ( unsigned int line = 0; line < size; line++ ) {
    for ( unsigned int col = 0; col < size; col++ ) {
      M += spin[line][col].get();
    }
  }

  // finite size corrections to the magnetization
  if ( ( fsize_correction_mode == 1 ) ||
       ( ( fsize_correction_mode == 2 ) && fsize_ordered_phase ) ) {
    M = abs( M );
  }
  return M;
}


double IsingModel2d::m() const
{
  // measures the system's magnetization per spin
  return double( M() ) / N;
}

vector<double> IsingModel2d::ss_corr() const
{
  // measure spin-spin correlations
  vector<double> result;
  vector<unsigned int> samples;
  result.resize( spin.size(), 0 );
  samples.resize( spin.size(), 0 );
  for ( int i = 0; i < size; i++ ) {
    for ( int j = 0; j < size; j++ ) {
      result[abs( i - j ) ] += spin[i][i] * spin[i][j];
      result[abs( i - j ) ] += spin[i][i] * spin[j][i];
      samples[abs( i - j ) ] += 2;
    }
  }
  for ( unsigned int d = 0; d < size; d++ ) {
    result[d] /= samples[d];
  }
  return result;
}
