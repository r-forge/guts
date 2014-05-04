/** GUTS Method calcSample
 *
 * @file    GUTS_calcSample.cpp
 * @author  soeren.vogel@posteo.ch
 * @author  carlo.albert@eawag.ch
 * @date    2012-04-30
 * @license GPL-2
 */

#ifndef guts_h
#include "GUTS.h"
#endif

using namespace std;

/*
 * calcSample()
 * error read: GER_PAR, GER_N, GER_DIST
 * error write: GER_Z
 */
void GUTS::calcSample()
{
  /*
     * Errors checked elsewhere
     */
    int i[] = { GER_PAR, GER_N, GER_DIST };
    for ( unsigned int j = 0; j < 3; ++j )
    {
      if ( mErrors.at( i[j] ) )
        {
          mErrors.at( GER_Z )     = true;
          mErrorMessages.at( GER_Z )   = "calculation of z failed: " + mErrorMessages.at( i[j] );
            mz.resize( 0 );
            return;
        }
    }
  
  /*
     * par, N and dist seem to be okay
     * Run specific checks and assignment
     */
    if ( mdist == "lognormal" )
    {
      if ( mpar.size() < 5 )
        {
          mErrors.at( GER_Z )      = true;
            mErrorMessages.at( GER_Z )   = "error in calcSample(), par needs 5 parameters";
        }
        else if ( mpar.at(3) <= 0.0f || mpar.at(4) <= 0.0f )
        {
          mErrors.at( GER_Z )      = true;
            mErrorMessages.at( GER_Z )   = "error in calcSample(), par 4 and 5 must be non-negative";
        }
        else
        {
          doCalcSampleLognormal( mpar.at(3), mpar.at(4), mN );
          mErrors.at( GER_Z )      = false;
            mErrorMessages.at( GER_Z )   = "";
            return;
        }
    }
} // end GUTS::doSample()

/*
 * Sample lognormal
 */
void GUTS::doCalcSampleLognormal( const double& mean, const double& sigma, const int& n )
{
    using namespace boost;

    // Mersenne twister RNG, seeded once with #seconds since 1970
    static boost::random::mt19937 rng(static_cast<unsigned> (std::time(0)));

    // Select Gaussian probability distribution with mean and sigma
    boost::lognormal_distribution<double> lognorm_dist( mean, sigma );

    // Bind RNG to distribution, forming a function
    variate_generator<boost::random::mt19937&, boost::lognormal_distribution<double> >
        lognormal_sampler(rng, lognorm_dist);

    // Sample from the distribution
    mz.resize( 0 );
    for ( int i = 0; i < n; ++i )
    {
        mz.push_back( lognormal_sampler() );
    }
} // end GUTS::doSampleLognormal()
