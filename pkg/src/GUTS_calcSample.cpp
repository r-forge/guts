/**
 * GUTS Method calcSample.
 *
 * @file    GUTS_calcSample.cpp
 * @author  soeren.vogel@posteo.ch
 * @author  carlo.albert@eawag.ch
 * @date    2015-04-01
 * @license GPL-2
 */





#ifndef guts_h
#include "GUTS.h"
#endif





using namespace std;





/**
 * calcSample().
 * Messages read: GER_PAR, GER_MODEL, GER_DIST, GER_N, GER_EMP.
 * Return if GWA_PAR is true.
 * Messages write: GER_Z.
 * Set mz as sorted sample.
 * Set mzw (weights) according to sample procedure.
 */
void GUTS::calcSample() {

	/*
	 * Check for messages from other vectors.
	 * If, overwrite GER_Z, even if already set from setSample().
	 */
	int er[] = { GER_PAR, GER_MODEL, GER_DIST, GER_N, GER_EMP };
	for ( unsigned int i = 0; i < 5; ++i ) {
		if ( mMessages.at( er[i] ) ) {
			mMessages.at( GER_Z ) = true; // "Calculation of z failed.";
			return;
		}
	}

	/*
	 * Do not sample if GWA_PAR == true.
	 */
	if ( mMessages.at(GWA_PAR) ) return;

	/*
	 * We have not returned. So reset messages, continue.
	 */
	mMessages.at( GER_Z ) = false;

	/*
	 * par, N, model and dist seem to be okay.
	 * Check for errors in par according to dist etc.
	 */
	if ( mExperiment < 15 ) {
		// lognormal: 11, 12
		if (  wpar.at(4) == 0.0  ) {
			// delta-like
			mz.assign(  mN, wpar.at(3) );
			mzw.assign( mN, 0.0 );
			return;
		} else {
			doSampleGrid();
			return;
		}
	} else if ( mExperiment == 15 ) {
		// lornormold: 15 – testing only
		if (  wpar.at(4) == 0.0  ) {
			// delta-like
			mz.assign(  mN, wpar.at(3) );
			mzw.assign( mN, 0.0 );
			return;
		} else {
			doSampleLognormold();
			return;
		}
	} else if ( mExperiment < 30 ) {
		// delta: 21, 22
		mz.assign(  mN, wpar.at(3) );
		mzw.assign( mN, 0.0 );
		return;
	}

} // End of GUTS::calcSample().





/**
 * doSampleGrid()
 * Grid, speed-up version.
 */
void GUTS::doSampleGrid() {

	/*
	 * Calculate grid intervals. Old Version.
	 */
	//double sigma       =  sqrt(log(   1.0  +  (pow(wpar.at(4), 2.0) / pow(wpar.at(3), 2.0))   ));
	//double mu          =  log(wpar.at(3))  -  (pow(sigma, 2.0) / 2.0);
	//double intervalLo  =  mu - (mInterval * sigma);
	//double intervalHi  =  mu + (mInterval * sigma) + DBL_EPSILON; // DBL_MIN, http://www.cplusplus.com/reference/cfloat/
	//double gridstep    =  (intervalHi - intervalLo) / (mN - 1);
	//mz.resize( mN );
	//mzw.resize( mN );
	//for ( unsigned int i = 0; i < mN; ++i ) {
	//	mz.at(i)    =   exp(   intervalLo  +  (gridstep * i)   );
	//	mzw.at(i)   =   -0.5   *   pow(mu - log(mz.at(i)), 2.0)   /   pow(sigma, 2.0);
	//}

	/*
	 * Make SV’s Beautiful New Grid :-)
	 */
	// R:
	// defaultGrid <- seq( from=-(N-1), to=(N-1), by=2 ) / (N-1)
	// y <- defaultGrid * sigmaD
	// z <- exp( y + mu )
	// w <- -1/2 * y^2 / sigma2
	double sigma2   =  log(   1.0  +  pow( (wpar.at(4) / wpar.at(3)), 2.0 )   );
	double mu       =  log(wpar.at(3))  -  (0.5 * sigma2);
	double sigmaD   =  sqrt(sigma2) * mInterval;
	for ( unsigned int i = 0; i < mN; ++i ) {
		mz.at(i)  =  exp( defaultGrid.at(i) * sigmaD + mu);
		mzw.at(i) =  -0.5  *  pow((defaultGrid.at(i) * sigmaD), 2.0)  /  sigma2;
	}

} // End of GUTS::doSampleGrid().





/**
 * doSampleLognormold( … ).
 * The old actual lognormal sampling engine.
 */
void GUTS::doSampleLognormold() {

	/*
	 * Codes taken from boost examples.
	 */
	using namespace boost;

	/*
	 * Mersenne twister RNG, seeded once with #seconds since 1970.
	 */
	static boost::random::mt19937 rng( static_cast<unsigned> (std::time(0)) );

	/*
	 * Select Gaussian probability distribution with mean and sigma.
	 */
	//boost::lognormal_distribution<double> lognorm_dist( m, sd );
	boost::lognormal_distribution<double> lognorm_dist( wpar.at(3), wpar.at(4) );

	/*
	 * Bind RNG to distribution, forming a function.
	 */
	variate_generator< boost::random::mt19937&, boost::lognormal_distribution<double> > lognormal_sampler(rng, lognorm_dist);

	/*
	 * Sample from the distribution, sort.
	 * Set weights to 0.0 with length of mz.
	 */
	for ( int i = 0; i < mN; ++i ) {
		mz.at(i) = lognormal_sampler();
	}
	sort( mz.begin(), mz.end() );
	mzw.assign( mN, 0.0 );

} // End of GUTS::doSampleLognormold()



