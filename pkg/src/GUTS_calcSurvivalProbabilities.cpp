/**
 * GUTS Method calcSurvivalProbabilities
 *
 * @file    GUTS_calcSurvivalProbabilities.cpp
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
 * Method calcSurvivalProbabilities
 * Overloaded!
 * Call blank to use values in my and myt.
 * Call with vector, use this as myt-Values.
 * Call with integer, use this as the length of a vector.
 * Messages read: GER_C, GER_CT, GER_YT, GER_PAR, GER_M, GER_Z
 * Waring read: GWA_PAR, if return.
 * Messages write: GER_S
 */





/**
 * Standard implementation.
 * @param   <none>
 * Rely on GUTS objects fields.
 * Note: All overloaded methods will refer to this one here.
 */
void GUTS::calcSurvivalProbabilities() {

	/*
	 * Do sample.
	 * Set mz, the sample vector,
	 * and mzw, the weight vector.
	 */
	calcSample();

	/*
	 * Messages checked in calcSample: { GER_PAR, GER_MODEL, GER_DIST, GER_N, GER_EMP } --> GER_Z
	 * Errors set in setSample: GER_EMP, also checked in calcSample.
	 * Errors needed here: { GER_Z, GER_YT, GER_C, GER_CT, GER_M }
	 */
	int er[] = { GER_Z, GER_YT, GER_C, GER_CT, GER_M };
	for ( unsigned int i = 0; i < 5; ++i ) {
		if ( mMessages.at( er[i] ) ) {
			mMessages.at(GER_S)  = true; // "Calculation of survival probabilities failed.";
			return;
		}
	}

	/*
	 * Return if par warning.
	 */
	if ( mMessages.at(GWA_PAR) ) return;

	/*
	 * We have not returned so far.
	 * Reset error on this method.
	 */
	mMessages.at(GER_S)  = false;

	/*
	 * We have not returned, and hence there is no errors, so far. Check now:
	 * (1) Both time vectors must be set up properly
	 * (2) The last time point of concentrations must be larger or equal
	 *     then the last time point of survivors: there must be information
	 *     about concentrations when there is no more information about
	 *     survivors.
	 */
//	if ( myt.size() <= 0 || mCt.size() <= 0 || myt.back() > mCt.back() ) {
	if ( mCt.back() < myt.back() ) {
		mMessages.at(GER_S)  = true; // "Calculation of survival probabilities failed, check yt and Ct.";
		return;
	}

	/*
	 * No new errros, calculate survival probabilities.
	 */
	doCalcSurvivalProbabilities();

} // End of GUTS::calcSurvivalProbabilities()





/**
 * With-vector implementation.
 * @param   St  vector of survival time points.
 */
void GUTS::calcSurvivalProbabilities( const vector<double>& St ) {
	
	/* Create a 0-vector of survivors. */
	vector<int> hy;
	hy.assign(St.size(), 0);
	
	/* Use the standard setter. */
	setSurvivors(hy, St);
	
	/* Employ the standard method above. */
	calcSurvivalProbabilities();
	
} // End of GUTS::calcSurvivalProbabilities( … )





/**
 * Vector-length implementation.
 * @param   Stlength    number of survivor time points
 */
void GUTS::calcSurvivalProbabilities( const unsigned int& Stlength ) {
	
	/* Create a 0-vector of survivors. */
	vector<int> hy;
	hy.assign( Stlength, 0 );
	
	/* Create an ascending surivor time points vector. */
	vector<double> hyt;
	for ( unsigned int i=0; i<Stlength; ++i ) {
		hyt.push_back(i);
	}
	
	/* Use the standard setter. */
	setSurvivors(hy, hyt);
	
	/* Employ the standard method above. */
	calcSurvivalProbabilities();
	
} // End of GUTS::calcSurvivalProbabilities( … )





/**
 * Calculate Survival Probabilities – The Work Horse.
 * Each method above calls this work horse.
 * All errors should have been checked!
 */
void GUTS::doCalcSurvivalProbabilities() {

	/*
	 * Method-relevant vectors and variables:
	 * mS:    Vector of survial probabilities, set in setSurvivors.
	 * D:     Vector of damage, set in setTimeGridPoints.
	 * ee:    TODO…, set in setSampleLength and setSample.
	 * ff:    TODO…, set in setSampleLength and setSample.
	 * dtau:  TODO: …, method-internal variable.
	 * E:     last ee
	 * F:     last ff
	 */
	D.assign( mM, 0.0 );
	ee.assign( mN, 0.0 );
	ff.assign( mN, 0 );
	double dtau = (myt.back() - myt.front()) / mM;
	double Scale;
	double tmp;
	double carlosKlammer;
	double E;
	double F;

	/*
	 * Variables:
	 * j: Iterator for D, ends at mM
	 * ii: iterator for D, is the last j
	 * zpos: current position in z
	 * tau: 0 -- max(yt), 10000 elements, incremented by value of dtau
	 * k: index for Ct
	 * i: iterator for yt
	 */
	int j       = 0;
	int ii      = 0;
	int zpos    = 0;
	double tau  = 0.0;
	int k       = 0;

	/*
	 * Loop over yt_i's.
	 */
	for ( unsigned int i = 0; i < myt.size(); ++i ) {

		/*
		 * Must control for tau.
		 * Must control for j to not exceed dim of D due to calculation inaccuracy.
		 */
		while ( tau < myt.at(i) && j < mM ) {

			/*
			 * tmp is needed for…                  TODO: Documentation.
			 */
			tmp = exp( -mpar.at(1) * (tau - mCt.at(k)) );

			/*
			 * Catch exception if the second parameter (k_r) is not greater than 0.
			 */
			if ( mpar.at(1) > 0.0 ) {
				carlosKlammer = (1.0 / mpar.at(1)) * (1.0-tmp);
			} else {
				carlosKlammer = tau - mCt.at(k);
			}

			/*
			 * Calculate D at j.
			 */
			D.at(j) = D.at(ii) * tmp
			        + mC.at(k) * (1.0-tmp) + (mC.at(k+1) - mC.at(k))
			        * (tau - mCt.at(k) - carlosKlammer) / (mCt.at(k+1) - mCt.at(k));

			/*
			 * Increment or decrement zpos.
			 */
			if ( D.at(j) < mz.at(zpos) ) {
				while ( D.at(j) < mz.at(zpos) && zpos > 0 ) {
					--zpos;
				}

				/*
				 * zpos may have changed, so mz.at(zpos) may also.
				 * D.at(j) can now be larger than mz.at(zpos)!
				 */
				if ( D.at(j) > mz.at(zpos) ) {
					++zpos;
				}
			} else {
				while ( D.at(j) > mz.at(zpos) && zpos < (mN-1) ) {
					++zpos;
				}
			}

			/*
			 * Update ee and ff.
			 */
			if ( D.at(j) > mz.at(mN-1) ) {
				ee.at(mN-1) += D.at(j);
				ff.at(mN-1)++;
			} else if ( D.at(j) > mz.at(0) ) {
				ee.at(zpos-1) += D.at(j);
				ff.at(zpos-1)++;
			}

			/*
			 * Increment or decrement j, tau, k
			 */
			j++;
			tau += dtau;
			if ( tau > mCt.at(k+1) ) {
				k++;
				ii = j-1;
			}

		} // End of while ( tau < myt.at(i) )

		/*
		 * Write mS at i.
		 * * E: last ee
		 * * F: last ff
		 * i is the loop through all yt elements
		 * hence mS appends until the length of yt
		 * mpar.at(2) * dtau = killing rate
		 */
		E = ee.back();
		F = ff.back();
		mS.at(i) = exp(   (mpar.at(2) * dtau * (mz.back() * F - E))  +  mzw.back()   );
		for ( unsigned int u=mN-1; u > 0; --u ) {
			E += ee.at(u-1);
			F += ff.at(u-1);
			mS.at(i) += exp(   (mpar.at(2) * dtau * (mz.at(u-1) * F - E))  +  mzw.at(u-1)   );
		}
		mS.at(i) *= exp( -mpar.front() * myt.at(i) ) / mN;

		/*
		 * Scale for S with 1/S_1.
		 * Build diffS from i==2 on. mS.at(i-1) is already scaled!
		 */
		if ( i == 0 ) {
			Scale = 1.0 / mS.at(0);
		} else {
			diffS.at(i-1) = mS.at(i-1) - Scale * mS.at(i);
		}
		mS.at(i) *= Scale;

	} // End of for ( unsigned int i = 0; i < myt.size (); i++ )

	/*
	 * Set the last diffS = mS.back() - 0 = mS.back().
	 */
	diffS.back() = mS.back();

	/*
	 * Notes:
	 *
	 * Multiply vector S with invert of S.at(0).
	 * Old Version.
	 * http://stackoverflow.com/questions/3885095/multiply-vector-elements-by-a-scalar-value-using-stl
	 * std::transform(mS.begin(), mS.end(), mS.begin(), std::bind1st( std::multiplies<double>(), (1.0/mS.at(0)) ));
	 */

} // End of GUTS::calcSurvivalProbabilities()


























