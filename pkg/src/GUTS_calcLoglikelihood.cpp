/**
 * GUTS Method calcLoglikelihood.
 *
 * @file    GUTS_calcLoglikelihood.cpp
 * @author  soeren.vogel@posteo.ch
 * @author  carlo.albert@eawag.ch
 * @date    2015-04-01
 * @license GPL-2
 *
 * Calculate the logarithm of the likelihood of data present in GUTS object.
 */





#ifndef guts_h
#include "GUTS.h"
#endif





using namespace std;





/*
 * Method: calcLoglikelihood
 * set: mLL (the logLikelihood or NaN)
 * Messages read: GER_Y, GER_S
 */
void GUTS::calcLoglikelihood() {

	/*
	 * We need survival probabilities (mS) first.
	 * calcSurvProb will also call calcSample().
	 */
	calcSurvivalProbabilities();

	/*
	 * Check for messages, if, return.
	 */
	int er[] = {GER_Y, GER_S};
	for ( unsigned int i = 0; i < 2; ++i ) {
		if ( mMessages.at( er[i] ) ) {
			mLL = GNAN_DOUBLE; // Error.
			return;
		}
	}
	if ( mMessages.at(GWA_PAR) ) {
		mLL = -GINF_DOUBLE;
		return;
	}

	/*
	 * We have not returned so far: continue calculation.
	 * Reset mLL (loglikelihood).
	 * Set mLL to -INF or calc mLL.
	 */
	mLL = 0.0;
	for ( unsigned int i=0; i < myt.size(); ++i ) {
		if ( diffS.at(i) <= 0.0 && diffy.at(i) != 0 ) {
			mLL   =   -GINF_DOUBLE;
			return;
		} else {
			mLL   +=  diffy.at(i) * log( diffS.at(i) );
		}
	}

} // End of GUTS::calcLoglikelihood()
