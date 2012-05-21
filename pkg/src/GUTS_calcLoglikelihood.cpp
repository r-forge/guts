/** GUTS method logLikelihood
 * @file        GUTS_logLikelihood.cpp
 * @author      soeren.vogel@uzh.ch
 * @author      carlo.albert@eawag.ch
 * @date        2012-05-07
 * @license     GPL-2
 *
 * Calculate the logarithm of the likelihood of data present in GUTS object.
 */

#ifndef guts_h
#include "GUTS.h"
#endif

using namespace std;

/*
 * set: mLL
 * error read: GER_Y, GER_S
 * error write: GER_LL
 */

double GUTS::calcLoglikelihood()
{
    /*
     * We need survival probabilities first
     * Will also call calcSample
     */
    calcSurvivalProbabilities();

	// Reset error on this method
    mErrors.at( GER_LL ) = false;
    mErrorMessages.at( GER_LL )	= "";

	/*
     * Check for errors
     */
    int i[] = { GER_Y, GER_S };
    for ( unsigned int j = 0; j < 2; ++j )
    {
        if ( mErrors.at( i[j] ) )
        {
            mErrors.at( GER_LL )		= true;
            mErrorMessages.at( GER_LL )	= "calcLoglikelihood failed: " + mErrorMessages.at( i[j] );
            mLL = GNAN_DOUBLE;
            return mLL;
        }
    }

	/*
     * We have not returned so far: continue calculation
     */
    
    /*
     * Reset mLL (loglikelihood).
     * Append 0 to mS, but don't forget to remove it at the end!
     * Append 0 to my, but don't forget to remove it at the end!
     */
    mLL = 0.0;
    mS.push_back( 0.0 );
    my.push_back( 0 );

	for ( unsigned int i=0; i < myt.size(); ++i )
    {
    	double diffS = mS.at(i) - mS.at(i+1);
        int diffy = my.at(i) - my.at(i+1);

		if ( diffS <= 0.0f && diffy != 0 )
        {
        	mLL = -INFINITY;
            return mLL;
        }
        else
        {
            mLL += diffy * log( diffS );
        }
    }

	/*
     * reset S and y to their original length
     */
    mS.pop_back();
    my.pop_back();

    return mLL;

} // end GUTS::calcLoglikelihood()