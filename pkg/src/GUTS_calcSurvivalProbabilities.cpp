/** GUTS method calcSurvivalProbabilities
 * @file		GUTS_calcSurvivalProbabilities.cpp
 * @author      soeren.vogel@uzh.ch
 * @author      carlo.albert@eawag.ch
 * @date        2012-04-30
 * @license     GPL-2
 */

#ifndef guts_h
#include "GUTS.h"
#endif

using namespace std;

/*
 * Method calcSurvivalProbabilities
 * Overloaded!
 * Call blank to use my- and myt-Values
 * Call with vector, use this as myt-Values
 * Call with integer, use this as the length of a vector
 * error read: GER_C, GER_CT, GER_YT, GER_PAR, GER_M, GER_Z
 * error write: GER_S
 */


/*
 * standard implementation
 * with no argument = rely on GUTS objects fields
 * all overloaded methods will refer to this one here
 */
void GUTS::calcSurvivalProbabilities()
{
    // Set mz
    calcSample();

	// Reset error on this method
    mErrors.at( GER_S ) = false;

    // Check for errors
    int i[] = { GER_C, GER_CT, GER_YT, GER_PAR, GER_M, GER_Z };
    for ( unsigned int j = 0; j < 6; ++j )
    {
    	if ( mErrors.at( i[j] ) )
        {
        	mErrors.at( GER_S )			= true;
            mErrorMessages.at( GER_S )	= "calculation of survival probabilities failed, check conditions";
	        mS.resize( 0 );
	        return;
        }
    }

	/*
     * We have not returned, and
     * hence check errors not checked so far
     */
    if ( myt.size() <= 0 || mCt.size() <= 0 || myt.back() > mCt.back() )
    {
        mErrors.at( GER_S )			= true;
        mErrorMessages.at( GER_S )	= "calculation of survival probabilities failed, check yt and Ct";
        mS.resize( 0 );
        return;
    }

	// If no new errors
	doCalcSurvivalProbabilities();

} // end GUTS::calcSurvivalProbabilities()


/*
 * St: a vector of survival time points
 */
void GUTS::calcSurvivalProbabilities( const vector<double>& St )
{
    vector<int> hy;
    hy.assign( St.size(), 0 );
    setSurvivors( hy, St);
	// refer to standard method above
    calcSurvivalProbabilities();
} // end GUTS::calcSurvivalProbabilities( const vector<double>& St )


/*
 * Stlength: the number of survivor time points
 */
void GUTS::calcSurvivalProbabilities( const unsigned int& Stlength )
{
    vector<int> hy;
    hy.assign( Stlength, 0 );
    vector<double> hyt;
    for ( unsigned int i=0; i<Stlength; ++i )
        hyt.push_back( i );
    setSurvivors( hy, hyt );
    // refer to standard method above
    calcSurvivalProbabilities();
} // end GUTS::calcSurvivalProbabilities( const int& Stlength )


/*
 * the work horse
 * all errors should have been checked
 */
void GUTS::doCalcSurvivalProbabilities()
{
	/*
     * We re-calculate mS
     */
    mS.resize( 0 );

    /*
     * Method internal vectors and variables
     */
    vector<double> D;
    D.assign( mM, 0.0 );
    vector<double> ee;
    ee.assign( mN, 0.0 );
    vector<int> ff;
    ff.assign( mN, 0 );
    double dtau = ( myt.back() - myt.front() ) / mM;

    /*
     * we need to copy z because we must sort
     */
    vector<double> mmz = mz;
    sort( mmz.begin(), mmz.end() );

    /*
     * Iterators:
     *
     * j: iterator for D, ends at mM
     * ii: iterator for D, is the last j
     * zpos: current position in z
     * tau: 0 -- max(yt), 10000 elements, incremented by value of dtau
     * k: index for C
     * i: iterator for yt
     */
    int j = 0;
    int ii = 0;
    int zpos = 0;
    double tau = 0.0;
    double dk = mpar.at(2) * dtau; // killing rate
    int k = 0;

    /*
     * Loop over yt_i's.
     */
    for ( unsigned int i = 0; i < myt.size(); ++i )
    {
        /*
         * Must control for tau.
         * Must control for j to not exceed dim
         * of D due to calculation inaccuracy.
         */
        while ( tau < myt.at(i) && j < mM )
        {
            /*
             * Calculate D
             */
            double tmp = exp( -mpar.at(1) * ( tau - mCt.at(k) ) );

            double carlosKlammer;
            if ( mpar.at(1) <= 0.0f )
                carlosKlammer = tau - mCt.at(k);
            else
                carlosKlammer = ( 1 / mpar.at(1) ) * ( 1-tmp );

            D.at(j) = D.at(ii) * tmp
                    + mC.at(k) * ( 1-tmp )
            		+ ( mC.at(k+1) - mC.at(k) )
            		* ( tau - mCt.at(k) - carlosKlammer )
            		/ ( mCt.at(k+1) - mCt.at(k) );
            
            /*
             * Increment or decrement zpos.
             */
            if ( D.at(j) < mmz.at(zpos) )
            {
                while ( D.at(j) < mmz.at(zpos) && zpos > 0 )
                {
                    --zpos;
                }
                
                // zpos may have changed, so mmz.at(zpos) may also.
                // D.at(j) can now be larger than mmz.at(zpos).
                if ( D.at(j) > mmz.at(zpos) )
                {
                    ++zpos;
                }
            }
            else
            {
                while ( D.at(j) > mmz.at(zpos) && zpos < (mN-1) )
                {
                    ++zpos;
                }
            }
            
            /*
             * Update ee and ff.
             */
            if ( D.at(j) > mmz.at(mN-1) ) //
            {
                ee.at(mN-1) += D.at(j);
                ff.at(mN-1)++;
            }
            else if ( D.at(j) > mmz.at(0) )
            {
                ee.at(zpos-1) += D.at(j);
                ff.at(zpos-1)++;
            }
            
            // Increment or decrement j, tau, k
            j++;
            tau += dtau;
            if ( tau > mCt.at(k+1) )
            {
                k++;
                ii = j-1;
            }
        } // end while ( tau < myt.at(i) )

        // Calculate S.at(i):
        // E = sum of all ee
        // F = sum of all ff
/*
        double E = accumulate( ee.begin(), ee.end(), 0.0);
        double F = accumulate( ff.begin(), ff.end(), 0.0);
        // Calculate mmS.at(i):
        mmS.at(i) = exp( dk * ( mmz.at(0) * F - E ) );
        for ( int u=1; u < mN; ++u )
        {
            E -= ee.at(u-1);
            F -= ff.at(u-1);
            mmS.at(i) += exp( dk * ( mmz.at(u) * F - E ) );
        }
*/
/**/    // new
        double F = ff.back();
        double E = ee.back();

        /*
         * write mS.at( i )
         * i is the loop through all yt elements
         * hence mS appends until the length of yt
         */
        mS.push_back( exp( dk * ( mmz.back() * F - E ) ) );
//        int tmpa = mmz.size();//debug

        for ( unsigned int u=mN-1; u > 0; --u )
        {
            E += ee.at(u-1);
            F += ff.at(u-1);
            mS.at(i) += exp( dk * ( mmz.at(u-1) * F - E ) );
        }
/**/

        mS.at(i) *= exp( -mpar.front() * myt.at(i) ) / mN;
        
    } // end for ( unsigned int i = 0; i < myt.size (); i++ )

} // end GUTS::calcSurvivalProbabilities()
