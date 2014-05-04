/** GUTS Method calcSurvivalProbabilities
 *
 * @file    GUTS_calcSurvivalProbabilities.cpp
 * @author  soeren.vogel@posteo.ch
 * @author  carlo.albert@eawag.ch
 * @date    2012-05-14
 * @license GPL-2
 */

#ifndef guts_h
#include "GUTS.h"
#endif

using namespace std;



/** Method calcSurvivalProbabilities
 * Overloaded!
 * Call blank to use my- and myt-Values
 * Call with vector, use this as myt-Values
 * Call with integer, use this as the length of a vector
 * error read: GER_C, GER_CT, GER_YT, GER_PAR, GER_M, GER_Z
 * error write: GER_S
 */



/** Standard implementation.
 * @param   <none>
 * Rely on GUTS objects fields.
 * Note: All overloaded methods will refer to this one here.
 */
void GUTS::calcSurvivalProbabilities() {

  // Set mz, the sample vector
  calcSample();

  // Reset error on this method
  mErrors.at(GER_S)         = false;
  mErrorMessages.at(GER_S)  = "";

  // Check for errors
  int i[] = {GER_C, GER_CT, GER_YT, GER_PAR, GER_M, GER_Z};
  for (unsigned int j = 0; j < 6; ++j) {
    if ( mErrors.at(i[j]) ) {
      mErrors.at(GER_S)         = true;
      mErrorMessages.at(GER_S)  = "calculation of survival probabilities failed, check conditions";
      mS.resize(0);
      return;
    }
  }

  /*
   * We have not returned, and hence check errors not checked so far:
   * (1) Both time vectors must be set up properly
   * (2) The last time point of concentrations must be larger or equal
   *     to the last time point of survivors: there must be information
   *     about concentrations when there is no more information about
   *     survivors.
   */
  if (myt.size() <= 0 || mCt.size() <= 0 || myt.back() > mCt.back()) {
    mErrors.at(GER_S)         = true;
    mErrorMessages.at(GER_S)  = "calculation of survival probabilities failed, check yt and Ct";
    mS.resize(0);
    return;
  }

  // If no new errors
  doCalcSurvivalProbabilities();

} // end GUTS::calcSurvivalProbabilities()



/** With-vector implementation.
 * @param   St  vector of survival time points
 */
void GUTS::calcSurvivalProbabilities(const vector<double>& St) {
  // Create a 0-vector of survivors.
  vector<int> hy;
  hy.assign(St.size(), 0);
  // Use the standard setter.
  setSurvivors(hy, St);
  // Employ the standard method above.
  calcSurvivalProbabilities();
} // end GUTS::calcSurvivalProbabilities( const vector<double>& St )

/** Vector-length implementation.
 * @param   Stlength    number of survivor time points
 */
void GUTS::calcSurvivalProbabilities(const unsigned int& Stlength) {
  // Create a 0-vector of survivors.
  vector<int> hy;
  hy.assign( Stlength, 0 );
  // Create an ascending surivor time points vector.
  vector<double> hyt;
  for (unsigned int i=0; i<Stlength; ++i)
    hyt.push_back(i);
  // Use the standard setter.
  setSurvivors(hy, hyt);
  // Employ the standard method above.
  calcSurvivalProbabilities();
} // end GUTS::calcSurvivalProbabilities( const int& Stlength )



/** The work horse.
 * Each method above calls this work horse.
 * All errors should have been checked!
 */
void GUTS::doCalcSurvivalProbabilities() {

  /*
   * We re-calculate mS. Hence we resize it
   *   to 0 first, and append later.
   */
  mS.resize(0);

  /*
   * mD is a member variable holding the damage.
   */
  mD.resize(mM);
  mD.assign(mM, 0.0);

  /** Method-internal vectors and variables:
   * ee:              TODO: missing documentation
   * ff:              TODO: missing documentation
   * dtau:            TODO: missing documentation
   */
  vector<double> ee;
  ee.assign(mN, 0.0);
  vector<int> ff;
  ff.assign(mN, 0);
  double dtau = (myt.back() - myt.front()) / mM;

  /*
   * We need a copy of mz because we must sort.
   */
  vector<double> mmz = mz;
  sort(mmz.begin(), mmz.end());

  /** Variables:
   * j: Iterator for mD, ends at mM
   * ii: iterator for mD, is the last j
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
  for (unsigned int i = 0; i < myt.size(); ++i) {
    /*
     * Must control for tau.
     * Must control for j to not exceed dim of mD
     * due to calculation inaccuracy.
     */
    while (tau < myt.at(i) && j < mM) {

      /*
       * Calculate mD at j
       */
      double tmp = exp( -mpar.at(1) * (tau - mCt.at(k)) );
      // Catch exception if the second parameter is not greater than 0
      double carlosKlammer;
      if (mpar.at(1) <= 0.0f)
        carlosKlammer = tau - mCt.at(k);
      else
        carlosKlammer = (1 / mpar.at(1)) * (1-tmp);
      mD.at(j) = mD.at(ii) * tmp
                + mC.at(k) * (1-tmp)
                + (mC.at(k+1) - mC.at(k))
                * (tau - mCt.at(k) - carlosKlammer)
                / (mCt.at(k+1) - mCt.at(k));

      /*
       * Increment or decrement zpos.
       */
      if (mD.at(j) < mmz.at(zpos)) {
        while (mD.at(j) < mmz.at(zpos) && zpos > 0) {
          --zpos;
        }

        /*
         * zpos may have changed, so mmz.at(zpos) may also.
         * mD.at(j) can now be larger than mmz.at(zpos)!
         */
        if (mD.at(j) > mmz.at(zpos)) {
          ++zpos;
        }
      } else {
        while (mD.at(j) > mmz.at(zpos) && zpos < (mN-1)) {
          ++zpos;
        }
      }

      /*
       * Update ee and ff.
       */
      if (mD.at(j) > mmz.at(mN-1)) {
        ee.at(mN-1) += mD.at(j);
        ff.at(mN-1)++;
      } else if ( mD.at(j) > mmz.at(0) ) {
        ee.at(zpos-1) += mD.at(j);
        ff.at(zpos-1)++;
      }

      /*
       * Increment or decrement j, tau, k
       */
      j++;
      tau += dtau;
      if (tau > mCt.at(k+1)) {
        k++;
        ii = j-1;
      }
    } // end while ( tau < myt.at(i) )
    
    /** Calculate S.at(i):
     * E: last ee
     * F: last ff
     */
    double E = ee.back();
    double F = ff.back();

    /** Write mS at i
     * i is the loop through all yt elements
     * hence mS appends until the length of yt
     * mpar.at(2) * dtau = killing rate
     */

    mS.push_back( exp( mpar.at(2) * dtau * (mmz.back() * F - E) ) );
    //        int tmpa = mmz.size();//debug

    for (unsigned int u=mN-1; u > 0; --u) {
      E += ee.at(u-1);
      F += ff.at(u-1);
      mS.at(i) += exp( mpar.at(2) * dtau * (mmz.at(u-1) * F - E) );
    }

    mS.at(i) *= exp( -mpar.front() * myt.at(i) ) / mN;

  } // end for ( unsigned int i = 0; i < myt.size (); i++ )

} // end GUTS::calcSurvivalProbabilities()
