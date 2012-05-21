/** Class GUTS.
 * @file        GUTS.h
 * @author      soeren.vogel@uzh.ch
 * @author      carlo.albert@eawag.ch, 
 * @date        2012-05-11
 * @license     GPL-3
 *
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model
 */

#ifndef guts_h
#define guts_h guts_h

#include <vector>
#include <valarray>
#include <numeric>
#include <algorithm>
#include <string>
#include <map>
#include <ctime>  // For random number generator.
#include <boost/random.hpp>
#include <boost/math/distributions.hpp>

using namespace std;

// This does not work, instead sets GNAN_INT to 0
//const int GNAN_INT = std::numeric_limits<int>::quiet_NaN();
const int         GNAN_INT    = 0;
const double      GNAN_DOUBLE = numeric_limits<double>::quiet_NaN();
const std::string GNAN_STRING = "NaN";


/*
 * The Class GUTS
 */
class GUTS {

public:

  GUTS();
  ~GUTS();

  /** Title of the Experiment
   */
  void setTitle(const string Title) {
    mTitle = Title;
  }

  /** Set time series vector of concentrations.
   * @param   C   Vector of concentrations.
   * @param   Ct  Vector of time points of concentrations (must start at 0).
   * Must be of equal length. Set mC and mCt.
   * Error: GER_C, GER_CT
   */
  void setConcentrations(const vector<double>& C, const vector<double>& Ct);

  /** Set time series vector of survivors.
   * @param   y   Vector of survivors.
   * @param   yt  Vector of time points of survivors (must start at 0).
   * Must be of equal length. Set my and myt.
   * Error: GER_Y, GER_YT
   */
  void setSurvivors(const vector<int>& y, const vector<double>& yt);

  /** Set parameter vector of the object.
   * @param   par     Vector of numerics.
   * Vector must be of size > 2. First three are required for the calculation
   * of the loglikelihood, the latter are used for the construction of sample
   * vector z. Set mpar.
   * Error: GER_PAR
   */
  void setParameters(const vector<double>& par);

  /** Set number of grid points on the time-axis.
   * @param   M   Integer.
   * Grid points of the time axis used for the integration. Set mM.
   * Error: GER_M
   */
  void setTimeGridPoints(const int& M);

  /** The (name of the) distribution to sample from.
   * @param:  dist    string
   * Set mdist.
   * Error: GER_DIST
   */
  void setDistribution(const string dist);

  /** Size of the sample.
   * @param   N   Integer.
   * Set mN.
   * Error: GER_N
   */
  void setSampleLength(const int& N);

  /** Sample vector.
   * @param   z    Vector of numerics.
   * Set mz.
   * Error: GER_Z
   */
  void setSample(const vector<double> &z);

  /** Construct the sample, sort, and save result in mz.
   * Set mz.
   * Error: GER_Z
   */
  void calcSample();

  /** Calculate survival probabilities and save result in mS.
   * Overloaded. Use the method with either of the following arguments:
   * @param   <none>    Standard invokation. Use local fields for calculation.
   * @param   St        Survival time points, defaults to myt.
   * @param   Stlength  Number of survival time point, defaults to myt.size().
   * Set mS.
   * Error: GER_S
   */
  void calcSurvivalProbabilities();
  void calcSurvivalProbabilities(const unsigned int& Stlength);
  void calcSurvivalProbabilities(const vector<double>& St);

  /** Calculate loglikelihood.
   * Set mLL.
   * Error: GER_LL
   */
  void calcLoglikelihood();

  /** Getters.
   * Simply return contents of field.
   */
  string          getTitle()  const { return mTitle; }
  vector<double>  getC()      const { return mC; }
  vector<double>  getCt()     const { return mCt; }
  vector<int>     gety()      const { return my; }
  vector<double>  getyt()     const { return myt; }
  vector<double>  getpar()    const { return mpar; }
  int             getM()      const { return mM; }
  string          getdist()   const { return mdist; }
  int             getN()      const { return mN; }
  vector<double>  getz()      const { return mz; }
  vector<double>  getD()      const { return mD; }
  vector<double>  getS()      const { return mS; }
  double          getLL()     const { return mLL; }
  vector<bool>    getErrors() const { return mErrors; }
  vector<string>  getErrorMessages()
    const { return mErrorMessages; }
  
  /** showObject()
   * Print complete object in a formatted way to cout.
   * This method is only used in the console version of GUTS.
   * The R version has its own print method.
   */
  void showObject();

protected:
  
  /** Attributes of a GUTS object.
   * See method description for details on the variables.
   */
  string mTitle;                  // Title of the experiment.
  vector<double> mC;              // Concentrations
  vector<double> mCt;             // Time points of concentratinos
  vector<int> my;                 // Survivors
  vector<double> myt;             // Survivor time points
  vector<double> mpar;            // Parameters
  int mM;                         // Time grid points
  string mdist;                   // Name of distribution
  int mN;                         // Length of z (sample)
  vector<double> mz;              // Actual sample
  vector<double> mD;              // Damage
  vector<double> mS;              // Survival probabilities
  double mLL;                     // The loglikelihood.

  /** do-Methods
   * Work horses for specific samples and the
   * calculation of the survival probabilities.
   */
  void doCalcSampleLognormal
    (const double& mean, const double& sigma, const int& n);
  void doCalcSurvivalProbabilities();

  /** Errors.
   * mErrors are true, if errors exist, false, if no errors exist!
   */
  vector<bool> mErrors;
  vector<string> mErrorMessages;
  enum
  {
    GER_C = 0,
    GER_CT,
    GER_Y,
    GER_YT,
    GER_PAR,
    GER_M,
    GER_DIST,
    GER_N,
    GER_Z,
    GER_D,
    GER_S,
    GER_LL,
    
    GER_COUNT = 12 // ***RE-COUNT if you add or remove error positions!!!***
  };

}; // end class GUTS

#endif
