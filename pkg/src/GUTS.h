/**
 * Class GUTS.
 *
 * @file    GUTS.h
 * @author  soeren.vogel@posteo.ch
 * @author  carlo.albert@eawag.ch
 * @date    2015-04-02
 * @license GPL-2
 *
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
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





/**
 * Constants.
 * NaN, Max, Inf
 * Only used expressions.
 */

/* NaN. */
// NaN for int is always 0. See: http://stackoverflow.com/questions/3949457/can-an-integer-be-nan-in-c
//const int         GNAN_INT    = 0;
const double      GNAN_DOUBLE = std::numeric_limits<double>::quiet_NaN();
const std::string GNAN_STRING = "NaN";

/* Max. */
const double      GMAX_DOUBLE = std::numeric_limits<double>::max();

/* Inf. */
const double      GINF_DOUBLE = std::numeric_limits<double>::infinity();

/*
 * More Constants, see
 * * http://www.cplusplus.com/reference/cfloat/
 * * http://en.cppreference.com/w/cpp/types/numeric_limits
 * DBL_MAX_EXP = largest double.
 * DBL_EPSILON
 * DBL_MIN
 */





/**
 * The Class GUTS.
 */
class GUTS {

public:

	GUTS();
	~GUTS();

	/**
	 * Set Time Series Vector of Concentrations.
	 * @param   C   Vector of concentrations.
	 * @param   Ct  Vector of time points of concentrations (must start at 0).
	 * Must be of equal length. Set mC and mCt.
	 * Error: GER_C, GER_CT
	 */
	void setConcentrations( const vector<double>& C, const vector<double>& Ct );

	/**
	 * Set Time Series Vector of Survivors.
	 * @param   y   Vector of survivors.
	 * @param   yt  Vector of time points of survivors (must start at 0).
	 * Must be of equal length. Set my and myt.
	 * Error: GER_Y, GER_YT
	 */
	void setSurvivors( const vector<int>& y, const vector<double>& yt );

	/**
	 * The (Name of the) Distribution to Sample From.
	 * @param:  dist    string
	 * Distributions:
	 * * lognormal: Sample from Grid.
	 * * delta: SD = 0.
	 * * empirical: Give vector z.
	 * * logold: For Testing only.
	 * Set mdist.
	 * Error: GER_DIST
	 */
	void setDistribution( const string dist );

	/**
	 * Set GUTS Model: proper, it.
	 * @param   model   string
	 */
	void setModel( const string model );

	/**
	 * Set Parameter Vector of the Object.
	 * @param   par     Vector of numerics.
	 * Vector must be of size > 2. First three are required for the calculation
	 * of the loglikelihood, the latter are used for the construction of sample
	 * vector z. Set mpar.
	 * Parameters:
	 * * h_b: background mortality rate
	 * * k_r: recovery rate
	 * * k_k: killing rate
	 * * m: sample mean
	 * * sd: sample standard deviation
	 * Error: GER_PAR
	 */
	void setParameters( const vector<double>& par );

	/**
	 * Size of the Sample.
	 * @param   N   Integer.
	 * Set mN.
	 * Error: GER_N
	 */
	void setSampleLength( const int& N );

	/**
	 * Set number of grid points on the time-axis.
	 * @param   M   Integer.
	 * Grid points of the time axis used for the integration. Set mM.
	 * Error: GER_M
	 */
	void setTimeGridPoints( const int& M );

	/**
	 * Sample Vector.
	 * Set sample manually and sort.
	 * @param   z    Vector of numerics.
	 * Set mz.
	 * Set weights mzw, length of mz, all 0.
	 * Error: GER_Z
	 */
	void setSample( const vector<double>& z );

	/**
	 * Calculate Survival Probabilities and Save Result in MS.
	 * Overloaded. Use the method with either of the following arguments:
	 * @param   <none>    Standard invokation. Use local fields for calculation.
	 * @param   St        Survival time points, defaults to myt.
	 * @param   Stlength  Number of survival time point, defaults to myt.size().
	 * Set mS.
	 * Error: GER_S
	 */
	void calcSurvivalProbabilities();
	void calcSurvivalProbabilities( const unsigned int& Stlength );
	void calcSurvivalProbabilities( const vector<double>& St );

	/**
	 * Calculate loglikelihood.
	 * Set mLL.
	 * Error: GER_LL
	 */
	void calcLoglikelihood();

	/**
	 * Getters.
	 */
	vector<double>  getC()             const { return mC; }
	vector<double>  getCt()            const { return mCt; }
	vector<int>     gety()             const { return my; }
	vector<double>  getyt()            const { return myt; }

	string          getdist()          const { return mdist; }
	string          getmodel()         const { return mmodel; }
	vector<double>  getpar()           const { return mpar; }

	int             getN()             const { return mN; }
	int             getM()             const { return mM; }

	vector<double>  getS()             const { return mS; }
	double          getLL()            const { return mLL; }

	vector<bool>    getMessages()      const { return mMessages; }

	/**
	 * Method showObject().
	 * Print complete object in a formatted way to cout.
	 * This method is only used in the console version of GUTS.
	 * The R version has its own print method.
	 */
	void showObject();

	/**
	 * TODO FOR TESTING ONLY.
	 */
	void           setInterval( const double& Interval );
	double         getInterval() const { return mInterval; }
	double         getmExperiment() const { return mExperiment; }
	vector<double> getz() const { return mz; }

protected:

	/**
	 * TODO FOR TESTING ONLY.
	 */
	double mInterval;
	void doSampleLognormold();

	/**
	 * Protected Fields.
	 */
	
	/*
	 * Concentrations, Survivors and Repective Time Points.
	 * Default/fallback resize(0).
	 */
	vector<double> mC;
	vector<double> mCt;
	vector<int>    my;
	vector<double> myt;

	/*
	 * Distribution Name.
	 * Values: "lognormal", "delta", "empirical".
	 * Default/fallback: "lognormal".
	 *
	 * Model.
	 * Values: "proper", "it".
	 * If dist=="lornormold" (testing), model==proper.
	 * Default/fallback: "proper".
	 */
	string mdist;
	string mmodel;

	/*
	 * Parameters.
	 * * h_b: background mortality rate
	 * * k_r: recovery rate
	 * * k_k: killing rate
	 * * m: sample mean
	 * * sd: sample standard deviation
	 * Default: assign(5, NaN).
	 * Fallback: last parameters provided.
	 *
	 * Positions of necessary parameters.
	 * Default/fallback: { 0, 1, 2, 3, 4 }.
	 *
	 * Working representation of parameters.
	 * Alway length 5.
	 * Default/fallback: { 0.0, 0.0, DBL_MAX_EXP, 0.0, 0.0 }
	 * Represent Inf on k_k with largest double.
	 */
	vector<double> mpar;
	vector<int>    mparPos;
	vector<double> wpar;

	/*
	 * Length of Sample Vector (z).
	 * Default: 1’000.
	 * Fallback: Last valid value.
	 *
	 * Time grid points.
	 * Default: 10’000.
	 * Fallback: Last valid value.
	 */
	int mN;
	int mM;

	/*
	 * Sample.
	 * Default: assign( mN, 0.0 ), where mN is 1000.
	 * Fallback: assign( mN, 0.0 ).
	 *
	 * Sample Weights.
	 * Default: assign( mN, 0.0 ), where mN is 1000.
	 * Fallback: assign( mN, 0.0 ).
	 *
	 * Resize on setSampleLength.
	 */
	vector<double> mz;
	vector<double> mzw;

	/*
	 * Survival Probabilities.
	 * Set in setSurvivors.
	 * Default/fallback: assign( yt.size(), 0.0 ).
	 */
	vector<double> mS;

	/*
	 * The loglikelihood.
	 * Default/fallback: GNAN_DOUBLE.
	 */
	double mLL;

	/*
	 * Default Grid.
	 * Constructed when Object constructed.
	 * Changed when sample length changes.
	 * Method makeDefaultGrid.
	 * Default/fallback: makeDefaultGrid depending on mN.
	 */
	vector<double> defaultGrid;

	/*
	 * Experiment.
	 * Experiment model as integer:
	 *               Proper   IT
	 * lognormal     11       12
	 * lognormold    15       15 // intentionally the same model
	 * delta         21       22
	 * empirical     91       92
	 * Default 11.
	 * Fallback 11.
	 */
	int mExperiment;

	/*
	 * Internal vectors in doCalcSurvivalProbabilities()
	 */
	vector<double> D;   // Vector of damages. Set in setTimeGridPoints. Default assign(mM, 0.0).
	vector<double> ee;  // Set in setSampleLength and setSample. Default assign(mN, 0.0).   TODO: What is this good for.
	vector<int>    ff;  // Set in setSampleLength and setSample. Default assign(mN, 0).     TODO: What is this good for.

	/*
	 * Internal vectors for calcLoglikelihood().
	 */
	vector<int>    diffy;  // Differences between survivors, set with y, defaults forward_diff(y), append 0.
	vector<double> diffS;  // Differences between S, set with y, defaults assign(yt.size(), 0.0), append 0.

	/**
	 * Make Sample Vector.
	 *
	 * Create Sample Vector from set parameters.
	 * Sort and save result in mz.
	 * Set weights mzw.
	 * Error: GER_Z
	 */
	void calcSample();

	/**
	 * Make default Grid.
	 * Set mExperiment and mnpar.
	 */
	void makeDefaultGrid();
	void makeExperiment();

	/**
	 * do-Methods.
	 * Work horses for specific samples and the
	 * calculation of the survival probabilities.
	 */
	void doSampleGrid();
	void doCalcSurvivalProbabilities();

	/**
	 * Messages.
	 * mErrors are true, if errors/warnings exist, false, if no errors/warnings exist!
	 * GWA_PAR catches inappropriate pars, this is not an error but LL is -Inf.
	 */
	vector<bool> mMessages;
	enum {
		GER_C = 0,
		GER_CT,
		GER_Y,
		GER_YT,
		GER_DIST,
		GER_MODEL,
		GER_PAR, // Error on par.
		GWA_PAR, // Waring on par.
		GER_N,
		GER_M,
		GER_EMP,
		GER_Z,
		GER_S,
		GER_COUNT = 13 // ***RE-COUNT if you add or remove error positions!!!***
	};

}; // End of class GUTS

#endif


