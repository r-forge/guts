/**
 * Class GUTS Rcpp.
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch
 * License GPL-2
 * 2015-04-12
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 */





#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include <Rcpp.h>





using namespace std;





/**
 * Constants.
 */
const double  GNAN_DOUBLE = std::numeric_limits<double>::quiet_NaN();
const double  GMAX_DOUBLE = std::numeric_limits<double>::max();
const double  GINF_DOUBLE = std::numeric_limits<double>::infinity();





/**
 * The Class GUTS.
 */
class GUTS {

public:



	/**
	 * Standard Constructor and Destructor.
	 */
	GUTS();
	~GUTS();



	/**
	 * Setters.
	 */
	void setConcentrations( const vector<double>& C, const vector<double>& Ct );
	void setSurvivors( const vector<int>& y, const vector<double>& yt );

	void setDistribution( string dist );
	void setModel( string model );

	void setParameters( const vector<double>& par ) {
		mpar = par; // take whatever the user provides
		mparMustBeRechecked = false; // because we just do this
		if ( par.size() != mparPos.size() ) {
			makeMessage(MSG_PAR, 51); // Wrong vector length.
			return;
		} else if ( *min_element( par.begin(), par.begin()+mparPos.size() ) < 0.0 ) {
			makeMessage(MSG_PAR_WARNING, 61); // Parameters must be non-negative values.
			return;
		} else {
			for ( unsigned int i = 0; i < mparPos.size(); ++i ) {
				wpar.at(mparPos.at(i)) = par.at(i); // assign to wpar at positions in mparPos
			}
			if ( wpar.at(3) == 0.0 && wpar.at(4) != 0.0 ) {
				makeMessage(MSG_PAR_WARNING, 62); // Sample sd = 0 only if m = 0.
				return;
			} else if ( std::isinf(wpar.at(2)) ) {
				wpar.at(2) = GMAX_DOUBLE;
			}
			mMessages.at(MSG_PAR_WARNING) = 0;
			mMessages.at(MSG_PAR) = 0;
		}
	}

	void setSampleLength( const int& N );
	void setTimeGridPoints( const int& M );

	void setSample( const vector<double>& z );



	/**
	 * Calculators.
	 */
	void calcSurvivalProbabilities();
	void calcSurvivalProbabilities( const unsigned int& Stlength );
	void calcSurvivalProbabilities( const vector<double>& St );
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

	/* Get message strings where set. */
	vector<string> getMessages() {
		vector<string> ret;
		for ( unsigned int i = 0; i < mMessages.size(); ++i ) {
			if ( mMessages.at(i) > 0 ) {
				ret.push_back( mMessageStrings.at(mMessages.at(i)) );
			}
		}
		return ret;
	}



	/**
	 * C++-only Methods.
	 */
	//void showObject();



	/**
	 * TODO Devel and Testing.
	 */
	void setInterval( const double& Interval ) {
		if ( Interval > 1.0 ) {
			mInterval = Interval;
		}
	}
	double         getInterval() const { return mInterval; }
	double         getmExperiment() const { return mExperiment; }
	vector<double> getz() const { return mz; }

protected:



	/**
	 * TODO Devel and Testing.
	 */
	double mInterval;
	void makeSampleLognormold();



	/**
	 * Fields.
	 */
	vector<double> mC;
	vector<double> mCt;
	vector<int>    my;
	vector<double> myt;

	string mdist;
	string mmodel;
	vector<double> mpar;
	vector<int>    mparPos;
	vector<double> wpar;
	bool           mparMustBeRechecked;

	int mN;
	int mM;

	vector<double> mS;
	double mLL;

	vector<double> mz;
	vector<double> mzw;

	vector<double> defaultGrid;

	int mExperiment;

	vector<double> D;
	vector<double> ee;
	vector<int>    ff;

	vector<int>    diffy;
	vector<double> diffS;



	/**
	 * Maker Methods.
	 */
	void makeSample();
	void makeDefaultGrid() {
		defaultGrid.resize( mN );
		for ( int i = 0; i < mN; ++i ) {
			defaultGrid.at(i) = (2.0 * i - mN + 1.0) / (mN-1.0);
		}
	}
	void makeExperiment();
	void makeSampleGrid();
	void makeSurvivalProbabilities();
	void makeMessage( const unsigned int& messPos, const unsigned int& messCode ) {
		mMessages.at(messPos) = messCode;
		makeReset();
		if ( messPos == MSG_PAR_WARNING ) {
			Rcpp::warning(mMessageStrings.at(messCode));
		} else {
			Rcpp::stop(mMessageStrings.at(messCode));
		}
	}
	void makeReset() {
		mS.assign( myt.size(), GNAN_DOUBLE );
		mLL = GNAN_DOUBLE;
	}



	/**
	 * Messages.
	 */
	vector<unsigned int> mMessages;
	vector<string>       mMessageStrings;
	enum {
		MSG_CONC = 0,
		MSG_SURV,
		MSG_DIST,
		MSG_MODEL,
		MSG_PAR,
		MSG_PAR_WARNING,
		MSG_N,
		MSG_M,
		MSG_EMP,
		MSG_Z,
		MSG_S,
		MSG_COUNT = 11 // recount if change, for filling mMessages.
	};

}; // End of class GUTS





/**
 * End of GUTS.h.
 */
