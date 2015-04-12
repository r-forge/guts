/**
 * Class GUTS Rcpp Implementation.
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch
 * License GPL-2
 * 2015-04-12
 */





#include "GUTS_Rcpp.h"





using namespace std;





/**
 * Constructor.
 */
GUTS::GUTS() : mInterval(4.0),
mC(), mCt(), my(), myt(),
mdist("lognormal"), mmodel("proper"),
mpar(5, 0.0), mparPos(5, 0), wpar(5, 0.0), mparMustBeRechecked(true),
mN(1000), mM(10000),
mS(0, GNAN_DOUBLE), mLL(GNAN_DOUBLE),
mz(1000, GNAN_DOUBLE), mzw(1000, GNAN_DOUBLE),
defaultGrid(1000, GNAN_DOUBLE), mExperiment(11),
D(10000, GNAN_DOUBLE),
ee(1000, GNAN_DOUBLE), ff(1000, 0),
diffy(), diffS(),
mMessages(MSG_COUNT, 0), mMessageStrings( 120, "" ) {

	/* Set messages for non-defaults. */
	int er [] = {MSG_CONC, MSG_SURV, MSG_PAR, MSG_Z, MSG_S}; // codes 11–15
	for ( unsigned int i = 0; i < (sizeof(er)/sizeof(*er)); ++i ) {
		mMessages.at(er[i]) = i + 11;
	}

	/* Set default wpar 3 to Model "it" (Inf, represented as maximum double). */
	wpar.at(2) = GMAX_DOUBLE;

	/* Default mparPos {0,1,2,3,4}. */
	for ( int i = 1; i < 5; i++ ) {
		mparPos.at(i) = i;
	}

	/* Create default grid according to N. */
	makeDefaultGrid();

	/* Message strings. */
	mMessageStrings.at(11) = "Concentrations not set.";
	mMessageStrings.at(12) = "Survivors not set.";
	mMessageStrings.at(13) = "Parameters not set.";
	mMessageStrings.at(14) = "No sample present.";
	mMessageStrings.at(15) = "No survivor probabilities present.";

	mMessageStrings.at(21) = "Vector C must be longer than 1.";
	mMessageStrings.at(22) = "Vector Ct must be longer than 1.";
	mMessageStrings.at(23) = "Vectors C and Ct must have the same length.";
	mMessageStrings.at(24) = "Vector Ct must start at 0.0.";
	mMessageStrings.at(25) = "Vector Ct must contain unique values in ascending order.";

	mMessageStrings.at(31) = "Vector y must be longer than 1.";
	mMessageStrings.at(32) = "Vector yt must be longer than 1.";
	mMessageStrings.at(33) = "Vectors y and yt must have the same length.";
	mMessageStrings.at(34) = "Vector yt must start at 0.0.";
	mMessageStrings.at(35) = "Values in vector y must not ascend.";
	mMessageStrings.at(36) = "Vector yt must contain unique values in ascending order.";

	mMessageStrings.at(41) = "Unknown distribution (valid: lognormal, delta, empirical).";
	mMessageStrings.at(46) = "Unknown model (valid: proper, it).";

	mMessageStrings.at(51) = "Wrong parameter vector length.";

	mMessageStrings.at(61) = "Parameters must be non-negative values, setting LL to -Inf.";
	mMessageStrings.at(62) = "Sample sd = 0 only if m = 0, setting LL to -Inf.";

	mMessageStrings.at(71) = "N must be greater than 2.";
	mMessageStrings.at(73) = "M must be greater than 1.";

	mMessageStrings.at(81) = "Sample must be longer than 2.";
	mMessageStrings.at(82) = "Sample must contain non-negative values.";

	mMessageStrings.at(91) = "Sampling failed.";
	mMessageStrings.at(101) = "Calculation of survival probabilities failed.";
	mMessageStrings.at(111) = "Calculation of survival probabilities failed, check yt and Ct.";

} // End of Constructor.





/**
 * Destructor.
 */
GUTS::~GUTS() {
	mC.clear();
	mCt.clear();
	my.clear();
	myt.clear();

	mpar.clear();
	mparPos.clear();
	wpar.clear();

	mS.clear();

	mz.clear();
	mzw.clear();

	defaultGrid.clear();

	D.clear();
	ee.clear();
	ff.clear();
	diffy.clear();
	diffS.clear();

	mMessages.clear();
	mMessageStrings.clear();
} // End of Destructor.





/**
 * Setters.
 */



/*
 * Setter for mC and mCt.
 */
void GUTS::setConcentrations( const vector<double>& C, const vector<double>& Ct ) {

	mMessages.at(MSG_CONC) = 0;
	if ( C.size() < 2 ) {
		makeMessage(MSG_CONC, 21); // Vector C must be longer than 1.
	} else if ( Ct.size() < 2 ) {
		makeMessage(MSG_CONC, 22); // Vector Ct must be longer than 1.
	} else if ( C.size() != Ct.size() ) {
		makeMessage(MSG_CONC, 23); // Vectors C and Ct must have the same length.
	} else if ( Ct.at(0) != 0.0 ) {
		makeMessage(MSG_CONC, 24); // Vector Ct must start at 0.0.
	}
	for ( unsigned int i = 1; i < Ct.size(); ++i ) {
		if ( (Ct.at(i) - Ct.at(i-1)) <= 0.0 ) {
			makeMessage(MSG_CONC, 25); // Vector Ct must contain unique values in ascending order.
		}
	}

	if ( mMessages.at(MSG_CONC) > 0 ) {
		mC.resize(0);
		mCt.resize(0);
	} else {
		mC = C;
		mCt = Ct;
		makeReset();
	}

} // End of GUTS::setConcentrations().



/*
 * Setter for my and myt.
 */
void GUTS::setSurvivors( const vector<int>& y, const vector<double>& yt ) {

	mMessages.at(MSG_SURV) = 0;
	if ( y.size() < 2 ) {
		makeMessage(MSG_SURV, 31); // Vector y must be longer than 1.
	} else if ( yt.size() < 2 ) {
		makeMessage(MSG_SURV, 32); // Vector yt must be longer than 1.
	} else if ( y.size() != yt.size() ) {
		makeMessage(MSG_SURV, 33); // Vectors y and yt must have the same length.
	} else if ( yt.at(0) != 0.0 ) {
		makeMessage(MSG_SURV, 34); // Vector yt must start at 0.0.
	}
	for ( unsigned int i = 1; i < y.size(); ++i ) {
		if ( (y.at(i-1) - y.at(i)) < 0 ) {
			makeMessage(MSG_SURV, 35); // Values in vector y must not ascend.
		} else if ( (yt.at(i) - yt.at(i-1)) <= 0.0 ) {
			makeMessage(MSG_SURV, 36); // Vector yt must contain unique values in ascending order.
		}
	}

	/*
	 * Assign y, yt if no errors.
	 * Assign diffy with y.diff and the last element of y.
	 * Blow up mS (survival propabilities).
	 * Blow up diffS.
	 */
	if ( mMessages.at(MSG_SURV) > 0 ) {
		my.resize( 0 );
		myt.resize( 0 );
		diffy.resize( 0 );
		diffS.resize( 0 );
	} else {
		my = y;
		diffy.assign( my.size(), my.back() );
		for ( unsigned int i = 1; i < my.size(); ++i ) {
			diffy.at(i-1) = my.at(i-1) - my.at(i);
		}
		myt = yt;
		//mS.assign( myt.size(), 0.0 ); // in makeReset().
		diffS.assign( myt.size(), 0.0 );
		makeReset();
	}

} // End of GUTS::setSurvivors().



/*
 * Setter for mdist.
 */
void GUTS::setDistribution( string dist ) {
	std::transform(dist.begin(), dist.end(), dist.begin(), ::tolower);
	if ( dist != "lognormal" && dist != "delta" && dist != "empirical" && dist != "lognormold" ) {
		makeMessage(MSG_DIST, 41); // Unknown distribution name.
		return;
	} else {
		mMessages.at(MSG_DIST) = 0;
		mdist = dist;
		makeReset();
		makeExperiment();
	}
} // End of GUTS::setDistribution().



/*
 * Setter for mmodel.
 */
void GUTS::setModel( string model ) {
	std::transform(model.begin(), model.end(), model.begin(), ::tolower);
	if ( model != "proper" && model != "it" ) {
		makeMessage(MSG_MODEL, 46); // Unknown model (valid: proper, it).
		return;
	} else {
		mMessages.at(MSG_MODEL) = 0;
		mmodel = model;
		makeReset();
		makeExperiment();
	}
} // End of GUTS::setModel().



/*
 * Make experiment (integer code) from mdist and mmodel.
 */
void GUTS::makeExperiment() {

	mparMustBeRechecked = true; // provoke re-check of parameters in calcSample()

	wpar.assign( 5, 0.0 );
	wpar.at(2) = GMAX_DOUBLE; // represent Inf with largest double
	mparPos.assign( 5, 0 ); // five zeros, fill rest with loop
	for ( int i = 1; i < 5; i++ ) {
		mparPos.at(i) = i;
	}

	/* Distribution. */
	if ( mdist == "lognormal" ) {
		mExperiment = 11;
	} else if ( mdist == "delta" ) {
		mExperiment = 21;
		mparPos.resize( 4 ); // remove 5th parameter for delta (sd)
	} else if ( mdist == "empirical" ) {
		mExperiment = 91;
		mparPos.resize( 3 ); // remove 4th and 5th
	}

	/* Model. */
	if ( mmodel == "it" ) {
		mExperiment = mExperiment + 1;
		mparPos.erase(mparPos.begin()+2); // remove third parameter
	}

	/* TODO: Testing only. */
	if ( mdist == "lognormold" ) {
		mmodel = "proper";
		mExperiment = 15;
		mparPos.assign( 5, 0 );
		for ( int i = 1; i < 5; i++ ) {
			mparPos.at(i) = i;
		}
	}

} // End of GUTS::makeExperiment().



/*
 * Sample length mN.
 */
void GUTS::setSampleLength( const int& N ) {
	if ( N < 3 ) {
		makeMessage(MSG_N, 71); // N must be greater than 2.
		return;
	} else {
		mMessages.at(MSG_N) = 0;
		mN = N;
		ee.assign( mN, 0.0 );
		ff.assign( mN, 0 );
		makeDefaultGrid();
		mz.resize( mN );
		mzw.resize( mN );
		makeReset();
	}
} // End of GUTS::setSampleLength().



/*
 * Time grid points mM.
 */
void GUTS::setTimeGridPoints( const int& M ) {
	if ( M < 2 ) {
		makeMessage(MSG_M, 73); // M must be greater than 1.
		return;
	} else {
		mMessages.at(MSG_M) = 0;
		mM = M;
		D.assign( mM, 0.0 );
		makeReset();
	}
} // End of GUTS::setTimeGridPoints().



/*
 * Direkt setter for mz.
 */
void GUTS::setSample( const vector<double>& z ) {
	if ( z.size() < 3 ) {
		makeMessage(MSG_EMP, 81); // Sample must be longer than 2.
		return;
	} else if ( *min_element( z.begin(), z.end() ) < 0.0 ) {
		makeMessage(MSG_EMP, 82); // Sample must contain non-negative values.
		return;
	} else {
		mMessages.at(MSG_EMP) = 0;
		setDistribution( "empirical" );
		setSampleLength( static_cast<int>( z.size() )); // mN, makeDefaultGrid(), z/zw resizes
		mz = z;
		sort( mz.begin(), mz.end() );
		mzw.assign( mz.size(), 0.0 );
	}
} // End of GUTS::setSample().





/**
 * Sampling.
 */



/*
 * Make sample from par and N.
 */
void GUTS::makeSample() {

	if ( mparMustBeRechecked ) setParameters( mpar ); // re-check mpars

	int er[] = {MSG_PAR, MSG_MODEL, MSG_DIST, MSG_N, MSG_EMP};
	for ( unsigned int i = 0; i < (sizeof(er)/sizeof(*er)); ++i ) {
		if ( mMessages.at(er[i]) > 0 ) {
			makeMessage(MSG_Z, 91); // "Calculation of z failed.";
			return;
		}
	}
	if ( mMessages.at(MSG_PAR_WARNING) > 0 ) return;

	/* We have not returned so far... */
	mMessages.at(MSG_Z) = 0;
	if ( mExperiment < 15 ) {
		// lognormal: 11, 12
		if (  wpar.at(4) == 0.0  ) {
			// delta-like
			mz.assign(  mN, wpar.at(3) );
			mzw.assign( mN, 0.0 );
			return;
		} else {
			makeSampleGrid();
			return;
		}
	} else if ( mExperiment == 15 ) {
		// TODO: lornormold: 15 – testing only
		if (  wpar.at(4) == 0.0  ) {
			// delta-like
			mz.assign(  mN, wpar.at(3) );
			mzw.assign( mN, 0.0 );
			return;
		} else {
			makeSampleLognormold();
			return;
		}
	} else if ( mExperiment < 30 ) {
		// delta: 21, 22
		mz.assign(  mN, wpar.at(3) );
		mzw.assign( mN, 0.0 );
		return;
	}
	
} // End of GUTS::makeSample().



/*
 * Make sample grid.
 */
void GUTS::makeSampleGrid() {
	double sigma2   =  log(   1.0  +  pow( (wpar.at(4) / wpar.at(3)), 2.0 )   );
	double mu       =  log(wpar.at(3))  -  (0.5 * sigma2);
	double sigmaD   =  sqrt(sigma2) * mInterval;
	for ( unsigned int i = 0; i < mN; ++i ) {
		mz.at(i)  =  exp( defaultGrid.at(i) * sigmaD + mu);
		mzw.at(i) =  -0.5  *  pow((defaultGrid.at(i) * sigmaD), 2.0)  /  sigma2;
	}
} // End of GUTS::makeSampleGrid().



/*
 * Make lognormold sample.
 */
void GUTS::makeSampleLognormold() {
	using namespace boost;
	static boost::random::mt19937 rng( static_cast<unsigned> (std::time(0)) );
	boost::lognormal_distribution<double> lognorm_dist( wpar.at(3), wpar.at(4) );
	variate_generator< boost::random::mt19937&, boost::lognormal_distribution<double> > lognormal_sampler(rng, lognorm_dist);
	for ( int i = 0; i < mN; ++i ) {
		mz.at(i) = lognormal_sampler();
	}
	sort( mz.begin(), mz.end() );
	mzw.assign( mN, 0.0 );
} // End of GUTS::makeSampleLognormold().





/**
 * Survival Probabilities.
 */



/*
 * Survival probabilities standard implementation.
 */
void GUTS::calcSurvivalProbabilities() {

	makeSample();

	int er[] = {MSG_Z, MSG_SURV, MSG_CONC, MSG_M};
	for ( unsigned int i = 0; i < (sizeof(er)/sizeof(*er)); ++i ) {
		if ( mMessages.at(er[i]) > 0 ) {
			makeMessage(MSG_S, 101); // Calculation of survival probabilities failed.
			return;
		}
	}
	if ( mMessages.at(MSG_PAR_WARNING) > 0 ) return;

	/* We have not returned so far... */
	mMessages.at(MSG_S) = 0;
	if ( mCt.back() < myt.back() ) {
		makeMessage(MSG_S, 111); // Calculation of survival probabilities failed, check yt and Ct.
		return;
	}

	/* We have not returned so far... */
	makeSurvivalProbabilities();

} // End of GUTS::calcSurvivalProbabilities().



/*
 * Survival probabilities implementation with vector of survival time points.
 */
void GUTS::calcSurvivalProbabilities( const vector<double>& St ) {
	vector<int> hy;
	hy.assign(St.size(), 0); // 0-vector of survivors
	setSurvivors(hy, St); // standard setter
	calcSurvivalProbabilities(); // standard method
} // End of GUTS::calcSurvivalProbabilities().



/*
 * Implementation with number of survival time points.
 */
void GUTS::calcSurvivalProbabilities( const unsigned int& Stlength ) {
	vector<int> hy;
	hy.assign( Stlength, 0 ); // 0-vector of survivors
	vector<double> hyt;
	for ( unsigned int i=0; i<Stlength; ++i ) {
		hyt.push_back(i); // ascending surivor time points
	}
	setSurvivors(hy, hyt); // standard setter
	calcSurvivalProbabilities(); // standard method
} // End of GUTS::calcSurvivalProbabilities().



/*
 * Maker of Survival Probabilities.
 */
void GUTS::makeSurvivalProbabilities() {

	/* Method-relevant vectors and variables. */
	D.assign( mM, 0.0 );
	ee.assign( mN, 0.0 );
	ff.assign( mN, 0 );
	double dtau = (myt.back() - myt.front()) / mM;
	double Scale = 0.0; // to prevent warnings
	double tmp;
	double carlosKlammer;
	double E;
	double F;

	/* Variables, iterators. */
	int dpos    = 0;
	int ii      = 0;
	int zpos    = 0;
	double tau  = 0.0;
	int k       = 0;

	/* Loop over yt_i's. */
	for ( unsigned int i = 0; i < myt.size(); ++i ) {

		while ( tau < myt.at(i) && dpos < mM ) {

			tmp = exp( -mpar.at(1) * (tau - mCt.at(k)) );
			if ( mpar.at(1) > 0.0 ) {
				carlosKlammer = (1.0 / mpar.at(1)) * (1.0-tmp);
			} else {
				carlosKlammer = tau - mCt.at(k);
			}

			D.at(dpos) = D.at(ii) * tmp
			           + mC.at(k) * (1.0-tmp)
			           + (mC.at(k+1) - mC.at(k))  *  (tau - mCt.at(k) - carlosKlammer)  /  (mCt.at(k+1) - mCt.at(k));

			/* Increment or decrement zpos. */
			if ( D.at(dpos) < mz.at(zpos) ) {
				while ( D.at(dpos) < mz.at(zpos) && zpos > 0 ) {
					--zpos;
				}
				if ( D.at(dpos) > mz.at(zpos) ) {
					++zpos;
				}
			} else {
				while ( D.at(dpos) > mz.at(zpos) && zpos < (mN-1) ) {
					++zpos;
				}
			}

			/* Update ee and ff. */
			if ( D.at(dpos) > mz.at(mN-1) ) {
				ee.at(mN-1) += D.at(dpos);
				ff.at(mN-1)++;
			} else if ( D.at(dpos) > mz.at(0) ) {
				ee.at(zpos-1) += D.at(dpos);
				ff.at(zpos-1)++;
			}

			/* Increment or decrement dpos, tau, k. */
			dpos++;
			tau += dtau;
			if ( tau > mCt.at(k+1) ) {
				k++;
				ii = dpos-1;
			}

		} // End of while ( tau < myt.at(i) )

		/* Write mS at i. */
		E = ee.back();
		F = ff.back();
		mS.at(i) = exp(   (mpar.at(2) * dtau * (mz.back() * F - E))  +  mzw.back()   );
		for ( unsigned int u=mN-1; u > 0; --u ) {
			E += ee.at(u-1);
			F += ff.at(u-1);
			mS.at(i) += exp(   (mpar.at(2) * dtau * (mz.at(u-1) * F - E))  +  mzw.at(u-1)   );
		}
		mS.at(i) *= exp( -mpar.front() * myt.at(i) ) / mN;

		/* Scale S, build diffS. */
		if ( i == 0 ) {
			Scale = 1.0 / mS.at(0);
		} else {
			diffS.at(i-1) = mS.at(i-1) - Scale * mS.at(i);
		}
		mS.at(i) *= Scale;
		
	} // End of for ( unsigned int i = 0; i < myt.size (); i++ )

	/* Set the last diffS = mS.back() - 0 = mS.back(). */
	diffS.back() = mS.back();

} // End of GUTS::makeSurvivalProbabilities().





/**
 * Loglikelihood.
 */
void GUTS::calcLoglikelihood() {

	calcSurvivalProbabilities();

	if ( mMessages.at(MSG_S) > 0 ) {
		mLL = GNAN_DOUBLE; // error
		return;
	}
	if ( mMessages.at(MSG_PAR_WARNING) > 0 ) {
		mLL = -GINF_DOUBLE; // warning
		return;
	}

	/* We have not returned so far... */
	mLL = 0.0;
	for ( unsigned int i=0; i < diffy.size(); ++i ) {
		if ( diffy.at(i) > 0 ) {
			if ( diffS.at(i) == 0.0 ) {
				mLL = -GINF_DOUBLE;
				return;
			}
			mLL += diffy.at(i) * log(diffS.at(i));
		}
	} // End of for ( unsigned int i=0; i < diffy.size(); ++i ).
	
} // End of GUTS::calcLoglikelihood().









/**
 * Rcpp Module modguts.
 */
using namespace Rcpp;

RCPP_MODULE(modguts) {

	class_<GUTS>( "GUTS" )

		.constructor()

		.method( "setConcentrations",   &GUTS::setConcentrations, "Set vector of concentrations and time points." )
		.method( "setSurvivors",        &GUTS::setSurvivors, "Set vector of survivors and time points." )

		.method( "setDistribution",     &GUTS::setDistribution, "Set sample distribution." )
		.method( "setModel",            &GUTS::setModel, "Set model." )
		.method( "setParameters",       &GUTS::setParameters, "Set parameter vector." )

		.method( "setTimeGridPoints",   &GUTS::setTimeGridPoints, "Set number of grid points on the time-axis." )
		.method( "setSampleLength",     &GUTS::setSampleLength, "Set length of sample (numerical exactness)." )

		.method( "setSample",           &GUTS::setSample, "Set sample vector manually." )

		.method( "calcSurvivalProbabilities", &GUTS::calcSurvivalProbabilities, "Calculate the survivor probabilities." )
		.method( "calcLoglikelihood",         &GUTS::calcLoglikelihood, "Calculate the logarithm of the likelihood." )

		.method( "setInterval",         &GUTS::setInterval, "Set interval for uniform dist (testing only)." )
		.property( "Interval",     &GUTS::getInterval, "Interval." )

		.property( "C",     &GUTS::getC, "Vector of concentrations." )
		.property( "Ct",    &GUTS::getCt, "Vector of concentration times." )
		.property( "y",     &GUTS::gety, "Vector of survivors." )
		.property( "yt",    &GUTS::getyt, "Vector of survivors times." )
		.property( "dist",  &GUTS::getdist, "Distribution." )
		.property( "model", &GUTS::getmodel, "Model." )
		.property( "par",   &GUTS::getpar, "Parameter vector." )
		.property( "M",     &GUTS::getM, "Grid points on time axis." )
		.property( "N",     &GUTS::getN, "Sample length (numerical exactness)." )

		.property( "S",     &GUTS::getS, "Survival probabilities." )
		.property( "LL",    &GUTS::getLL, "Loglikelihood." )

		.property( "Messages", &GUTS::getMessages, "Messages from GUTS object." )
		;
}





/*
 * End of GUTS_Rcpp.cpp
 */
