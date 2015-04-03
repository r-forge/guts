/**
 * GUTS Class Implementation.
 *
 * @file    GUTS.cpp
 * @author  soeren.vogel@posteo.ch
 * @author  carlo.albert@eawag.ch
 * @date    2015-04-02
 * @license	GPL-2
 */





#ifndef guts_h
#include "GUTS.h"
#endif





using namespace std;





/**
 * Constructor.
 */
GUTS::GUTS() : mInterval(4.0),
		mC(), mCt(), my(), myt(),
		mdist("lognormal"), mmodel("proper"),
		mpar(5, GNAN_DOUBLE), mparPos(5, 0), wpar(5, 0.0),
		mN(1000), mM(10000),
		mz(1000, 0.0), mzw(1000, 0.0),
		mS(), mLL(GNAN_DOUBLE),
		defaultGrid(1000, 0.0), mExperiment(11),
		D(10000, 0.0),
		ee(1000, 0.0), ff(1000, 0),
		diffy(), diffS(),
		mMessages(GER_COUNT, true) {

	/*
	 * Reset errors for fields with valid defaults.
	 */
	int er[] = { GER_DIST, GER_MODEL, GWA_PAR, GER_N, GER_M, GER_EMP };
	for ( unsigned int i = 0; i < 6; ++i ) {
		mMessages.at( er[i] ) = false;
	}

	/*
	 * Set default wpar 3 to Model "it" (Inf, represented as maximum double).
	 */
	wpar.at(2) = GMAX_DOUBLE; // largest double

	/*
	 * Default mparPos {0,1,2,3,4}.
	 */
	for ( int i = 1; i < 5; i++ ) {
		mparPos.at(i) = i;
	}

	/*
	 * Create default grid according to N.
	 */
	makeDefaultGrid();

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

	mz.clear();
	mzw.clear();

	mS.clear();

	defaultGrid.clear();

	D.clear();
	ee.clear();
	ff.clear();
	diffy.clear();
	diffS.clear();

	mMessages.clear();
} // End of Destructor.




/**
 * Make Default Grid.
 * Internal Method.
 * Used on Initialize(), and on setSampleLength().
 */
void GUTS::makeDefaultGrid() {
	defaultGrid.resize( mN );
	for ( int i = 0; i < mN; ++i ) {
		defaultGrid.at(i) = (2.0 * i - mN + 1.0) / (mN-1.0);
	}
}





/**
 * mC and mCt
 * Message read: none
 * Messages write: GER_C, GER_CT
 */
void GUTS::setConcentrations( const vector<double>& C, const vector<double>& Ct ) {

	/*
	 * Reset messages on this method.
	 */
	int er[] = {GER_C, GER_CT};
	for ( unsigned int i = 0; i < 2; ++i ) {
		mMessages.at(er[i])   = false; // false means okay
	}

	/*
	 * New messages on C and Ct.
	 */
	if ( C.size() < 2 ) {
		mMessages.at(GER_C)   = true; // "Vector C must be longer than 1.";
	}
	if ( Ct.size() < 2 ) {
		mMessages.at(GER_CT)  = true; // "Vector Ct must be longer than 1.";
	}
	if ( C.size() != Ct.size() ) {
		mMessages.at(GER_C)   = true; // "Vectors C and Ct must have the same length.";
		mMessages.at(GER_CT)  = true; // "Vectors C and Ct must have the same length.";
	}
	if ( Ct.at(0) != 0.0 ) {
		mMessages.at(GER_CT)  = true; // "Vector Ct must start at 0.0.";
	}
	// Prevent division by 0.
	for ( unsigned int i = 1; i < Ct.size(); ++i ) {
		if ( (Ct.at(i) - Ct.at(i-1)) <= 0.0 ) {
			mMessages.at(GER_CT)  = true; // "Vector Ct must contain unique values in ascending order.";
		}
	}

	/*
	 * Assign C and Ct, if no errors, false means no errors!
	 */
	if ( mMessages.at(GER_C) || mMessages.at(GER_CT) ) {
		mC.resize(0);
		mCt.resize(0);
	} else {
		mC = C;
		mCt = Ct;
	}

} // End of GUTS::setConcentrations()





/**
 * my and myt.
 * Messages read: none
 * Messages write: GER_Y, GER_YT
 */
void GUTS::setSurvivors( const vector<int>& y, const vector<double>& yt ) {

	/*
	 * Reset messages on this method.
	 */
	int er[] = {GER_Y, GER_YT};
	for ( unsigned int i = 0; i < 2; ++i ) {
		mMessages.at( er[i] )  =  false; // false means okay
	}

	/*
	 * New messages on y and yt.
	 */
	if ( y.size() < 2 ) {
		mMessages.at( GER_Y )  = true; // "Vector y must be longer than 1.";
	}
	if ( yt.size() < 2 ) {
		mMessages.at( GER_YT )  = true; // "Vector yt must be longer than 1.";
	}
	if ( y.size() != yt.size() ) {
		mMessages.at( GER_Y )   = true; // "Vectors y and yt must have the same length.";
		mMessages.at( GER_YT )  = true; // "Vectors yt and y must have the same length.";
	}
	if ( yt.at(0) != 0.0 ) {
		mMessages.at( GER_YT )  = true; // "Vector yt must start at 0.0.";
	}
	// Prevent division by 0
	for ( unsigned int i = 1; i < yt.size(); ++i ) {
		if ( (yt.at(i) - yt.at(i-1)) <= 0.0 ) {
			mMessages.at( GER_YT )   = true; // "Vector yt must contain unique values in ascending order.";
		}
	}

	/*
	 * Assign y if okay.
	 * Assign diffy with y.diff and the last element of y.
	 */
	if ( ! mMessages.at( GER_Y ) ) {
		my = y;
		diffy.assign( my.size(), my.back() );
		for ( unsigned int i = 1; i < my.size(); ++i ) {
			diffy.at(i-1) = my.at(i-1) - my.at(i);
		}
	} else {
		my.resize( 0 );
		diffy.resize( 0 );
	}

	/*
	 * Assign yt if okay.
	 * Blow up mS (survival propabilities).
	 * Blow up diffS.
	 */
	if ( ! mMessages.at( GER_YT ) ) {
		myt = yt;
		mS.assign( myt.size(), 0.0 );
		diffS.assign( myt.size(), 0.0 );
	} else {
		myt.resize( 0 );
		mS.resize( 0 );
		diffS.resize( 0 );
	}

} // End of GUTS::setSurvivors()





/*
 * Set string dist.
 * Arguments: lognormal, delta, empirical, lognormold.
 * Lognormold is the old lornormal with sampler.
 * Messages write: GER_DIST.
 */
void GUTS::setDistribution( const string dist ) {

	/*
	 * Local copy (to edit), lowercase.
	 */
	string hdist = dist;
	std::transform(hdist.begin(), hdist.end(), hdist.begin(), ::tolower);

	/*
	 * If dist changes.
	 * If dist doesn’t change, these check are already passed.
	 * Default dist is "lognormal".
	 */
	if ( hdist != mdist ) {

		/*
		 * New messages on dist.
		 * Else assign and make experiment.
		 */
		if ( hdist != "lognormal" && hdist != "delta" && hdist != "empirical" && hdist != "lognormold" ) {
			mMessages.at( GER_DIST )  = true; // "Unknown distribution name.";
			return;
		} else {
			mMessages.at( GER_DIST )  = false;
			mdist = hdist;
			makeExperiment();
		}

	} // End of if ( hdist != mdist ).

} // End of GUTS::setDistribution( dist ).





/**
 * Method setModel
 * @param   model (default: "proper", other: it).
 * Set mmodel.
 * Messages write: GER_MODEL.
 */
void GUTS::setModel( const string model ) {

	/*
	 * Local copy (to edit), lowercase.
	 */
	string hmodel = model;
	std::transform(hmodel.begin(), hmodel.end(), hmodel.begin(), ::tolower);

	/*
	 * If model changes:
	 * If model doesn’t change, these check are already passed.
	 * Default model is "proper".
	 */
	if ( hmodel != mmodel ) {

		/*
		 * New messages on model.
		 * Else assign and make experiment.
		 */
		if ( hmodel != "proper" && hmodel != "it" ) {
			mMessages.at( GER_MODEL )  = true; // "Unknown model (valid: proper, it).";
		} else {
			mMessages.at( GER_MODEL )  =  false;
			mmodel = hmodel;
			makeExperiment();
		}

	} // End of if ( hmodel != mmodel )

} // End of GUTS::setModel( model ).





/**
 * Method: makeExperiment().
 * dist or model changed.
 * Set mExperiment and mnpar according to dist and model.
 * Invoke setParameters( mpar ) to provoke reseting of parameters, for messages.
 */
void GUTS::makeExperiment() {

	/*
	 * Default wpar { 0.0, 0.0, DBL_MAX_EXP, 0.0, 0.0 }.
	 * Default mparPos { 0, 1, 2, 3, 4 }. Resize needed, since sizes changes!
	 */
	wpar.assign( 5, 0.0 );
	wpar.at(2) = GMAX_DOUBLE; // represent Inf with largest double
	mparPos.assign( 5, 0 ); // five zeros, fill rest with loop.
	for ( int i = 1; i < 5; i++ ) {
		mparPos.at(i) = i;
	}

	/*
	 * Distribution.
	 */
	if ( mdist == "lognormal" ) {
		mExperiment = 11;
		//mparPos = {0, 1, 2, 3, 4};
	} else if ( mdist == "delta" ) {
		mExperiment = 21;
		//mparPos = {0, 1, 2, 3};
		//mparPos.erase(mparPos.begin()+4); // remove 5th parameter for delta (sd).
		mparPos.resize( 4 );
	} else if ( mdist == "empirical" ) {
		mExperiment = 91;
		//mparPos = {0, 1, 2};
		//mparPos.erase(mparPos.begin()+3, mparPos.end()); // remove 4th and 5th
		mparPos.resize( 3 );
	}

	/*
	 * Model.
	 */
	if ( mmodel == "it" ) {
		mExperiment = mExperiment + 1;
		mparPos.erase(mparPos.begin()+2); // remove third parameter
		//mnpar = hnpar - 1;
	}

	/*
	 * Testing only.
	 */
	if ( mdist == "lognormold" ) {
		mmodel = "proper";
		mExperiment = 15;
		mparPos.assign( 5, 0 ); // five zeros, fill rest with loop.
		for ( int i = 1; i < 5; i++ ) {
			mparPos.at(i) = i;
		}
	}

	/*
	 * Provoke set parameters with current mpar, for errors.
	 */
	setParameters( mpar );

} // End of GUTS::makeExperiment().





/**
 * Method setParameters:
 * * h_b: background mortality rate
 * * k_r: recovery rate
 * * k_k: killing rate
 * * ...: mean and sd of the lognormal distribution, if needed.
 * Messages write: GER_PAR, GWA_PAR
 */
void GUTS::setParameters( const vector<double>& par ) {

	/*
	 * Assign everything the user provides.
	 */
	mpar = par;

	/*
	 * New errors on par; warning on par if par < 0.
	 * If, return immediately with error/warning set.
	 * Else, reset errors on par.
	 */
	if ( par.size() != mparPos.size() ) {
		// Error.
		mMessages.at( GER_PAR ) = true; // "Vector of parameters must be of length " + std::to_string(mparPos.size()) + ".";
		return;
	} else if ( *min_element( par.begin(), par.begin()+mparPos.size() ) < 0.0 ) {
		// Warning.
		mMessages.at( GER_PAR ) = false; // this is no error, but a warning
		mMessages.at( GWA_PAR ) = true; // "Parameters must be non-negative values.";
		return;
	} else {
		/* Assign to wpar at positions in mparPos. */
		for ( unsigned int i = 0; i < mparPos.size(); ++i ) {
			wpar.at(mparPos.at(i)) = par.at(i);
		}
		/* Check for dist mean and sd. */
		if ( wpar.at(3) == 0.0 && wpar.at(4) != 0.0 ) {
			mMessages.at( GER_PAR ) = false;
			mMessages.at( GWA_PAR ) = true; // "Sample sd = 0 only if m = 0.";
			return;
		} else if ( std::isinf(wpar.at(2)) ) {
			/* Check for inf at 3rd par. */
			wpar.at(2) = GMAX_DOUBLE;
		}
		mMessages.at( GER_PAR ) = false;
		mMessages.at( GWA_PAR ) = false;
	}

} // End of GUTS::setParameters( par ).





/*
 * mN.
 * Length of mz.
 * Length of mzw.
 * Messages write: GER_N
 */
void GUTS::setSampleLength( const int& N ) {
	
	/*
	 * New messages on N.
	 * Else reset and assign.
	 */
	if ( N < 3 ) {
		mMessages.at(GER_N)  = true; // "N must be greater than 2.";
		return;
	} else {
		mMessages.at( GER_N )  = false;
		mN = N;

		ee.assign( mN, 0.0 );
		ff.assign( mN, 0 );

		makeDefaultGrid();
		mz.resize( mN );
		mzw.resize( mN );
	}

} // End of GUTS::setSampleLength()





/*
 * mM.
 * Messages write: GER_M
 */
void GUTS::setTimeGridPoints( const int& M ) {

	/*
	 * New message on M.
	 * Else reset messages, assign and resize D.
	 */
	if ( M < 2 ) {
		mMessages.at( GER_M )  = true; // "M must be greater than 1.";
		// Do not reset D, if M has not changed, because M has default.
		return;
	} else {
		mMessages.at( GER_M )  = false;
		mM = M;
		D.assign( mM, 0.0 );
	}

} // End of GUTS::setTimeGridPoints()





/**
 * setSample( … ).
 * Sample Vector Implementation.
 * Set mz, mzw.
 * Messages write: GER_EMP.
 */
void GUTS::setSample( const vector<double>& z ) {

	/*
	 * New messages on z.
	 * Else reset. Assign. Sort. Set weights.
	 */
	if ( z.size() < 3 ) {
		mMessages.at( GER_EMP )  = true; // "Sample must be longer than 2,";
		return;
	} else if ( *min_element( z.begin(), z.end() ) < 0.0 ) {
		mMessages.at( GER_EMP )  = true; // "Sample must contain non-negative values,";
		return;
	} else {
		mMessages.at( GER_EMP )  = false;
		
		mz = z;
		sort( mz.begin(), mz.end() );
		mzw.assign( mz.size(), 0.0 ); // Set all weights to 0.0.

		setDistribution( "empirical" );
		setSampleLength( mz.size() ); // Set N, mN, makeDefaultGrid(), z/zw resizes.
	}

} // End of GUTS::setSample()





/*
 * TODO: mInterval.
 * TODO: For testing purposes only.
 */
void GUTS::setInterval( const double& Interval ) {
	if ( Interval > 1.0 ) {
		mInterval = Interval;
	}
} // End of GUTS::setInterval()





/* End of GUTS Class Implementation. */
