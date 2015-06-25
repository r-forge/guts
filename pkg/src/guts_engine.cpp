/**
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 * Function guts_engine(GUTS-Object, Parameters).
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch
 * License GPL-2
 * 2015-06-25
 */

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void guts_engine( Rcpp::List gobj, Rcpp::NumericVector par ) {

	/*
	 * Check if this is an object of class GUTS.
	 */
	if ( !gobj.inherits("GUTS") ) {
		Rcpp::stop( "No GUTS object. Use `guts_setup()` to create or modify objects." );
	}

	/*
	 * Get object elements that are necessary until here.
	 */
	NumericVector       par_pos  =  gobj.attr("par_pos");
	NumericVector       wpar     =  gobj.attr("wpar");
	NumericVector       opar     =  gobj["par"];
	NumericVector       oS       =  gobj["S"];
	std::vector<double> par_NA(opar.size(), NA_REAL);
	std::vector<double> S_NA(oS.size(), NA_REAL);

	/*
	 * Check parameters.
	 * Error: Stop and do nothing. Warning: Assign NAs.
	 * Assign working parameters.
	 * Assign parameters to object.
	 */
	if ( par.size() != par_pos.size() ) {
		// Error.
		Rcpp::stop( "Vector of parameters has wrong length." );
	} else if ( *std::min_element( par.begin(), par.end() ) < 0.0 ) {
		// Warning.
		Rcpp::warning( "Parameters must be non-negative values -- ignored." );
		gobj["par"] = par_NA;
		gobj["S"]   = S_NA;
		gobj["LL"]  =  -std::numeric_limits<double>::infinity();
		return;
	} else {
		// assign
		for ( int i = 0; i < par_pos.size(); ++i ) {
			wpar[(par_pos[i] - 1)] = par[i]; // R positions to C++ position.
		}
		// adjust wpar at 3.
		if ( std::isinf(wpar[2]) ) {
			wpar[2] = std::numeric_limits<double>::max();
		}
		// check for sd
		if ( wpar[3] == 0.0 && wpar[4] != 0.0 ) {
			Rcpp::warning( "Parameter sd = 0 only if mn = 0 -- ignored." );
			gobj["par"] = par_NA;
			gobj["S"]   = S_NA;
			gobj["LL"]  =  -std::numeric_limits<double>::infinity();
			return;
		}
		// assign if not returned earlier
		gobj["par"] = par;
		gobj.attr("wpar") = wpar;
		
	}

	/*
	 * Experiment type.
	 * Vectors.
	 * S.
	 * LL.
	 */
	std::vector<double> C  = gobj["C"];
	std::vector<double> Ct = gobj["Ct"];
	std::vector<int>    y  = gobj["y"];
	std::vector<double> yt = gobj["yt"];
	int N = gobj["N"];
	int M = gobj["M"];
	int experiment = gobj.attr("experiment");
	std::vector<double> S(yt.size(), 0.0);
	double LL = 0.0;

	/*
	 * Temporary variables.
	 */
	std::vector<double>  D( M, 0.0 );  // Damages.
	std::vector<double> ee( N, 0.0 );
	std::vector<int>    ff( N, 0 );
	double Scale = 0.0;          // To prevent warnings.
	double tmp;
	double summand3;
	double E;
	int    F;

	/*
	 * tau, dtau, tauit, wpar[2] * dtau.
	 */
	double tau     = 0.0;
	double dtau    = (yt.back() - yt.front()) / M;
	int    tauit   = 0;
	double kkXdtau = wpar[2] * dtau;

	/*
	 * Representations of z, zw.
	 */
	double sigma2  =  0.0;
	double mu      =  0.0;
	double sigmaD  =  0.0;
	std::vector<double> z(N);
	std::vector<double> zw(N);
	if ( experiment < 20 && wpar[4] != 0 ) {
		sigma2   =  log(   1.0  +  pow( (wpar[4] / wpar[3]), 2.0 )   );
		mu       =  log(wpar[3])  -  (0.5 * sigma2);
		sigmaD   =  sqrt(sigma2) * 4; // Interval fixed to 4.
		double ztmp;
		for ( unsigned int i = 0; i < N; ++i ) {
			ztmp = (2.0 * i - N + 1.0) / (N-1.0);
			z.at(i)  =  exp( ztmp * sigmaD + mu );
			zw.at(i) =  -0.5  *  pow((ztmp * sigmaD), 2.0)  /  sigma2;
		}
	} else {
		z.assign(N, wpar[3]);
		zw.assign(N, 0.0);
	}

	/* 
	 * Iterators, positions.
	 */
	int dpos    = 0;
	int ii      = 0;
	int zpos    = 0;
	int k       = 0;

	/*
	 * Diffs.
	 */
	std::vector<double>   diffS(yt.size());
	std::vector<int>      diffy(yt.size());
	std::vector<double>   diffC(C.size());
	std::vector<double>  diffCt(Ct.size());
	std::vector<double> diffCCt(Ct.size());
	for ( int i = 1; i < C.size(); ++i ) {
		  diffC.at(i-1) =  C.at(i) -  C.at(i-1);
		 diffCt.at(i-1) = Ct.at(i) - Ct.at(i-1);
		diffCCt.at(i-1) = diffC.at(i-1) / diffCt.at(i-1);
	}

	/*
	 * Loop over yt.
	 */
	for ( unsigned int ytpos = 0; ytpos < yt.size(); ++ytpos ) {

		while ( tau < yt.at(ytpos) && dpos < M ) {

			tmp = exp( -wpar[1] * (tau - Ct.at(k)) );
			if ( wpar[1] > 0.0 ) {
				summand3 = (tau - Ct.at(k) - (1.0-tmp)/wpar[1])  *  diffCCt.at(k);
			} else {
				summand3 = 0.0;
			}
			D.at(dpos) =  tmp * (D.at(ii) - C.at(k))  +  C.at(k)  +  summand3;

			/*
			 * Find zpos from D.at(i).
			 * Update ee and ff.
			 */
			while ( D.at(dpos) < z.at(zpos) && zpos > 0 ) {
				--zpos;
			}
			while ( D.at(dpos) > z.at(zpos) && zpos < (N-1) ) {
				++zpos;
			}
			if ( D.at(dpos) > z.at(N-1) ) {
				ee.at(N-1) += D.at(dpos);
				ff.at(N-1)++;
			} else if ( D.at(dpos) > z.at(0) ) {
				ee.at(zpos-1) += D.at(dpos);
				ff.at(zpos-1)++;
			}

			/*
			 * Increment or decrement dpos, tau, k.
			 */
			++dpos;
			tau = dtau * (++tauit);
			if ( tau > Ct[k+1] ) {
				++k;
				ii = dpos-1;
			}

		} // End of while ( tau < yt[ytpos] && dpos < M )

		/*
		 * Write S at ytpos.
		 */
		E = 0.0;
		F = 0;
		for ( int u=N-1; u >= 0; --u ) {
			E += ee.at(u);
			F += ff.at(u);
			S.at(ytpos) += exp(   (kkXdtau * (z.at(u) * F - E))  +  zw.at(u)   );
		}
		S.at(ytpos) *= exp( -wpar[0] * yt[ytpos] ) / N;

		/*
		 * Scale S.
		 * Calculate diffy, diffS.
		 */
		if ( ytpos == 0 && S.at(0) > 0.0 ) {
			Scale = 1.0 / S.at(0);
		} else {
			diffS.at(ytpos-1) = S.at(ytpos-1) - Scale * S.at(ytpos);
			diffy.at(ytpos-1) = y.at(ytpos-1) - y.at(ytpos);
		}
		S.at(ytpos) *= Scale;

	} // End of for ( unsigned int ytpos = 0; ytpos < yt.size(); ytpos++ )

	/*
	 * Last diffy, diffS.
	 */
	diffS.back() = S.back();
	diffy.back() = y.back();

	/*
	 * Update object fields.
	 */
	gobj["S"]   = S;

	/*
	 * Calculate Loglikelihood.
	 */
	for ( unsigned int i=0; i < diffy.size(); ++i ) {
		if ( diffy.at(i) > 0 ) {
			if ( diffS.at(i) == 0.0 ) {
				gobj["LL"] = -std::numeric_limits<double>::infinity();
				return;
			}
			LL += diffy.at(i) * log(diffS.at(i));
		}
	} // End of for ( unsigned int i=0; i < diffy.size(); ++i ).

	/*
	 * Update object fields.
	 */
	gobj["LL"]  = LL;

} // End of void guts_engine( Rcpp::List gobj, NumericVector par ).
