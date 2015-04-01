/**
 * GUTS Rcpp module.
 *
 * @file    GUTS_Rcpp.cpp
 * @author  soeren.vogel@posteo.ch
 * @author  carlo.albert@eawag.ch
 * @date    2015-03-26
 * @license GPL-2
 */



#ifndef guts_h
#include "GUTS.h"
#endif
#include <Rcpp.h>



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

		.property( "Messages", &GUTS::getMessages, "Bool vector of messages." )
		;
}
