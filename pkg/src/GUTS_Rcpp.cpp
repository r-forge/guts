/** GUTS Rcpp module.
 * @file        GUTS_Rcpp.cpp
 * @author      soeren.vogel@uzh.ch
 * @author      carlo.albert@eawag.ch
 * @date        2012-05-03
 * @license     GPL-2
 *
 * Class GUTS in R (RGUTS)
 */

#ifndef guts_h
#include "GUTS.h"
#endif
#include <Rcpp.h>

using namespace Rcpp;

RCPP_MODULE(modguts)
{
    class_<GUTS>( "GUTS" )
        .constructor()
        
        .method( "setConcentrations",   &GUTS::setConcentrations, "Set time series vector of concentrations." )
        .method( "setSurvivors",        &GUTS::setSurvivors, "Set time series vector of survivors." )
        .method( "setParameters",       &GUTS::setParameters, "Set parameter vector of the object." )
        .method( "setTimeGridPoints",   &GUTS::setTimeGridPoints, "Set number of grid points on the time-axis." )
        .method( "setDistribution",     &GUTS::setDistribution, "Set distribution to sample from." )
        .method( "setSampleLength",     &GUTS::setSampleLength, "Set length of sample (numerical exactness)." )

        .method( "setSample",           &GUTS::setSample, "Set sample vector." )

        .method( "calcSample",          &GUTS::calcSample, "Calculate a sample." )
        .method( "calcSurvivalProbabilities", &GUTS::calcSurvivalProbabilities, "Calculate the survivor probabilities." )
        .method( "calcLoglikelihood",   &GUTS::calcLoglikelihood, "Calculate the logarithm of the likelihood." )

        .property( "C",     &GUTS::getC, "Vector of concentrations." )
        .property( "Ct",    &GUTS::getCt, "Vector of concentration times." )
        .property( "y",     &GUTS::gety, "Vector of survivors." )
        .property( "yt",    &GUTS::getyt, "Vector of survivors times." )
        .property( "par",   &GUTS::getpar, "Parameter vector." )
        .property( "M",     &GUTS::getM, "Grid points on time axis." )
        .property( "dist",  &GUTS::getdist, "Distribution." )
        .property( "N",     &GUTS::getN, "Sample length (numerical exactness)." )
        .property( "z",     &GUTS::getz, "Sample." )

        .property( "D",     &GUTS::getD, "Vector of damages." )
        .property( "S",     &GUTS::getS, "Survival probabilities." )
        .property( "LL",    &GUTS::getLL, "Loglikelihood." )

        .property( "errors",    &GUTS::getErrorMessages, "Errors." )
    ;
}
