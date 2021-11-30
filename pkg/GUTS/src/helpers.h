/**
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch, alexander singer@rifcon.de, oliver.jakoby@rifcon.de, dirk.nickisch@rifcon.de
 * License GPL-2
 * 2021-11-30 
 */

#ifndef HELPERS_H
#define HELPERS_H

#include <Rcpp.h>

inline int back(const Rcpp::IntegerVector& vec) {return vec.at(vec.size()-1);}
inline int front(const Rcpp::IntegerVector& vec) {return vec.at(0);}
inline double back(const Rcpp::NumericVector& vec) {return vec.at(vec.size()-1);}
inline double front(const Rcpp::NumericVector& vec) {return vec.at(0);}
inline double back(const std::vector<double >& vec) {return vec.back();}
inline double front(const std::vector<double >& vec) {return vec.front();}

#endif