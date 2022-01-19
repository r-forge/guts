/**
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch, alexander singer@rifcon.de, oliver.jakoby@rifcon.de, dirk.nickisch@rifcon.de
 * License GPL-2
 * 2017-10-09 
 * updated: 2019-01-29
 * updated: 2021-11-30
 * updated: 2022-01-17 
 */


#ifndef RANDOM_DISTRIBUTIONS_H
#define RANDOM_DISTRIBUTIONS_H

#include <limits>
#include <vector>
#include <cmath>
#include <random>
#include <stdexcept>

#include "external_data.h"

struct distribution {
	virtual ~distribution() {}
	virtual double CDF(const double x) const = 0;
};

struct lognormal_parameters {
  lognormal_parameters() :
  mn(std::numeric_limits<double>::quiet_NaN()),
  sd(std::numeric_limits<double>::quiet_NaN())
  {}
  virtual ~lognormal_parameters() {}
  inline void set_threshold_mean(const double new_mn) {mn = new_mn;}
  inline void set_threshold_sd(const double new_sd) {sd = new_sd;}
  inline double get_threshold_mean() const {return mn;}
  inline double get_threshold_sd() const {return sd;}
    protected:
    double mn;
    double sd;
};

struct lognormal :
		virtual public distribution,
		virtual public lognormal_parameters
		{
	virtual ~lognormal() {}
	inline double CDF(const double x) const final {
	  double sigma_square = std::log((sd * sd) / mn / mn + 1);
	  double mu = std::log(mn) - sigma_square / 2;
		return 0.5 + std::erf( (std::log(x)-mu)/std::sqrt(2 * sigma_square) ) / 2;
	}
};

class loglogistic_parameters {
public:
  loglogistic_parameters() :
  alpha(std::numeric_limits<double>::quiet_NaN()),
  beta(std::numeric_limits<double>::quiet_NaN())
{}
  virtual ~loglogistic_parameters() {}
  inline void set_threshold_alpha(const double new_alpha) {alpha = new_alpha;}
  inline void set_threshold_beta(const double new_beta) {beta = new_beta;}
  inline double get_threshold_alpha() const {return alpha;}
  inline double get_threshold_beta() const {return beta;}
protected:
  double alpha;
  double beta;
};

struct loglogistic :
		virtual public distribution,
		virtual public loglogistic_parameters
		{
	virtual ~loglogistic() {}
	inline double CDF(const double x) const final {
		return 1/(1+std::pow(x/alpha,-beta));
	}
};

class delta_parameters {
public:
  delta_parameters() :
  z_val(std::numeric_limits<double>::quiet_NaN())
{}
  virtual ~delta_parameters() {}
  inline void set_threshold(const double threshold){
    z_val = threshold;
  }
  inline double get_threshold() const {return z_val;}
protected:
  double z_val;
};

#endif //TD_SAMPLERS_H
