/**
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch, alexander singer@rifcon.de, oliver.jakoby@rifcon.de, dirk.nickisch@rifcon.de
 * License GPL-2
 * 2017-10-09 
 * updated: 2019-01-29
 * updated: 2021-11-30 
 */

#ifndef SAMPLERS_H
#define SAMPLERS_H

#include <limits>
#include <vector>
#include <cmath>
#include <stdexcept>

#include "random_distributions.h"

class importance_sampler {
public:
	typedef std::vector<double > sample_type;
  importance_sampler(const std::size_t sample_size = 0) : z(sample_size), zw(sample_size) {}
  virtual ~importance_sampler() {}
  inline virtual void initialize(const std::size_t sample_size) {
	  z.assign(sample_size, 0.0);
	  zw.assign(sample_size, 0.0);
  }
  virtual void calc_sample() = 0;
  inline double variate_at(const size_t i) const {return z.at(i);}
  inline double weight_at(const size_t i) const {return zw.at(i);}
  inline double variate_back() const {return z.back();}
  inline std::size_t sample_size() const {return z.size();}
protected:
  std::vector<double > z; 
  std::vector<double > zw;
};

class imp_lognormal : public importance_sampler, public lognormal_parameters {
public:
  imp_lognormal(
    const std::size_t sample_size = 0, 
    const double importance_sampling_rate = 4.0
  ) : 
  importance_sampler(sample_size),
  lognormal_parameters(),
  R(importance_sampling_rate) {}
  virtual ~imp_lognormal() {}
  virtual void calc_sample();
    protected:
    double R;
};

class imp_loglogistic : public importance_sampler, public loglogistic_parameters {
public:
  imp_loglogistic(
    const std::size_t sample_size = 0, 
    const double importance_sampling_rate = 50.0
  ) : 
  importance_sampler(sample_size),
  loglogistic_parameters(),
  R(importance_sampling_rate) {}
  virtual ~imp_loglogistic() {}
  virtual void calc_sample();
protected:
  double R;
};

class imp_delta : public importance_sampler, public delta_parameters {
public:
  imp_delta() :
	  importance_sampler(1),
	  delta_parameters()
  {}
  virtual ~imp_delta() {}
  inline virtual void initialize() {importance_sampler::initialize(1);}
  void calc_sample() override;
};

template<typename tz >
class random_sample  {
public:
	typedef tz sample_type;
  random_sample() : z() {}
  virtual ~random_sample() {}
  inline double variate_at(const size_t i) const {return z.at(i);}
  inline void set_variates(const tz& variates) {z = variates;}
	void set_variates(typename tz::const_iterator begin, typename tz::const_iterator end) {
		z.assign(begin, end);
	}
  tz get_variates() const {return z;}
  double variate_back() const {return *(z.end()-1);}
  std::size_t sample_size() const {return z.size();}
protected:
  tz z;
};

#endif //SAMPLERS_H
