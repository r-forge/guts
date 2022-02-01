/**
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch, alexander singer@rifcon.de, oliver.jakoby@rifcon.de, dirk.nickisch@rifcon.de
 * License GPL-2
 * 2017-10-09 
 * updated: 2019-01-29
 * updated: 2021-11-30 
 * updated: 2022-01-17
 * updated: 2022-02-01
 */

#ifndef TD_IT_H
#define TD_IT_H

#include <vector>
#include <numeric>
#include <iostream>

#include "TD_base.h"
#include "samplers.h"
/**
 * @class TD IT
 * 
 * @brief accumulates damage above the threshold and executes respective mortality
 */
template< typename sampler >
class TD_IT_base : public TD_base, public background_mortality {
public:
	void initialize_from_parameters() override {}
  virtual ~TD_IT_base() {}
  bool is_still_gathering() const override {return zit != samp.end();}
  void update_to_next_survival_measurement() const override {
    // zit not reset, as lowest z above D can only increase over time
  }
    /**
     *\brief gather an effect from known damage
        * \param[in] D damage
        */
        inline void gather_effect(const double D) const override {
          zit = std::lower_bound(this->zit, samp.end(), D);
          // zit points to lowest z >= D
        }
  public:
  	mutable sampler samp;
  protected:
    ///lowest threshold with $z[zpos] \geq D$
    mutable typename sampler::sample_type::const_iterator zit;
};

double sumExp(const double& a, const double& b) {return a + std::exp(b);}

template< typename sampler >
class TD<sampler, 'I' > : public TD_IT_base<sampler > {
public:
  TD() : TD_IT_base<sampler >() {}
  virtual ~TD() {}
  inline virtual void initialize(const std::size_t sample_size) {
		this->samp.initialize(sample_size);
    Sj.resize(sample_size);
  }
  template<typename tTDdata >
  inline void initialize(const tTDdata& TDdata) {
    initialize(TDdata.N);
  }
  inline void set_start_conditions() const override {
  	this->samp.calc_sample();
  	this->zit = this->samp.z.begin();
	  std::vector<double >::reverse_iterator it = Sj.rbegin();
	  this->zw.back() = exp(this->zw.back());
	  std::partial_sum(this->zw.rbegin(), this->zw.rend(), it, sumExp);
  }
  inline double calculate_current_survival(const double yt) const override {
    return this->zit == this->z.end() ? 0 : Sj.at(this->zit - this->z.begin())  / this->samp.sample_size() * exp( -this->hb * yt );
  }
private:
  mutable std::vector<double > Sj;
};


struct TD_IT_CDF :
		virtual public TD_base,
		virtual public background_mortality {
	virtual ~TD_IT_CDF() {}
	  template<typename tTDdata >
	  inline void initialize(const tTDdata&) {M=0;}
	  void initialize_from_parameters() override {}
	  inline void set_start_conditions() const override {M=0;}
	  inline bool is_still_gathering() const override {return M<1;}
	  inline void update_to_next_survival_measurement() const override {};
	  inline double calculate_current_survival(const double yt) const override {
	    return (1-M) * std::exp( -this->hb * yt );
	  }
protected:
  mutable double M;
};

template<>
class TD<loglogistic, 'I' > :
public TD_IT_CDF
{
public:
  virtual ~TD() {}
  inline void gather_effect(const double D) const override {
	M = std::max(M, samp.CDF(D));
  }
	loglogistic samp;
};

template<>
class TD<lognormal, 'I' > :
public TD_IT_CDF
{
public:
  virtual ~TD() {}
  inline void gather_effect(const double D) const override {
	M = std::max(M, samp.CDF(D));
  }
	lognormal samp;
};


template<typename tz >
struct TD<random_sample<tz >, 'I' >: public TD_IT_base<random_sample<tz > > {
	TD() : TD_IT_base<random_sample<tz > >() {}
	virtual ~TD() {}
	template<typename tTDdata > inline void initialize([[gnu::unused]] const tTDdata& TDdata) {}
	void set_start_conditions() const override {
		this->zit = this->samp.begin();
	}
  inline double calculate_current_survival(const double yt) const override {
    std::size_t S = this->samp.end() - this->zit;
    return static_cast<double>(S) * exp( -this->hb * yt ) / 
    	static_cast<double>(this->samp.sample_size());
  }
};

#endif //TD_IT_H
