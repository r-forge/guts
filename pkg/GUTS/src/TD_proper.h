/**
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch, alexander singer@rifcon.de, oliver.jakoby@rifcon.de, dirk.nickisch@rifcon.de
 * License GPL-2
 * 2017-10-09 
 * updated: 2019-01-29
 * updated: 2021-11-30 
 * updated: 2022-01-17
 */


#ifndef TD_PROPER_H
#define TD_PROPER_H

#include <vector>

#include <iostream>
#include <algorithm>

#include "TD_base.h"
#include "samplers.h"
/**
 * @class abstract TD interface
 * 
 * @brief accumulates damage above the threshold and executes respective mortality
 */
template< typename sampler >
class TD_proper_base : public TD_base {
public:
	TD_proper_base() : TD_base(), samp(), ee(), ff(), zpos(0),
	kk(std::numeric_limits<double>::quiet_NaN()),
	dtau(std::numeric_limits<double>::quiet_NaN()),
	kkXdtau(std::numeric_limits<double>::quiet_NaN()),
	hb(std::numeric_limits<double>::quiet_NaN())
{}
	virtual ~TD_proper_base() {}
	bool is_still_gathering() const override {return true;}
	inline void set_killing_rate(const double new_kk) {
		kkXdtau = new_kk * dtau;
		kk = new_kk;
	}
	inline void set_background_mortality(const double new_hb) {hb = new_hb;}
	inline double get_killing_rate() const {return kk;}
	inline double get_background_mortality() const {return hb;}
	void update_to_next_survival_measurement() const override {}
	/**
	 * @brief gather an effect from known damage
	 * @param[in] D damage
	 */
	inline void gather_effect(const double D) const override {
		if ( D > samp.variate_back() ) {
			// damage higher than the largest value in threshold distribution
			ee.back() += D;
			ff.back() ++;
			return;
		}
		if ( D > samp.variate_at(0) ) {
			// damage within threshold distribution
			// Search quantile in threshold distribution that covers the damage.
			while ( zpos > 0 && D < samp.variate_at(zpos)) {
				--zpos;
			}
			while ( zpos < (samp.sample_size() - 1) && D > samp.variate_at(zpos) ) {
				++zpos;
			}
			ee.at(zpos-1) += D;
			ff.at(zpos-1)++;
		}
	}

	inline void set_start_conditions() const override {
		std::fill(ee.begin(), ee.end(), 0.0);
		std::fill(ff.begin(), ff.end(), 0);
		zpos = samp.sample_size()/2;
	}
protected:
	void initialize_threshold_distribution(const std::size_t sample_size) {
		ee.assign(sample_size, 0.0);
		ff.assign(sample_size, 0);
	}
	void initialize_time_discretization(const double new_dtau) {
		dtau = new_dtau;
	}
public:
	///the sampler
	mutable sampler samp;
protected:
	///brief gathered damage
	mutable std::vector<double > ee;
	///brief frequency distribution of damage == threshold
	mutable std::vector<unsigned > ff;
	mutable std::size_t zpos;
	///killing rate
	double kk;
	///length discrete time step
	double dtau;
	///killing rate times discrete time step
	double kkXdtau;
	///background mortality
	double hb;
};

template<typename sampler >
struct TD_proper_impsampling : public TD_proper_base<sampler > {
	TD_proper_impsampling() : TD_proper_base<sampler >() {}
	void set_start_conditions() const override {
		TD_proper_base<sampler >::set_start_conditions();
		this -> samp.calc_sample();
	}
	double calculate_current_survival(const double yt) const override {
		double E = 0.0;
		unsigned F = 0;
		double S = 0;
		std::size_t N = this->samp.sample_size();
		for (std::size_t u = N; u > 0; --u) {
			F += this->ff.at(u-1);
			E += this->ee.at(u-1);
			S += F == 0 ? exp(this->samp.weight_at(u-1) ) : exp((this->kkXdtau * (this->samp.variate_at(u-1) * F - E)) + this->samp.weight_at(u-1) );
		}
		return S * exp( -this->hb * yt ) / static_cast<double>(this->samp.sample_size());
	}
	virtual ~TD_proper_impsampling() {}
protected:
	void initialize_from_parameters() override {}
	inline void initialize(const double new_dtau, const std::size_t sample_size) {
		this -> samp.initialize(sample_size);
		TD_proper_impsampling<sampler >::initialize_threshold_distribution(sample_size);
		TD_proper_impsampling<sampler >::initialize_time_discretization(new_dtau);
	}
};

template< typename sampler >
class TD<sampler, 'P' > : public TD_proper_impsampling<sampler > {
public:
	TD() : TD_proper_impsampling<sampler >() {}
	template<typename tTDdata >
	inline void initialize(const tTDdata& TDdata) {
		this->samp.initialize(TDdata.N);
		this->initialize_threshold_distribution(TDdata.N);
		this->initialize_time_discretization(TDdata.calculate_dtau());
	}
	virtual ~TD() {}
};

template<>
struct TD<imp_delta, 'P' > : public TD_proper_impsampling<imp_delta > {
	TD() : TD_proper_impsampling<imp_delta >() {}
	template<typename tTDdata >
	inline void initialize(const tTDdata& TDdata) {
		this->samp.initialize();
		this->initialize_threshold_distribution(1);
		this->initialize_time_discretization(TDdata.calculate_dtau());
	}
	inline void set_threshold(const double new_z) {samp.set_threshold(new_z);}
	inline double get_threshold() const {return samp.get_threshold();}
	virtual ~TD() {}
};

template<typename tz >
class TD<random_sample<tz >, 'P' >: public TD_proper_base<random_sample<tz > > {
public:
	TD() : TD_proper_base<random_sample<tz > >() {  }
	template<typename tTDdata >
	inline void initialize(const tTDdata& TDdata) {
		TD_proper_base<random_sample<tz > >::initialize_time_discretization(TDdata.calculate_dtau());
	}
	void initialize_from_parameters() override {
		TD_proper_base<random_sample<tz > >::initialize_threshold_distribution(
				this->samp.sample_size()
		);
	}
	virtual ~TD() {}
	inline double calculate_current_survival(const double yt) const override {
		double E = 0.0;
		unsigned F = 0;
		double S = 1;
		std::size_t N = this->samp.sample_size();
		for (std::size_t u = N; u > 0; --u) {
			E += this->ee.at(u - 1);
			F += this->ff.at(u - 1);
			S += exp(   (this->kkXdtau * (this->samp.variate_at(u - 1) * F - E)) );
		}
		return S * exp( -this->hb * yt ) / static_cast<double>(N);
	}
};
#endif //TD_PROPER_H
