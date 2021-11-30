/**
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch, alexander singer@rifcon.de, oliver.jakoby@rifcon.de, dirk.nickisch@rifcon.de
 * License GPL-2
 * 2017-10-09 
 * updated: 2019-01-29
 * updated: 2021-11-30 
 */

#ifndef GUTS_BASE_H_
#define GUTS_BASE_H_

#include <cstddef>
#include <memory>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <functional>


#include "helpers.h"

/** 
 * \brief Defines the public combination of a TK and a TD object 
 * \tparam TK Type of TK model
 * \tparam TD Type of TD model
 * \tparam tData object with model relevant data
 * \detail This generic definition should be parent of any guts_model
 */
template<typename TK_mod, typename TD_mod >
struct guts_model: 
  virtual public TK_mod, virtual public TD_mod {
  virtual ~guts_model() {};
  template<typename tData >
  inline void initialize(const tData& data) {
    TK_mod::initialize(data);
    TD_mod::initialize(data);
  }
  void initialize_from_parameters() override {
  	TK_mod::initialize_from_parameters();
  	TD_mod::initialize_from_parameters();
  	}
  void set_start_conditions() override {
    TK_mod::set_start_conditions();
    TD_mod::set_start_conditions();
  }
  guts_model() : TK_mod(), TD_mod() {}
};

template<typename tModel, typename tt, typename tSurvival >
struct guts_projector_base : public tModel {
  typedef tSurvival tProjection;
  virtual ~guts_projector_base() {}
  inline void set_start_conditions() override {
  	tModel::set_start_conditions();
  }
  inline void get_survival_projection(tProjection& proj) const {proj = p;}
  void project_survival () const {
    p.assign(yt->size(), 0);
    
    p.at(0) = tModel::TD_mod::calculate_current_survival(0);
    if ( p.at(0) <= 0.0 ) {
      // should never happen with well defined parameters
      throw std::underflow_error("Numeric underflow: Survival cannot be calculated for given parameter values." );
    }
    auto ytpos = 1; //index yt
    while (ytpos < yt->size() && p.at(ytpos-1) > 0) {
      tModel::TD_mod::update_to_next_survival_measurement();
      gather_effect_per_time_step(yt->at(ytpos));
      p.at(ytpos) = tModel::TD_mod::calculate_current_survival(yt->at(ytpos)) / p.at(0);
      ++ytpos;
    }
    p.at(0) = 1;
  }
  virtual std::vector<double > get_damage() const = 0;
  virtual std::vector<double > get_damage_time() const = 0;
  template<typename tData >
  inline void initialize(const tData& data) {
    yt = data.yt;
    p.assign(yt->size(), std::numeric_limits<double>::quiet_NaN());
    tModel::initialize(data);
  }
  //tData exp_dat;
protected:
  std::shared_ptr<const tt > yt;
  virtual void gather_effect_per_time_step(const double) const = 0;
private:
  mutable tSurvival p;
};

template<typename tModel, typename tt, typename tSurvival >
struct guts_projector: public guts_projector_base<tModel, tt, tSurvival > {
public:
	typedef tSurvival tProjection;
	typedef guts_projector_base<tModel, tt, tSurvival > parent;
	virtual ~guts_projector() {}
	template<typename tData >
	inline void initialize(const tData& data) {
	  M = data.M;
	  dtau = data.calculate_dtau(); 
	  parent::initialize(data);
	}
	inline virtual void set_start_conditions() override {
		tauit = 0; //index discrete time
		k = 0;     //index Ct
		D.assign(M, std::numeric_limits<double>::quiet_NaN());
		parent::set_start_conditions();
	}
	std::vector<double > get_damage() const override {return D;}
	std::vector<double > get_damage_time() const override {
		std::vector<double > damage_time(M, std::numeric_limits<double>::quiet_NaN());
		damage_time[0] = 0;
		for (
				std::vector<double >::iterator it = damage_time.begin() + 1; 
      it != damage_time.begin() + tauit;
      ++it) { 
       *it = *(it-1) + dtau;
		}
		return damage_time;
	}
protected:
	std::size_t M;
	double dtau;
	mutable std::vector<double > D;
private:
	mutable std::size_t tauit; //index discrete time
	mutable std::size_t k;     //index Ct
	void gather_effect_per_time_step (const double yt) const override {
		double tau = dtau * tauit;		 //discrete absolute time
		while ( tauit < M && tau < yt && tModel::TD_mod::is_still_gathering() ) {
			D.at(tauit) = tModel::TK_mod::calculate_damage(k, tau);
			tModel::TD_mod::gather_effect(D[tauit]);
			tau = dtau * (++tauit);
			if (tau > tModel::TK_mod::Ct->at(k+1)) {
				++k; // concentration index
				tModel::TK_mod::update_to_next_concentration_measurement();
			}
		}
	}
};

template<typename tModel, typename tt, typename tSurvival >
struct guts_projector_fastIT: 
  public guts_projector_base<tModel, tt, tSurvival > {
public:
	typedef tSurvival tProjection;
	typedef guts_projector_base<tModel, tt, tSurvival > parent;
	virtual ~guts_projector_fastIT() {}
	inline void set_start_conditions() override {
		k = 0;
		Dk = 0;
		damage.resize(0);
		damage_time.resize(0);
		damage_time.push_back(0);
		damage.push_back(0);
		parent::set_start_conditions();
	}
	std::vector<double > get_damage() const override {
		// ensure that the function is not called repeatedly. 
		// Note: survival calculations automatically increase Dk.
		if (Dk != 0) {
			tModel::TK_mod::set_start_conditions();
			extend_damage_values();
		}
		return damage;
		}
	std::vector<double > get_damage_time() const override {
		// ensure that the function is not called repetedly. 
		// Note: survival calculations automatically increase Dk. 
		if (Dk != 0) {
			tModel::TK_mod::set_start_conditions();
			extend_damage_values();
		}
		return damage_time;
	}
	
private:
	mutable std::size_t k;
	mutable std::size_t Dk;
	mutable std::vector<double > damage_time;
	mutable std::vector<double > damage;
	void gather_effect_per_time_step (const double yt) const override {
		double te;
		std::size_t Dk_old = Dk;
		while (this->Ct->at(k+1) < yt && this->is_still_gathering() ) {
			// check damage at local maximum
			if (this->is_maximum_damage(k)) {
				//a maximum exists (at an extreme point)
				te = this->calculate_time_of_extreme_damage(k);
				if (te < yt) {
					if (te < this->Ct->at(k+1)) {
						// the maximum is within the current concentration measurement interval
						damage_time.push_back(te);
						damage.push_back(this->calculate_damage(k, te));
						++Dk;
					}
				}
			}
		  // check damage at concentration measurement times (i.e. boundaries)
		  	damage_time.push_back(this->Ct->at(k+1));
		  	damage.push_back(this->calculate_damage(k, back(damage_time)));
        ++Dk;
        ++k;
        this->update_to_next_concentration_measurement();
		  }
		damage_time.push_back(yt);
		damage.push_back(this->calculate_damage(k, yt));
		++Dk;
		this->gather_effect( 
				*(std::max_element(damage.begin() + Dk_old, damage.end())) 
			);
	}
	
	void extend_damage_values(std::size_t num_extra_evals_per_time_interval = 10) const {
		double dtau;
		double max_time = *(std::max_element(damage_time.begin(), damage_time.end()));
		double cur_time;
		k = 0;
		Dk = 0;
		while (this->Ct->at(k) < max_time) {
				dtau = (this->Ct->at(k+1) - this->Ct->at(k)) / (num_extra_evals_per_time_interval);
				cur_time = this->Ct->at(k) + dtau;
				do {
					damage_time.push_back(cur_time);
					damage.push_back(this->calculate_damage(k, cur_time));
					cur_time += dtau;
				} while (cur_time < this->Ct->at(k+1) && cur_time < max_time);
				this->calculate_damage(k, this->Ct->at(k+1));
				++k;
				this->update_to_next_concentration_measurement();
		}
	}
};

template<typename tProjection, typename tmeasured_survivors >
  double calculate_loglikelihood(const tProjection& p, const tmeasured_survivors& y) {
    std::size_t diffy;
    double diffS;
    double loglik;
    if (back(y) > 0) {
      if (back(p) == 0.0) {
        return -std::numeric_limits<double >::infinity(); 
      } else {
        loglik = back(y) * std::log(back(p));
      }
    } else {
      loglik = 0;
    }
    for (auto i=1; i < y.size(); ++i ) {
      diffy = y.at(i-1) - y.at(i);
      if (diffy > 0) {
        diffS = p.at(i-1) - p.at(i);
        if (diffS == 0.0) {
          return -std::numeric_limits<double >::infinity();
        }
        loglik += diffy * std::log(diffS);
      }
    } 
    return loglik;
  }

template<typename tProjection, typename tmeasured_survivors >
  double calculate_SPPE(const tProjection& p, const tmeasured_survivors& y) {
    return (static_cast<double>(back(y)) / static_cast<double>(front(y)) - back(p)) * 100.0;
  }

template<typename tProjection, typename tmeasured_survivors >
  double calculate_sum_of_squares(const tProjection& p, const tmeasured_survivors& y) {
    double sum_of_squares = 0.0;
    double diff;
    
    double y0 = static_cast<double>(front(y));
    for (auto i=0; i < y.size(); ++i ) {
      diff = static_cast<double>(y.at(i)) - y0 * p.at(i);
      sum_of_squares += diff * diff;
    }
    
    return sum_of_squares;
  }

#endif //GUTS_BASE_H_
