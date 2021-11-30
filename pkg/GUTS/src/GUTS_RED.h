/**
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch, alexander singer@rifcon.de, oliver.jakoby@rifcon.de, dirk.nickisch@rifcon.de
 * License GPL-2
 * 2017-10-09 
 * updated: 2019-01-29
 * updated: 2021-11-30 
 */

#ifndef GUTS_RED_H
#define GUTS_RED_H

//#include <typeinfo>

#include "GUTS_base.h"
#include "TK_RED.h"
#include "TD.h"
#include "helpers.h"


template<typename tt, typename tc, typename TD_mod, typename tparam >
struct guts_RED_base :
  public guts_model<TK_RED<tt, tc >, TD_mod >
{
  typedef TK_RED<tt, tc > TK_mod;
  enum class position : std::size_t {hb = 0, kd = 1, kk = 2, t1 = 3, t2 = 4};
  
  virtual tparam  get_parameters() const = 0;
  virtual void set_parameters(const tparam& param) = 0;
};

template<typename tparam, typename parameterized_model >
void get_parameters_hb_kd(const parameterized_model& model, tparam& param = tparam(5)) {
  param[static_cast<std::size_t >(parameterized_model::position::hb)] = model.get_background_mortality();
  param[static_cast<std::size_t >(parameterized_model::position::kd)] = model.get_dominant_rate_constant();
}
template<typename tparam, typename parameterized_model >
void set_parameters_hb_kd(parameterized_model& model, const tparam& param) {
  model.set_background_mortality(param[static_cast<std::size_t >(parameterized_model::position::hb)]);
  model.set_dominant_rate_constant(param[static_cast<std::size_t >(parameterized_model::position::kd)]);
}
template<typename tparam, typename parameterized_model >
void get_parameters_kk(const parameterized_model& model, tparam& param = tparam(5)) {
  param[static_cast<std::size_t >(parameterized_model::position::kk)] = model.get_killing_rate();
}
template<typename tparam, typename parameterized_model >
void set_parameters_kk(parameterized_model& model, const tparam& param) {
  model.set_killing_rate(param[static_cast<std::size_t >(parameterized_model::position::kk)]);
}
template<typename tparam, typename parameterized_model >
void get_parameters_z(const parameterized_model& model, tparam& param = tparam(5)) {
  param[static_cast<std::size_t >(parameterized_model::position::t1)] = model.get_threshold();
}
template<typename tparam, typename parameterized_model >
void set_parameters_z(parameterized_model& model, const tparam& param) {
  model.set_threshold(param[static_cast<std::size_t >(parameterized_model::position::t1)]);
}
template<typename tparam, typename parameterized_model >
void get_lognormal_threshold_parameters(const parameterized_model& model, tparam& param) {
  param[static_cast<std::size_t >(parameterized_model::position::t1)] = model.get_threshold_mean();
  param[static_cast<std::size_t >(parameterized_model::position::t2)] = model.get_threshold_sd();
}
template<typename tparam, typename parameterized_model >
void set_lognormal_threshold_parameters(parameterized_model& model, const tparam& param) {
  model.set_threshold_mean(param[static_cast<std::size_t >(parameterized_model::position::t1)]);
  model.set_threshold_sd(param[static_cast<std::size_t >(parameterized_model::position::t2)]);
}
template<typename tparam, typename parameterized_model >
void get_loglogistic_threshold_parameters(const parameterized_model& model, tparam& param) {
  param[static_cast<std::size_t >(parameterized_model::position::t1)] = model.get_threshold_alpha();
  param[static_cast<std::size_t >(parameterized_model::position::t2)] = model.get_threshold_beta();
}
template<typename tparam, typename parameterized_model >
void set_loglogistic_threshold_parameters(parameterized_model& model, const tparam& param) {
  model.set_threshold_alpha(param[static_cast<std::size_t >(parameterized_model::position::t1)]);
  model.set_threshold_beta(param[static_cast<std::size_t >(parameterized_model::position::t2)]);
}

template<typename tt, typename tc, typename TD_mod, typename tparam >
struct guts_RED
{
  guts_RED() = delete;
};

template<typename tt, typename tC, typename tparam >
struct guts_RED<tt, tC, TD_proper_lognormal, tparam > :
  public guts_RED_base<tt, tC, TD_proper_lognormal, tparam  > {
  typedef TD_proper_lognormal TD_mod;
  tparam get_parameters() const override {
    tparam param(5);
    get_parameters_hb_kd(*this, param);
    get_parameters_kk(*this, param);
    get_lognormal_threshold_parameters(*this, param);
    return param;
  }
  void set_parameters(const tparam& param) override {
    set_parameters_hb_kd(*this, param);
    set_parameters_kk(*this, param);
    set_lognormal_threshold_parameters(*this, param);
  }
};

template<typename tt, typename tC, typename tparam >
struct guts_RED<tt, tC, TD_proper_loglogistic, tparam > :
  public guts_RED_base<tt, tC, TD_proper_loglogistic, tparam  > {
  typedef TD_proper_loglogistic TD_mod;
  tparam get_parameters() const override {
    tparam param(5);
    get_parameters_hb_kd(*this, param);
    get_parameters_kk(*this, param);
    get_loglogistic_threshold_parameters(*this, param);
    return param;
  }
  void set_parameters(const tparam& param) override {
    set_parameters_hb_kd(*this, param);
    set_parameters_kk(*this, param);
    set_loglogistic_threshold_parameters(*this, param);
  }
};

template<typename tt, typename tC, typename tparam >
struct guts_RED<tt, tC, TD_proper_delta, tparam  > :
  public guts_RED_base<tt, tC, TD_proper_delta, tparam  > {
  typedef TD_proper_delta TD_mod;
  tparam get_parameters() const override {
    tparam param(5);
    get_parameters_hb_kd(*this, param);
    get_parameters_kk(*this, param);
    get_parameters_z(*this, param);
    return param;
  }
  void set_parameters(const tparam& param) override {
    set_parameters_hb_kd(*this, param);
    set_parameters_kk(*this, param);
    set_parameters_z(*this, param);
  }
};

template<typename tt, typename tC, typename tparam >
struct guts_RED<tt, tC, TD<random_sample<tparam >, 'P' >, tparam > :
	public guts_RED_base<tt, tC, TD<random_sample<tparam >, 'P' >, tparam  > {
	typedef TD<random_sample<tparam >, 'P' > TD_mod;
	tparam get_parameters() const override {
		tparam param(3);
		get_parameters_hb_kd(*this, param);
		get_parameters_kk(*this, param);
		return param;
	}
	void set_parameters(const tparam& param) override {
		set_parameters_hb_kd(*this, param);
		set_parameters_kk(*this, param);
		TD_mod::set_variates(param.begin() + 3, param.end());
	}
};

template<typename tt, typename tC, typename tparam >
struct guts_RED<tt, tC, TD_SD, tparam  > :
  public guts_RED_base<tt, tC, TD_SD, tparam  > {
  typedef TD_SD TD_mod;
  tparam get_parameters() const override {
    tparam param(5);
    get_parameters_hb_kd(*this, param);
    get_parameters_kk(*this, param);
    get_parameters_z(*this, param);
    return param;
  }
  void set_parameters(const tparam& param) override {
    set_parameters_hb_kd(*this, param);
    set_parameters_kk(*this, param);
    set_parameters_z(*this, param);
  }
};

template<typename tt, typename tC, typename lognormal_sampler, typename tparam >
struct guts_RED_IT_lognormal :
public guts_RED_base<tt, tC, TD<lognormal_sampler, 'I' >, tparam  > {
  typedef TD<lognormal_sampler, 'I' > TD_mod;
  tparam get_parameters() const override {
    tparam param(5);
    get_parameters_hb_kd(*this, param);
    get_lognormal_threshold_parameters(*this, param);
    return param;
  }
  void set_parameters(const tparam& param) override {
    set_parameters_hb_kd(*this, param);
    set_lognormal_threshold_parameters(*this, param);
  }
};

template<typename tt, typename tC, typename tparam >
struct guts_RED<tt, tC, TD_IT_lognormal, tparam > :
  public guts_RED_IT_lognormal<tt, tC, lognormal, tparam  > {
};

template<typename tt, typename tC, typename tparam >
struct guts_RED<tt, tC, TD_IT_imp_lognormal, tparam > :
  public guts_RED_IT_lognormal<tt, tC, imp_lognormal, tparam  > {
};

template<typename tt, typename tC, typename loglogistic_sampler, typename tparam >
struct guts_RED_IT_loglogistic :
  public guts_RED_base<tt, tC, TD<loglogistic_sampler, 'I' >, tparam  > {
  typedef TD<loglogistic_sampler, 'I' > TD_mod;
  tparam get_parameters() const override {
    tparam param;
    get_parameters_hb_kd(*this, param);
    get_loglogistic_threshold_parameters(*this, param);
    return param;
  }
  void set_parameters(const tparam& param) override {
    set_parameters_hb_kd(*this, param);
    set_loglogistic_threshold_parameters(*this, param);
  }
};

template<typename tt, typename tC, typename tparam >
struct guts_RED<tt, tC, TD_IT_loglogistic, tparam > :
  public guts_RED_IT_loglogistic<tt, tC, loglogistic, tparam  > {
};

template<typename tt, typename tC, typename tparam >
struct guts_RED<tt, tC, TD_IT_imp_loglogistic, tparam > :
  public guts_RED_IT_loglogistic<tt, tC, imp_loglogistic, tparam  > {
};

template<typename tt, typename tC, typename tparam >
struct guts_RED<tt, tC, TD<random_sample<tparam >, 'I' >, tparam > :
	public guts_RED_base<tt, tC, TD<random_sample<tparam >, 'I' >, tparam  > {
	typedef TD<random_sample<tparam >, 'I' > TD_mod;
	tparam get_parameters() const override {
		tparam param(2);
		get_parameters_hb_kd(*this, param);
		return param;
	}
	void set_parameters(const tparam& param) override {
		set_parameters_hb_kd(*this, param);
		TD_mod::set_variates(param.begin() + 2, param.end());
	}
};

template<typename tProjector, typename tParameters >
typename tProjector::tProjection project(
    tProjector& projector,
    const tParameters& parameters) {
  typename tProjector::tProjection result;
  projector.set_parameters(parameters);
  projector.initialize_from_parameters();
  projector.set_start_conditions();
  projector.project_survival();
  projector.get_survival_projection(result);
  return result;
}

#endif //GUTS_TD_H
