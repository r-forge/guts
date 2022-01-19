/**
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch, alexander singer@rifcon.de, oliver.jakoby@rifcon.de, dirk.nickisch@rifcon.de
 * License GPL-2
 * 2017-10-09 
 * updated: 2019-01-29
 * updated: 2021-11-30 
 * updated: 2022-01-17
 */

#ifndef TK_RED_H
#define TK_RED_H

#include <vector>
#include <cmath>
#include <limits>

#include <iostream>

#include "TK_single_concentration.h"

/**
 * @class TK-RED: assuming the concentration is linearly interpolated between measurements.
 */
template<typename tCt, typename tC >
class TK_RED : public TK_single_concentration<tCt, tC > {
	typedef TK_single_concentration<tCt, tC > parent;
public:
	TK_RED (): parent(), ke(std::numeric_limits<double>::quiet_NaN()) {}
	virtual ~TK_RED () {}
	inline virtual void set_dominant_rate_constant(const double new_ke) {
		ke = new_ke;
		ke_times_SVR = ke * SVR;
	}
	void initialize(
			const std::shared_ptr<const tCt > new_Ct,
			const std::shared_ptr<const tC > new_C,
			const double new_SVR
	) {
		parent::initialize(new_Ct, new_C);
		SVR = new_SVR;
	}
	template<typename tTDdata >
	inline void initialize(const tTDdata& TDdata) {
		initialize(TDdata.Ct, TDdata.C, TDdata.SVR);
	}
	inline double get_dominant_rate_constant() const {return ke;}
	/**
	 * @brief Solve differential damage equation at time $t$
	 *
	 * @details Solves the differential TK equation (e.g. eq. 1) in Albert et al. (2016).
	 * External concentration $C(t)$ is linearly interpolated between measurement time steps $Ct$.
	 * @param[in] t time at which to calculate the damage
	 * @param[in] k index of concentration measurement interval. The index defines the boundary (starting) conditions and must point to the concentration measurement interval in which t lies (i.e. Ct[k] <= t < Ct[k+1])
	 */
	inline double calculate_damage(const std::size_t k, const double t) const override {
		double tmp = exp( -ke_times_SVR * (t - this->Ct->at(k)) );
		double summand3 =
			ke_times_SVR > 0.0  ? (t - this->Ct->at(k) - (1.0-tmp)/ke_times_SVR)  *  this->diffCCt[k] : 0.0;
		this -> D = tmp * (this->D_k - this->C->at(k)) + this->C->at(k) + summand3;
		return this -> D;
	}
	/**
	 * @returns the time $te$ at which the damage assumes an extreme value
	 *
	 * @details solution of the first derivative of damage is 0 ($\frac{dD}{dt} = 0$).
	 * @param[in] k index of concentration measurement interval (points to the beginning of the interval).
	 */
	inline double calculate_time_of_extreme_damage(const std::size_t k) const {
		return log((this->D_k - this->C->at(k))*ke_times_SVR/this->diffCCt.at(k) + 1) / ke_times_SVR + this->Ct->at(k);
	}
	/**
	 * @returns the extreme value of the damage
	 *
	 * @details  Damage at time $te$
	 * @param[in] te time at which the damage assumes an extreme value (as calculated in calculate_time_of_extreme_damage(const std::size_t))
	 * @param[in] k index of concentration measurement interval (points to the beginning of the interval)
	 */
	inline double calculate_extreme_damage(const double te, const std::size_t k) const {
		return this->diffCCt.at(k) * (te - this->Ct->at(k)) + this->C->at(k);
	}
	/**
	 * @returns true if an extreme value at $te$ is a maximum
	 *
	 * @details evaluates if the second derivative of damage at $te$ (as calculated in calculate_time_of_extreme_damage(const std::size_t)) is below or equal 0 ($\frac{d^{2}D}{dt^{2}} \leq 0$). Equal is admitted to ensure that the extreme value is evaluated in case of doubt.
	 * @param[in] k index of concentration measurement interval (points to the beginning of the interval).
	 */
	inline bool is_maximum_damage(const std::size_t k) const {
		return this->D_k < this->Ct->at(k) - this->diffCCt.at(k) / ke_times_SVR;
	}
protected:
	double ke;
	double SVR;
	double ke_times_SVR;
};
#endif //TK_RED_H
