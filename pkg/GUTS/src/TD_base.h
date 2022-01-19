/**
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch, alexander singer@rifcon.de, oliver.jakoby@rifcon.de, dirk.nickisch@rifcon.de
 * License GPL-2
 * 2017-10-09 
 * updated: 2019-01-29
 * updated: 2021-11-30
 * updated: 2022-01-17 
 */

#ifndef TD_BASE_H
#define TD_BASE_H

/**
 * @class abstract TD interface
 * 
 * @brief accumulates damage above the threshold and executes respective mortality
 */
class TD_base {
public:
	virtual ~TD_base() {}
  /**
   * @brief gather an effect from known damage
   * @param[in] D damage
   */
  virtual void gather_effect(const double D) const = 0;
  /**
   * @brief calculate survival rate at time yt
   * @param[in] yt time
   * @returns the survival probability
   */
  virtual double calculate_current_survival(const double yt) const = 0;
  /**
   * @brief simulate the number of survivors from the number of survivors in the previous time step
   * @param[in] y_previous: number of survivors in previous time step
   * @results simulated number of survivors
   */
  //virtual std::size_t simulate_current_number_of_survivors(const std::size_t y_previous) = 0;
  /**
   * @returns true if damage has not been gathered for all individuals/threshold values 
   */
  virtual bool is_still_gathering() const = 0;
  virtual void update_to_next_survival_measurement() const = 0;
  virtual void set_start_conditions() const = 0;
  virtual void initialize_from_parameters() = 0;
};

/**
 * \brief abstract template needs specialization
 * \tparam sampler sampler to generate the threshold distribution
 * \tparam TD_type one of three char types
 *    - P proper
 *    - I IT
 *    - S SD
 *  \details this template must be specialized
 */
template< typename sampler, char TD_type >
class TD {TD();};

struct background_mortality {
	background_mortality() :
		hb(std::numeric_limits<double>::quiet_NaN())
	{}
	virtual ~background_mortality() {}
	inline void set_background_mortality(const double new_hb) {hb = new_hb;}
	inline double get_background_mortality() const {return hb;}
protected:
    ///background mortality
    double hb;
};
#endif //TD_BASE_H
