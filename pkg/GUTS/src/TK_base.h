/**
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch, alexander singer@rifcon.de, oliver.jakoby@rifcon.de, dirk.nickisch@rifcon.de
 * License GPL-2
 * 2017-10-09 
 * updated: 2019-01-29
 * updated: 2021-11-30 
 * updated: 2022-01-17
 */

#ifndef TK_BASE_H
#define TK_BASE_H

#include <cstddef>
#include <vector>

//#include "external_data.hpp"

/**
 * @class abstract TK interface
 * 
 * @brief Solver for the TK differential equation to calculate the damage
 * @details The method double calculate_current_D() solves the TK equation at the next discretization step and returns the damage.
 * The function is called within a while-loop that iterates until bool is_timestep_in_range(const double) fails.
 * The discretization iterator is automatically updated.
 */
struct TK {
  TK() {}
  virtual ~TK() {}
  virtual double calculate_damage(
      const std::size_t k, const double tau) const = 0;
  virtual void set_start_conditions() const = 0;
  virtual void initialize_from_parameters() = 0; 
protected:
  virtual void update_to_next_concentration_measurement() const = 0;
};
#endif //TK_BASE_H
