/**
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch, alexander singer@rifcon.de, oliver.jakoby@rifcon.de, dirk.nickisch@rifcon.de
 * License GPL-2
 * 2017-10-09 
 * updated: 2019-01-29
 * updated: 2021-11-30 
 */

#ifndef EXTERNAL_DATA_H_
#define EXTERNAL_DATA_H_

#include "external_data_def_and_checks.h"

/** 
 * \brief Define data structure for projections
 * 
 * The object contains the
 *   - times of concentration measurements (Ct)
 *   - measured concentrations C
 *   - times at which to project survival (yt)
 * it further holds projection parameters such as 
 *   - the number of time steps for time discretization of the ODEs (M)
 *   - the number of steps to discretize the threshold distribution (N)
 *   - the surface-volume ratio (SVR)
 * 
 * Important: Template specializations either include or omit discretizations (M,N),
 * which allows implementation of different solvers for specific GUTS-RED flavors.
 *   
 * For details see external_data_def_and_checks.hpp
 * 
 * \tparam tt time vector type
 * \tparam tC concentration vector type
 * \tparam add_time_discretization if true, M is used in the data set
 * \tparam add_distribution_sample_size if true, N is used in the data set
 */
template<
  typename tt, 
  typename tC, 
  bool add_time_discretization, 
  bool add_distribution_sample_size 
  >
struct external_data {
  external_data() = delete;
};

template<typename tt, typename tC >
struct external_data<tt, tC, false, false > :
  public external_data_base<tt, tC >
{
  using parent = external_data_base<tt, tC >;
  using exposure_times = typename parent::exposure_times;
  using exposures = typename parent::exposures;
  using survival_times = typename parent::survival_times;
};

template<typename tt, typename tC >
struct external_data<tt, tC, true, false > :
  public external_data_base<tt, tC >,
  public num_discretization_time_steps
{
  using parent = external_data_base<tt, tC >;
  using exposure_times = typename parent::exposure_times;
  using exposures = typename parent::exposures;
  using survival_times = typename parent::survival_times;
  void set_data (
      const exposure_times& concentration_measurement_times,
      const exposures& concentration_measurements,
      const survival_times& survival_measurement_times,
      const std::size_t new_num_discretization_time_steps,
      const double new_SVR
  ) {
    num_discretization_time_steps::set_data(
      new_num_discretization_time_steps
    );
    external_data_base<tt, tC >::set_data(
        concentration_measurement_times,
        concentration_measurements,
        survival_measurement_times,
        new_SVR
    );
  }
    void set_data_unchecked (
        const exposure_times& concentration_measurement_times,
        const exposures& concentration_measurements,
        const survival_times& survival_measurement_times,
        const std::size_t new_num_discretization_time_steps,
        const double new_SVR
    ) {
      num_discretization_time_steps::set_data_unchecked(
        new_num_discretization_time_steps
      );
      external_data_base<tt, tC >::set_data_unchecked(
          concentration_measurement_times,
          concentration_measurements,
          survival_measurement_times,
          new_SVR
      );
  }
  inline double calculate_dtau() const noexcept final {
    return num_discretization_time_steps::calculate_dtau(
      this->experiment_duration()
    );
  }
};

template<typename tt, typename tC >
struct external_data<tt, tC, false, true > :
  public external_data_base<tt, tC >,
  public distribution_sample_size
{
  using parent = external_data_base<tt, tC >;
  using exposure_times = typename parent::exposure_times;
  using exposures = typename parent::exposures;
  using survival_times = typename parent::survival_times;
  void set_data (
      const exposure_times& concentration_measurement_times,
      const exposures& concentration_measurements,
      const survival_times& survival_measurement_times,
      const std::size_t distribution_sample_size,
      const double new_SVR
  ) {
    distribution_sample_size::set_data(
      distribution_sample_size
    );
    parent::set_data(
      concentration_measurement_times,
      concentration_measurements,
      survival_measurement_times,
      new_SVR
    );
  }
  void set_data_unchecked (
      const exposure_times& concentration_measurement_times,
      const exposures& concentration_measurements,
      const survival_times& survival_measurement_times,
      const std::size_t distribution_sample_size,
      const double new_SVR
  ) {
    distribution_sample_size::set_data_unchecked(
      distribution_sample_size
    );
    parent::set_data_unchecked(
      concentration_measurement_times,
      concentration_measurements,
      survival_measurement_times,
      new_SVR
    );
  }
};

template<typename tt, typename tC >
struct external_data<tt, tC, true, true > :
  public external_data_base<tt, tC >,
  public num_discretization_time_steps,
  public distribution_sample_size
{
  using parent = external_data_base<tt, tC >;
  using exposure_times = typename parent::exposure_times;
  using exposures = typename parent::exposures;
  using survival_times = typename parent::survival_times;
  void set_data (
      const exposure_times& concentration_measurement_times,
      const exposures& concentration_measurements,
      const survival_times& survival_measurement_times,
      const std::size_t num_discretization_time_steps,
      const std::size_t distribution_sample_size,
      const double new_SVR
  ) {
    num_discretization_time_steps::set_data(
      num_discretization_time_steps
    );
    distribution_sample_size::set_data(
      distribution_sample_size
    );
    parent::set_data(
      concentration_measurement_times,
      concentration_measurements,
      survival_measurement_times,
      new_SVR
    );
  }
  void set_data_unchecked (
      const exposure_times& concentration_measurement_times,
      const exposures& concentration_measurements,
      const survival_times& survival_measurement_times,
      const std::size_t num_discretization_time_steps,
      const std::size_t distribution_sample_size,
      const double new_SVR
  ) {
    num_discretization_time_steps::set_data_unchecked(
      num_discretization_time_steps
    );
    distribution_sample_size::set_data_unchecked(
      distribution_sample_size
    );
    parent::set_data_unchecked(
      concentration_measurement_times,
      concentration_measurements,
      survival_measurement_times,
      new_SVR
    );
  }
  inline double calculate_dtau() const noexcept final {
    return num_discretization_time_steps::calculate_dtau(this->experiment_duration());
  }
};

#endif /* EXTERNAL_DATA_H_ */
