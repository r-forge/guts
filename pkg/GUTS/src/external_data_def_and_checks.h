/**
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch, alexander singer@rifcon.de, oliver.jakoby@rifcon.de, dirk.nickisch@rifcon.de
 * License GPL-2
 * 2021-11-30
 */

#ifndef EXTERNAL_DATA_DEF_AND_CHECKS_H_
#define EXTERNAL_DATA_DEF_AND_CHECKS_H_

#include<cstddef>
#include<algorithm>
#include<memory>
#include<cmath>
#include<string>
#include<exception>
#include "helpers.h"

struct num_discretization_time_steps {
  std::size_t M;
  void set_data(
      const std::size_t new_num_discretization_time_steps
  ) {
    if ((std::isnan(new_num_discretization_time_steps)) || (new_num_discretization_time_steps < 2)) {
      throw std::invalid_argument(
          std::string("The number of discretization time steps M must be an integer >= 2.") + 
            std::string(" Value '") + 
            (std::isnan(new_num_discretization_time_steps) ? std::string("NaN' or 'NA") : std::to_string(new_num_discretization_time_steps)) + 
            std::string("' is not allowed")
      );
    }
    set_data_unchecked(new_num_discretization_time_steps);
  }
  void set_data_unchecked(
      const std::size_t new_num_discretization_time_steps
  ) {
    M = new_num_discretization_time_steps;
  }
  virtual double calculate_dtau() const noexcept = 0;
protected: 
  inline double calculate_dtau(const double experiment_duration) const noexcept {
    return experiment_duration / static_cast<double >(M);
  }
};

struct distribution_sample_size {
  std::size_t N;
  void set_data(
      const std::size_t new_distribution_sample_size
  ) {
    if ((std::isnan(new_distribution_sample_size)) || (new_distribution_sample_size < 3)) {
      throw std::invalid_argument(
          std::string("The distribution sample size N must be an integer >= 3.") + 
            std::string(" Value '") + 
            std::string(std::isnan(new_distribution_sample_size) ? std::string("NaN' or 'NA") : std::to_string(new_distribution_sample_size)) + 
            std::string("' is not allowed")
      );
    }
    set_data_unchecked(new_distribution_sample_size);
  }
  void set_data_unchecked(
      const std::size_t new_distribution_sample_size
  ) {
    N = new_distribution_sample_size;
  }
};

struct svr {
  double SVR;
  void set_data(
      const double new_SVR
  ) {
    if ((std::isnan(new_SVR)) || (new_SVR < 2)) {
      throw std::invalid_argument(
          std::string("The surface volume ration must be a positive value.") + 
            std::string(" Value '") + 
            (std::isnan(new_SVR) ? std::string("NaN' or 'NA") : std::to_string(new_SVR) + 
            std::string("' is not allowed")
            )
      );
    }
    set_data_unchecked(new_SVR);
  }
  void set_data_unchecked(
      const double new_SVR
  ) {
    SVR = new_SVR;
  }
};

inline void throw_invalid_argument(const std::string& specifier = "", const std::string& msg = "") {
  throw std::invalid_argument(
      specifier + std::string(specifier != "" ? ": " : "") +
        msg
  );
}

template<typename C >
inline void throw_invalid_argument_if_contains_nan_or_is_below_0(const C& c, const std::string& specifier = "") {
  for (auto e : c) {
    if (std::isnan(e) | (e < 0) ) {
      throw_invalid_argument(
        specifier, std::string("only positive values allowed.\n") +
          std::string("  Found: ") +
          std::string(std::isnan(e) ? std::string("NaN' or 'NA") : std::to_string(e))
      );
    }
  }
}

template<typename C >
inline void throw_invalid_argument_if_contains_less_than_two_elements(const C& c, const std::string& specifier = "") {
  if (c.size() < 2) {
    throw_invalid_argument(
      specifier, std::string("at least two elements needed.\n") +
        std::string("  Number of elements: ") + std::to_string(c.size())
    );
  }
}

template<typename C >
inline void throw_invalid_argument_if_first_not_0(const C& c, const std::string& specifier = "") {
  if (c.at(0) != 0) {
    throw_invalid_argument(
      specifier, std::string("first value must be zero.\n") +
        std::string("  First element: ") + std::to_string(c.at(0))
    );
  }
}


template<typename C >
inline void throw_invalid_argument_if_not_input_vector(const C& c, const std::string& lable) {
  throw_invalid_argument_if_contains_less_than_two_elements(c, lable);
  throw_invalid_argument_if_contains_nan_or_is_below_0(c, lable);
}

template<typename ttimes >
inline void throw_invalid_argument_if_not_time_vector(const ttimes& times, const std::string& lable) {
  throw_invalid_argument_if_not_input_vector(times, lable);
  if (!std::is_sorted(std::begin(times), std::end(times), std::less_equal<double >() )) {
    throw_invalid_argument(
      lable, std::string("times must be in ascending order.\n")
    );
  }
  throw_invalid_argument_if_first_not_0(times, lable);
}

template<typename ttimes, typename tvalues >
inline void throw_invalid_argument_if_not_time_series(const ttimes& times, const tvalues& values, const std::string& lable) {
  throw_invalid_argument_if_not_time_vector(times, lable + std::string(" measurement times"));
  throw_invalid_argument_if_not_input_vector(times, lable + std::string(" measurements"));
  if (times.size() != values.size()) {
    throw_invalid_argument(
      lable,
      std::string("In time series the number of time points must be similar to the number of values. \n") +
        std::string("Provided are: \n") +
        std::string("  ") + std::to_string(times.size()) + std::string(" time steps \n") +
        std::string("  ") + std::to_string(values.size()) + std::string(" values")
    );
  }
}

template<typename tt, typename tC >
struct exposure {
  using values = tC;
  using times = tt;
 std::shared_ptr<times > Ct;
 std::shared_ptr<values > C;
 inline std::size_t Ct_size() const noexcept {return Ct->size();}
 exposure() : 
   Ct{std::make_shared<times >()}, C{std::make_shared<values >()} {}
 void set_data(
     const exposure::times& new_times,
     const exposure::values& new_values
 ) {
   throw_invalid_argument_if_not_time_series(
     new_times, new_values, 
     "Concentration"
   );
   set_data_unchecked(new_times, new_values);
 }
 void set_data_unchecked(
     const exposure::times& new_times,
     const exposure::values& new_values
 ) {
   Ct = std::make_shared<times >(new_times);
   C = std::make_shared<values >(new_values);
 }
};

template<typename tt >
struct survival_projection {
  using times = tt;
   std::shared_ptr<times > yt;
   inline std::size_t yt_size() const noexcept {return yt->size();}
   inline double experiment_duration() const noexcept {return back(*survival_projection<tt >::yt);}
   survival_projection() : 
     yt{std::make_shared<times >()} {}
   virtual ~survival_projection() {}
   void set_data(
       const survival_projection::times& new_times
   ) {
     throw_invalid_argument_if_not_time_vector(
       new_times,
       "Survival projection"
     );
     set_data_unchecked(new_times);
   }
   void set_data_unchecked(
       const survival_projection::times& new_times
   ) {
     yt = std::make_shared<times >(new_times);
   }
};

template<typename tt, typename tS >
struct survival : public survival_projection<tt > {
  using times = typename survival_projection<tt >::times;
  struct values : public tS {using tS::tS;};
  std::shared_ptr<values > y;
  survival() : survival_projection<tt >{}, y{std::make_shared<values >()} {}                           
  void set_data(
      const survival::times& new_times,
      const survival::values& new_values
  ) {
    throw_invalid_argument_if_not_time_series(
      new_times, new_values, 
      "Survival"
    );
    set_data_unchecked(new_times, new_values);
  }
  void set_data_unchecked(
      const survival::times& new_times,
      const survival::values& new_values
  ) {
    survival_projection<tt >::set_data_unchecked(new_times);
    y = std::make_shared<values >(new_values);
  }
};

template<typename tt, typename tC >
inline void throw_invalid_argument_if_survivals_end_later_than_exposures (
    const typename survival_projection<tt >::times& survival_times, 
    const typename exposure<tt, tC >::times& exposure_times, 
    const std::string& lable
) {
  if (back(survival_times) > back(exposure_times)) {
    throw_invalid_argument(
      lable,
      std::string("Exposure must not end earlier than Survival. \n") +
        std::string("Latest time step: \n") +
        std::string("  Survival: ") + std::to_string(back(survival_times)) + std::string("\n") +
        std::string("  Exposure: ") + std::to_string(back(exposure_times))
    );
  }
}


template<typename tt, typename tC >
struct external_data_base : 
  public exposure<tt, tC >, 
  public survival_projection<tt >,
  public svr
{
  using exposure_times = typename exposure<tt, tC >::times;
  using exposures = typename exposure<tt, tC >::values;
  using survival_times = typename survival_projection<tt >::times;
  void set_data (
      const exposure_times& concentration_measurement_times,
      const exposures& concentration_measurements,
      const survival_times& survival_projection_times,
      const double new_SVR
  ) {
    exposure<tt, tC >::set_data(
        concentration_measurement_times,
        concentration_measurements
    );
    survival_projection<tt >::set_data(
        survival_projection_times
    );
    throw_invalid_argument_if_survivals_end_later_than_exposures<tt, tC >(
        survival_projection_times, concentration_measurement_times, 
        std::string("Projection data")
    );
    svr::set_data(new_SVR);
  }
  void set_data_unchecked (
      const exposure_times& concentration_measurement_times,
      const exposures& concentration_measurements,
      const survival_times& survival_projection_times,
      const double new_SVR
  ) {
    exposure<tt, tC >::set_data_unchecked(
        concentration_measurement_times,
        concentration_measurements
    );
    survival_projection<tt >::set_data_unchecked(
        survival_projection_times
    );
    svr::set_data_unchecked(new_SVR);
  }
protected:
  external_data_base() : 
  exposure<tt, tC >{}, survival_projection<tt >{} {}
};


#endif /* EXTERNAL_DATA_DEF_AND_CHECKS_H_ */
