/**
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch, alexander singer@rifcon.de, oliver.jakoby@rifcon.de, dirk.nickisch@rifcon.de
 * License GPL-2
 * 2017-10-09 
 * updated: 2019-01-29
 * updated: 2021-11-30 
 * updated: 2022-01-17
 */

#ifndef TD_SD_H
#define TD_SD_H

#include "TD_base.h"
#include "external_data.h"

template<>
class TD<double, 'S' > : public TD_base {
public:
  TD() : TD_base(), E(), dtau(), kk(), kkXdtau(), hb(), z() {}
  virtual ~TD() {}
  template<typename tTDdata >
  inline void initialize(const tTDdata& TDdata) {
	  dtau = TDdata.calculate_dtau();
  }
  void initialize_from_parameters() override {}
  inline void set_start_conditions() const override {E = 0.0;}
  bool is_still_gathering() const override {return true;}
  /**
  * @returns true if there are still survivors
  */
  inline void set_killing_rate(const double new_kk) {
    kkXdtau = new_kk * dtau;
    kk = new_kk;
  }
  inline void set_background_mortality(const double new_hb) {hb = new_hb;}
  inline void set_threshold(const double new_z) {z = new_z;}
  inline double get_killing_rate() const {return kk;}
  inline double get_background_mortality() const {return hb;}
  inline double get_threshold() const {return z;}
  inline void update_to_next_survival_measurement() const override {}
  /**
   *\brief gather an effect from known damage
   * \param[in] D damage
   */
  inline void gather_effect(const double D) const override {
    if ( D > z ) E += z - D;
  }
  /**
   * \returns  calculate survival at time yt
   * \param[in] yt survival measurement time
   */
  inline double calculate_current_survival(const double yt) const override {
    return std::exp(kkXdtau * E - hb * yt);
  }
  
protected:
  ///internally accumulated effect
  mutable double E;
  ///duration of discretization time step
  double dtau;
  ///killing rate
  double kk;
  ///killing rate times discrete time step
  double kkXdtau;
  ///background mortality
  double hb;
  ///threshold value
  double z;
};

#endif //TD_SD_H
