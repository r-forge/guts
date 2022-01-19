/**
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch, alexander singer@rifcon.de, oliver.jakoby@rifcon.de, dirk.nickisch@rifcon.de
 * License GPL-2
 * 2017-10-09 
 * updated: 2019-01-29
 * updated: 2021-11-30
 * updated: 2022-01-17 
 */

#ifndef TK_single_concentration_H
#define TK_single_concentration_H

#include <memory>
#include <vector>


#include "TK_base.h"


/**
 * @class abstract TK interface
 * 
 * @brief Solver for the TK differential equation to calculate the damage
 * @details The method double calculate_current_D() solves the TK equation at the next discretization step and returns the damage.
 * The function is called within a while-loop that iterates until bool is_timestep_in_range(const double) fails.
 * The discretized iterator is automatically updated.
 */
template<typename tCt, typename tC >
class TK_single_concentration : public TK  {
public:
  template<typename tTDdata >
    inline void initialize(const tTDdata& TDdata) {
  	  initialize(TDdata.Ct, TDdata.C);
    }
  void initialize_from_parameters() override {}
  void set_start_conditions() const override {
	  D = 0;
	  D_k = 0;
  }
  virtual ~TK_single_concentration() {}
protected:
  inline void update_to_next_concentration_measurement() const override {D_k = D;}
	void initialize(
			const std::shared_ptr<const tCt > new_Ct,
			const std::shared_ptr<const tC > new_C
	) {
		Ct = new_Ct;
		C = new_C;
		diffCCt.resize(new_Ct->size()-1);
		differentiateC();
	}
  ///brief time steps of concentration measurements
  std::shared_ptr<const tCt > Ct;
  ///brief concentration measurements
  std::shared_ptr<const tC > C;
  ///brief differential of concentrations C at measurement times Ct
  std::vector<double > diffCCt;
  ///brief current damage
  mutable double D;
  ///brief damage at last concentration measurement time step
  mutable double D_k;
private:
  /**
   * @brief Differentiate the external concentration C
   * 
   * @details External concentration C is numerically differentiated by linearly 
   * interpolating between measurement times Ct.
   * @ diffCCt numeric differential (size = size(Ct) - 1) 
   */
  virtual void differentiateC();
};

template<typename tCt, typename tC > 
void TK_single_concentration<tCt, tC >::differentiateC() {
  for ( auto i = 1; i < Ct->size(); ++i ) {
    diffCCt.at(i-1) = (C->at(i) - C->at(i-1)) /
    		(Ct->at(i) - Ct->at(i-1));
  }
}

#endif //TK_single_concentration_H

