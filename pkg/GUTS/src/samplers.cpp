/**
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch, alexander singer@rifcon.de, oliver.jakoby@rifcon.de, dirk.nickisch@rifcon.de
 * License GPL-2
 * 2017-10-09 
 * updated: 2019-01-29
 * updated: 2021-11-30 
 */

#include "samplers.h"

void imp_lognormal::calc_sample() {
  if ( mn == 0.0 && sd != 0 ) {
    throw std::domain_error( "mn = 0 and sd != 0 -- incomplete lognormal model ignored." );
  }
  double sigma2   =  log(   1.0  +  pow( (sd / mn), 2.0 )   );
  double mu       =  log(mn)  -  (0.5 * sigma2);
  double sigmaD   =  sqrt(sigma2) * R;
  
  
  if (sigmaD + mu > 700) {
    throw std::overflow_error( "Approximating lognormal distribution: infinite variates. Please check parameter values." );
  }
  
  std::size_t N = this->z.size();
  double ztmp;
  for ( std::size_t i = 0; i < N; ++i ) {
    ztmp = (2.0 * i - N + 1.0) / (N-1.0);
    this->z[i] = exp( ztmp * sigmaD + mu );
    this->zw[i] = -0.5 * ztmp * ztmp * R * R;
  }
}

void imp_loglogistic::calc_sample() {
  // if scale (wpar3]) <= 0 or shape (wpar[4]) <= 0:
  // the loglogistic distribution is undefined.
  // These cases are excluded.
  if (alpha <= 0) {
    throw std::domain_error( "Loglogistic distribution undefined for scale parameter <= 0. \nPlease check parameter values." );
  }
  if (beta <= 0) {
    throw std::domain_error( "Loglogistic distribution undefined for shape parameter <= 0. \nPlease check parameter values." );
  } else {
    /* if shape (wpar4]) <=1: the loglogistic mode = 0 and mean undefined
    * To avoid loglogistic distributions that peak at 0, wpar[4] <= 1 throws a warning
    * Excluding this distribution shape still allows approximation of a concentration threshold of 0,
    * by setting scale \approx 0
    */
    if (beta <= 1) {
      throw std::domain_error( "Approximating loglogistic distribution: \nShape parameter should be above 1 to avoid an unrealistic concentration threshold distribution that peaks at 0. A concentration threshold close to 0 is better described by a scale parameter that approximates 0. \nNummeric approximation might be wrong. Please check parameter values." );
    }
  }
  
  // parameters are given as alpha = scale and beta = shape
  // transform parameters to mu and s
  double mu  = log(alpha);
  double s   =  1 / beta;
  // if s * R + mu is above 700, z(N-1) -> infty; returning nan for S and LL
  if (s * R + mu > 700) {
    throw std::domain_error( "Approximating loglogistic distribution: infinite variates. \nPlease check parameter values." );
  }
  
  double ztmp;
  double N = this->z.size();
  
  // loglogistic weights formula
  // zw.at(i) =  log(R / 2.0)  - 2.0 *  log( cosh( (log(z.at(i)) - mu) / 2.0 / s ) );
  for ( std::size_t i = 0; i < N; ++i ) {
    ztmp = (2.0 * i - N + 1.0) / (N-1.0);
    this->z[i] = exp( ztmp * s * R + mu );
    this->zw[i] = - 2.0 *  log( cosh( ztmp * R / 2.0 ) );
  }
}


void imp_delta::calc_sample() {
  this->z.assign(this->z.size(), z_val);
  this->zw.assign(this->z.size(), 0.0);
}
