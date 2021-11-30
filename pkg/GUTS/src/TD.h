/**
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch, alexander singer@rifcon.de, oliver.jakoby@rifcon.de, dirk.nickisch@rifcon.de
 * License GPL-2
 * 2017-10-09 
 * updated: 2019-01-29
 * updated: 2021-11-30 
 */

#ifndef TD_H
#define TD_H

#include "samplers.h"
#include "TD_base.h"
#include "TD_proper.h"
#include "TD_SD.h"
#include "TD_IT.h"


typedef TD<imp_lognormal, 'P' > TD_proper_lognormal;
typedef TD<imp_loglogistic, 'P' > TD_proper_loglogistic;
typedef TD<imp_delta, 'P' > TD_proper_delta;

typedef TD<imp_lognormal, 'I' > TD_IT_imp_lognormal;
typedef TD<imp_loglogistic, 'I' > TD_IT_imp_loglogistic;
typedef TD<lognormal, 'I' > TD_IT_lognormal;
typedef TD<loglogistic, 'I' > TD_IT_loglogistic;

typedef TD<double, 'S' > TD_SD;

#endif //TD_H
