/**
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 * Function guts_engine
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch, alexander singer@rifcon.de, oliver.jakoby@rifcon.de, dirk.nickisch@rifcon.de
 * License GPL-2
 * 2017-10-09 
 * updated: 2019-01-29
 * updated: 2021-11-30 
 * updated: 2022-01-17
 * updated: 2022-02-01
 */

#include <Rcpp.h>
#include <cctype>
#include <iterator>
#include <vector>
#include "GUTS_RED.h"
#include "external_data.h"

typedef Rcpp::NumericVector ttime;
typedef Rcpp::NumericVector tconc;
typedef Rcpp::NumericVector tpara;
typedef std::vector<double > tsurv;
typedef Rcpp::IntegerVector tobssurv;
typedef R_xlen_t vec_size_t; 


enum TD_type {
  PROPER = 0,
  IT = 1,
  SD = 2
};

// RCPP_EXPOSED_ENUM_NODECL(TD_type)

enum dist_type {
  LOGLOGISTIC = 0,
  LOGNORMAL = 1,
  DELTA = 2,
  EXTERNAL = 3
};

// RCPP_EXPOSED_ENUM_NODECL(Dist_type)

template<typename TD_mod >
struct Rcpp_fast_projector : 
    public guts_projector_fastIT<guts_RED<ttime, tconc, TD_mod, tpara >, ttime, tsurv > {
    typedef TK_RED<ttime, tconc > TK_mod;
    typedef guts_model<TK_mod, TD_mod > guts_RED;
    typedef guts_projector_fastIT<guts_RED, ttime, tsurv > parent;
    template<bool add_distribution_sample_size >
    void add_data(
            const external_data<ttime, tconc, false, add_distribution_sample_size >& data) {
        this->initialize(data);
    }
    Rcpp::NumericVector predict(const tpara& parameters) {
        tsurv survival_probabilities = project(*this, parameters);
        return  Rcpp::wrap(survival_probabilities);
    }
    std::vector<double > get_D() const {
    	return this -> get_damage();
    }
    std::vector<double > get_Dt() const {
    	return this -> get_damage_time();
    }
};

template<typename TD_mod >
struct Rcpp_projector : 
    public guts_projector<guts_RED<ttime, tconc, TD_mod, tpara >, ttime, tsurv > {
    typedef TK_RED<ttime, tconc > TK_mod;
    typedef guts_model<TK_mod, TD_mod > guts_RED;
    typedef guts_projector<guts_RED, ttime, tsurv > parent;
    template<bool add_distribution_sample_size >
    void add_data(
            const external_data<ttime, tconc, true, add_distribution_sample_size >& data) {
        this->initialize(data);
    }
    Rcpp::NumericVector predict(const tpara& parameters) {
        tsurv survival_probabilities = project(*this, parameters);
        return  Rcpp::wrap(survival_probabilities);
    }
    std::vector<double > get_D() const {
      return this -> get_damage();
    }
    std::vector<double > get_Dt() const {
      return this -> get_damage_time();
    }
};

typedef external_data<ttime, tconc, true, true > ext_dat_timediscrete_thresholddistdiscrete;
typedef external_data<ttime, tconc, true, false > ext_dat_timediscrete;
typedef external_data<ttime, tconc, false, true > ext_dat_thresholddistdiscrete;
typedef external_data<ttime, tconc, false, false > ext_dat;

// Appends the random distribution of threshold values to the parameter vector
// 
// \code{z_dist} is transformed to \code{Rcpp::NumericVector} and sorted.
// If \code{z_dist == NULL} an error is thrown
//
// @param par GUTS-RED parameters
// @param z_dist unsorted random distribution of threshold values
// 
// @return the combined vector
Rcpp::NumericVector combine_par_and_external_distribution(
    Rcpp::NumericVector par,
    Rcpp::Nullable<Rcpp::NumericVector > z_dist
    ) {
  if (z_dist.isNull()) Rcpp::stop("dist = external: Need threshold sample");
  Rcpp::NumericVector zd(Rcpp::clone(z_dist));
  zd.sort();
  Rcpp::NumericVector all_values(par.size() + zd.size());
  std::copy(par.begin(), par.end(), all_values.begin());
  std::copy(zd.begin(), zd.end(), all_values.begin() + par.size());
  return all_values;
}

template<typename tProjector, typename tData, typename tPara >
void project_to_gobj(Rcpp::List gobj, tProjector& proj, const tData& dat, const tPara& par) {
  proj.add_data(dat);
  gobj["S"] = proj.predict(par);
  gobj["D"] = proj.get_D();
  gobj["Dt"] = proj.get_Dt();
}

// [[Rcpp::export]]
void guts_engine( Rcpp::List gobj, Rcpp::NumericVector par, Rcpp::Nullable<Rcpp::NumericVector > z_dist = R_NilValue) {
  if (!gobj.inherits("GUTS")) {
    Rcpp::stop( "No GUTS object. Use `guts_setup()` to create or modify objects." );
  }
  tpara par_obj = gobj["par"];
  vec_size_t par_len = par_obj.length();
  //unsigned par_len = static_cast<unsigned >(gobj.attr("par_len"));
  switch (static_cast<unsigned >(gobj.attr("TD_type"))) {
  case TD_type::IT : {
    ext_dat dat;
    dat.set_data_unchecked(gobj["Ct"], gobj["C"], gobj["yt"], gobj["SVR"]);
    switch (static_cast<unsigned >(gobj.attr("dist_type"))) {
    case dist_type::LOGLOGISTIC : {
      if (par.size() != par_len) Rcpp::stop("IT-loglogistic: Need parameters hb, kd, mn and beta"); 
      Rcpp_fast_projector<TD_IT_loglogistic > proj;
      project_to_gobj(gobj, proj, dat, Rcpp::NumericVector::create(par[0], par[1], NA_REAL, par[2], par[3]));
      break;
    }
    case dist_type::LOGNORMAL : {
      if (par.size() != par_len) Rcpp::stop("IT-lognormal: Need parameters hb, kd, mn and sd"); 
      Rcpp_fast_projector<TD_IT_lognormal > proj;
      project_to_gobj(gobj, proj, dat, Rcpp::NumericVector::create(par[0], par[1], NA_REAL, par[2], par[3]));
      break;
    }
    case dist_type::EXTERNAL : {
      if (par.size() != par_len) Rcpp::stop("IT-external: Need parameters hb and kd"); 
      Rcpp_fast_projector<TD<random_sample<tpara > , 'I' > > proj;
      project_to_gobj(gobj, proj, dat, 
                      combine_par_and_external_distribution(par, z_dist)
      );
      break;
    }
    default :
      Rcpp::stop("model 'IT' needs one of the distributions 'loglogistic', 'lognormal' or 'external'");
      break;
    }
    break;
  }
  case TD_type::SD : {
    if (par.size() != par_len) Rcpp::stop("SD: Need parameters hb, kd, kk and mn"); 
    ext_dat_timediscrete dat;
    dat.set_data_unchecked(gobj["Ct"], gobj["C"], gobj["yt"], gobj["M"], gobj["SVR"]);
    Rcpp_projector<TD_SD > proj;
    project_to_gobj(gobj, proj, dat, par);
    break;
  }
  case TD_type::PROPER : {
    switch (static_cast<unsigned >(gobj.attr("dist_type"))) {
    case dist_type::LOGLOGISTIC : {
      if (par.size() != par_len) Rcpp::stop("Proper-loglogistic: Need parameters hb, kd, kk, mn and beta"); 
      ext_dat_timediscrete_thresholddistdiscrete dat;
      dat.set_data_unchecked(gobj["Ct"], gobj["C"], gobj["yt"], gobj["M"], gobj["N"], gobj["SVR"]);
      Rcpp_projector<TD_proper_loglogistic > proj;
      project_to_gobj(gobj, proj, dat, par);
      break;
    } 
    case dist_type::LOGNORMAL : {
      if (par.size() != par_len) Rcpp::stop("Proper-lognormal: Need parameters hb, kd, kk, mn and sd"); 
      ext_dat_timediscrete_thresholddistdiscrete dat;
      dat.set_data_unchecked(gobj["Ct"], gobj["C"], gobj["yt"], gobj["M"], gobj["N"], gobj["SVR"]);
      Rcpp_projector<TD_proper_lognormal > proj;
      project_to_gobj(gobj, proj, dat, par);
      break;
    }
    case dist_type::DELTA : {
      if (par.size() != par_len) Rcpp::stop("Proper-delta: Need parameters hb, kd, kk and mn"); 
      ext_dat_timediscrete dat;
      dat.set_data_unchecked(gobj["Ct"], gobj["C"], gobj["yt"], gobj["M"], gobj["SVR"]);
      Rcpp_projector<TD_proper_delta > proj;
      project_to_gobj(gobj, proj, dat, par);
      break;
    } 
    case dist_type::EXTERNAL : {
      if (par.size() != par_len) Rcpp::stop("Proper-external: Need parameters hb, kd and kk"); 
      ext_dat_timediscrete dat;
      dat.set_data_unchecked(gobj["Ct"], gobj["C"], gobj["yt"], gobj["M"], gobj["SVR"]);
      Rcpp_projector<TD<random_sample<tpara >, 'P' > > proj;
      project_to_gobj(
        gobj, proj, dat, combine_par_and_external_distribution(par, z_dist)
      );
      break;
    }
    default :
      Rcpp::stop("model 'Proper' needs one of the distributions 'loglogistic', 'lognormal', 'delta' or 'external'");
      break;
    }
  break;
  }
  default : 
    Rcpp::stop("model needs to be one of 'Proper', 'IT' or 'SD'");
    break;
  }

  gobj["par"] = par;
  gobj["external_dist"] = z_dist;
  gobj["LL"] = calculate_loglikelihood<tsurv, tobssurv >(gobj["S"], gobj["y"]);
  gobj["SPPE"] = calculate_SPPE<tsurv, tobssurv >(gobj["S"], gobj["y"]);
  gobj["squares"] = calculate_sum_of_squares<tsurv, tobssurv >(gobj["S"], gobj["y"]);
}