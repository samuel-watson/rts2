#pragma once

#include <glmmr/general.h>
#include <glmmr/modelbits.hpp>
#include <glmmr/randomeffects.hpp>
#include <glmmr/modelmatrix.hpp>
#include <glmmr/modelmcmc.hpp>
#include <glmmr/family.hpp>
#include <glmmr/modelextradata.hpp>
#include <glmmr/calculator.hpp>
#include <glmmr/formula.hpp>
#include "ar1covariance.h"
#include "nngpcovariance.h"
#include "regionlinearpredictor.h"
#include "regiondata.h"

namespace rts {

using namespace Eigen;
using namespace glmmr;


template<typename cov, typename linpred>
class rtsModelBits {
  rtsModelBits(){};
  ~rtsModelBits() = default;
};

template<>
class rtsModelBits<rts::ar1Covariance, LinearPredictor> {
public:
  glmmr::Formula formula;
  rts::ar1Covariance covariance;
  LinearPredictor linear_predictor;
  glmmr::ModelExtraData data;
  glmmr::Family family;
  glmmr::calculator calc;
  glmmr::calculator vcalc;
  bool weighted = false;
  
  rtsModelBits(const std::string& formula_,
               const ArrayXXd& data_,
               const strvec& colnames_,
               std::string family_, 
               std::string link_,
               int T) : 
  formula(formula_), 
  covariance(formula_,data_,colnames_, T),
  linear_predictor(formula,data_,colnames_),
  data(data_.rows()),
  family(family_,link_) { setup_calculator(); };
  
  int n(){return linear_predictor.n();};
  ArrayXd xb(){return linear_predictor.xb() + data.offset;};
  void setup_calculator(){
    dblvec yvec(n(),0.0);
    calc = linear_predictor.calc;
    glmmr::linear_predictor_to_link(calc,family.link);
    glmmr::link_to_likelihood(calc,family.family);
    calc.y = yvec;
    calc.variance.conservativeResize(yvec.size());
    calc.variance = data.variance;
    vcalc = linear_predictor.calc;
    glmmr::re_linear_predictor(vcalc,covariance.Q());
    glmmr::linear_predictor_to_link(vcalc,family.link);
    glmmr::link_to_likelihood(vcalc,family.family);
    vcalc.y = yvec;
    vcalc.variance.conservativeResize(yvec.size());
    vcalc.variance = data.variance;
  }
  
};

template<>
class rtsModelBits<rts::nngpCovariance, LinearPredictor> {
public:
  glmmr::Formula formula;
  rts::nngpCovariance covariance;
  LinearPredictor linear_predictor;
  glmmr::ModelExtraData data;
  glmmr::Family family;
  glmmr::calculator calc;
  glmmr::calculator vcalc;
  bool weighted = false;
  
  rtsModelBits(const std::string& formula_,
               const ArrayXXd& data_,
               const strvec& colnames_,
               std::string family_, 
               std::string link_,
               int T, int m) : 
  formula(formula_), 
  covariance(formula_,data_,colnames_, T, m),
  linear_predictor(formula,data_,colnames_),
  data(data_.rows()),
  family(family_,link_) { setup_calculator(); };
  
  int n(){return linear_predictor.n();};
  ArrayXd xb(){return linear_predictor.xb() + data.offset;};
  void setup_calculator(){
    dblvec yvec(n(),0.0);
    calc = linear_predictor.calc;
    glmmr::linear_predictor_to_link(calc,family.link);
    glmmr::link_to_likelihood(calc,family.family);
    calc.y = yvec;
    calc.variance.conservativeResize(yvec.size());
    calc.variance = data.variance;
    vcalc = linear_predictor.calc;
    glmmr::re_linear_predictor(vcalc,covariance.Q());
    glmmr::linear_predictor_to_link(vcalc,family.link);
    glmmr::link_to_likelihood(vcalc,family.family);
    vcalc.y = yvec;
    vcalc.variance.conservativeResize(yvec.size());
    vcalc.variance = data.variance;
  }
};

// need to implement a calculator for the region model
// to enable the gradient functions.

template<>
class rtsModelBits<rts::ar1Covariance, rts::regionLinearPredictor> {
public:
  glmmr::Formula formula_region;
  glmmr::Formula formula;
  rts::ar1Covariance covariance;
  rts::regionLinearPredictor linear_predictor;
  glmmr::ModelExtraData data;
  glmmr::Family family;
  glmmr::calculator calc;
  glmmr::calculator vcalc;
  bool weighted = false;
  
  rtsModelBits(const std::string& form_region,
               const std::string& form_grid,
               const Eigen::ArrayXXd &data_region,
               const Eigen::ArrayXXd &data_grid,
               const strvec& colnames_region,
               const strvec& colnames_grid,
               std::string family_, 
               std::string link_,
               int T,
               rts::RegionData& region) : 
  formula_region(form_region),
  formula(form_grid), 
  covariance(form_grid,data_grid,colnames_grid, T),
  linear_predictor(formula_region,formula,data_region,data_grid,colnames_region,colnames_grid,region),
  data(data_region.rows()),
  family(family_,link_) {};
  
  int n(){return linear_predictor.n();};
  ArrayXd xb(){return linear_predictor.xb() + data.offset;};
  
};

template<>
class rtsModelBits<rts::nngpCovariance, rts::regionLinearPredictor> {
public:
  glmmr::Formula formula_region;
  glmmr::Formula formula;
  rts::nngpCovariance covariance;
  rts::regionLinearPredictor linear_predictor;
  glmmr::ModelExtraData data;
  glmmr::Family family;
  glmmr::calculator calc;
  glmmr::calculator vcalc;
  bool weighted = false;
  
  rtsModelBits(const std::string& form_region,
               const std::string& form_grid,
               const Eigen::ArrayXXd &data_region,
               const Eigen::ArrayXXd &data_grid,
               const strvec& colnames_region,
               const strvec& colnames_grid,
               std::string family_, 
               std::string link_,
               rts::RegionData& region,
               int T, int m) : 
  formula_region(form_region),
  formula(form_grid), 
  covariance(form_grid,data_grid,colnames_grid, T, m),
  linear_predictor(formula_region,formula,data_region,data_grid,colnames_region,colnames_grid,region),
  data(data_region.rows()),
  family(family_,link_) {};
  
  int n(){return linear_predictor.n();};
  ArrayXd xb(){return linear_predictor.xb() + data.offset;};
 
};

}

