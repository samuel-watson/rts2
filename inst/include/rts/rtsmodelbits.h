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

class rtsModelBitsBase {
public:
  glmmr::Formula formula;
  glmmr::ModelExtraData data;
  glmmr::Family family;
  glmmr::calculator calc;
  glmmr::calculator vcalc;
  bool weighted = false;
  
  rtsModelBitsBase(const glmmr::Formula& formula_, 
                 const glmmr::ModelExtraData& data_,
                 const glmmr::Family& family_) : formula(formula_),
    data(data_), family(family_) {};
  
  rtsModelBitsBase(const std::string& formula_, 
                   const ArrayXXd& data_,
                   const std::string& family_,
                   const std::string& link_) : formula(formula_),
                   data(data_.rows()), family(family_,link_) {};
  
  rtsModelBitsBase(const rts::rtsModelBitsBase& bits) : formula(bits.formula),
    data(bits.data), family(bits.family) {};
  
  virtual int n(){ return 0; };
  virtual ArrayXd xb(){return ArrayXd::Zero(1);};
  virtual void setup_calculator(){};
  ~rtsModelBitsBase() = default;
};

template<typename cov, typename linpred>
class rtsModelBits : public rtsModelBitsBase {
  cov covariance;
  linpred linear_predictor;
  rtsModelBits(){};
  ~rtsModelBits() = default;
  int n() override;
  ArrayXd zb() override;
  void setup_calculator() override;
};

template<>
class rtsModelBits<rts::ar1Covariance, LinearPredictor> : public rtsModelBitsBase {
public:
  rts::ar1Covariance covariance;
  LinearPredictor linear_predictor;
  
  rtsModelBits(const std::string& formula_,
               const ArrayXXd& data_,
               const strvec& colnames_,
               std::string family_, 
               std::string link_,
               int T) : 
  rtsModelBitsBase(formula_,data_,family_,link_),
  covariance(formula_,data_,colnames_, T),
  linear_predictor(formula,data_,colnames_) { setup_calculator(); };
  
  rtsModelBits(const rts::rtsModelBits<rts::ar1Covariance, LinearPredictor>& bits) :
    rtsModelBitsBase(bits.formula,bits.data,bits.family),
    covariance(bits.covariance), linear_predictor(bits.linear_predictor) { setup_calculator(); };
  
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
class rtsModelBits<rts::nngpCovariance, LinearPredictor> : public rtsModelBitsBase {
public:
  rts::nngpCovariance covariance;
  LinearPredictor linear_predictor;
  
  rtsModelBits(const std::string& formula_,
               const ArrayXXd& data_,
               const strvec& colnames_,
               std::string family_, 
               std::string link_,
               int T, int m) : 
    rtsModelBitsBase(formula_,data_,family_,link_),
  covariance(formula_,data_,colnames_, T, m),
  linear_predictor(formula,data_,colnames_) { setup_calculator(); };
  
  rtsModelBits(const rts::rtsModelBits<rts::nngpCovariance, LinearPredictor>& bits) : 
    rtsModelBitsBase(bits.formula,bits.data,bits.family),
    covariance(bits.covariance), linear_predictor(bits.linear_predictor) { setup_calculator(); };
  
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
class rtsModelBits<rts::ar1Covariance, rts::regionLinearPredictor> : public rtsModelBitsBase {
public:
  glmmr::Formula formula_grid;
  rts::ar1Covariance covariance;
  rts::regionLinearPredictor linear_predictor;
  
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
  rtsModelBitsBase(form_region,data_region,family_,link_),
  formula_grid(form_grid),
  covariance(form_grid,data_grid,colnames_grid, T),
  linear_predictor(formula_grid,formula,data_region,data_grid,colnames_region,colnames_grid,region) {};
  
  rtsModelBits(const rts::rtsModelBits<rts::ar1Covariance, rts::regionLinearPredictor>& bits) : 
    rtsModelBitsBase(bits.formula,bits.data,bits.family),
    formula_grid(bits.formula_grid),
    covariance(bits.covariance), linear_predictor(bits.linear_predictor) {};
  
  int n(){return linear_predictor.n();};
  ArrayXd xb(){return linear_predictor.xb() + data.offset;};
};

template<>
class rtsModelBits<rts::nngpCovariance, rts::regionLinearPredictor> : public rtsModelBitsBase {
public:
  glmmr::Formula formula_grid;
  rts::nngpCovariance covariance;
  rts::regionLinearPredictor linear_predictor;
  
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
  rtsModelBitsBase(form_region,data_region,family_,link_),
  formula_grid(form_grid),
  covariance(form_grid,data_grid,colnames_grid, T, m),
  linear_predictor(formula_grid,formula,data_region,data_grid,colnames_region,colnames_grid,region) {};
  
  rtsModelBits(const rts::rtsModelBits<rts::nngpCovariance, rts::regionLinearPredictor>& bits) : 
    rtsModelBitsBase(bits.formula,bits.data,bits.family),
    formula_grid(bits.formula_grid),
    covariance(bits.covariance), linear_predictor(bits.linear_predictor) {};
  
  int n(){return linear_predictor.n();};
  ArrayXd xb(){return linear_predictor.xb() + data.offset;};
  
};

}

typedef rts::rtsModelBits<rts::ar1Covariance, glmmr::LinearPredictor> BitsAR;
typedef rts::rtsModelBits<rts::nngpCovariance, glmmr::LinearPredictor> BitsNNGP;
typedef rts::rtsModelBits<rts::ar1Covariance, rts::regionLinearPredictor> BitsARRegion;
typedef rts::rtsModelBits<rts::nngpCovariance, rts::regionLinearPredictor> BitsNNGPRegion;


