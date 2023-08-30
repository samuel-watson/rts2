#pragma once

#include <glmmr/general.h>
#include <glmmr/modelbits.hpp>
#include <glmmr/randomeffects.hpp>
#include <glmmr/modelmatrix.hpp>
#include <glmmr/family.hpp>
#include <glmmr/modelextradata.hpp>
#include <glmmr/calculator.hpp>
#include <glmmr/formula.hpp>
#include "rtsmodeloptim.h"
#include "rtsmodelbits.h"

namespace rts {

using namespace Eigen;
using namespace glmmr;

template<typename modeltype>
class rtsModelBase {
public:
  glmmr::RandomEffects<modeltype> re;
  glmmr::ModelMatrix<modeltype> matrix;
  rts::rtsModelOptim<modeltype> optim;
  
  rtsModelBase(modeltype& model_) : re(model_), matrix(model_,re), optim(model_,matrix,re) {};
  
  rtsModelBase(const rts::rtsModelBase<modeltype>& mod) : re(mod.re), matrix(mod.matrix), optim(mod.optim) {};
  
};

template<typename modeltype>
class rtsModel : public rtsModelBase<modeltype> {
public:
  rtsRegionModel();
};

template<>
class rtsModel<BitsAR> : public rtsModelBase<BitsAR> {
  public:
    BitsAR model;
    rtsModel(const std::string& formula_,
             const ArrayXXd& data_,
             const strvec& colnames_,
             std::string family_, 
             std::string link_,
             int T) : model(formula_,data_,colnames_,family_,link_,T), rtsModelBase<BitsAR>(model) {};
    void set_offset(const VectorXd& offset_){
      model.data.set_offset(offset_);
    }
    void set_weights(const ArrayXd& weights_){
      model.data.set_weights(weights_);
      if((weights_ != 1.0).any()){
        model.weighted = true;
      }
    }
    void set_y(const VectorXd& y_){
      model.data.update_y(y_);
    }
    void update_beta(const dblvec &beta_){
      model.linear_predictor.update_parameters(beta_);
    }
    void update_theta(const dblvec &theta_){
      model.covariance.update_parameters(theta_);
      re.zu_ = model.covariance.ZLu(re.u_);
    }
    void update_u(const MatrixXd &u_){
      if(u_.cols()!=re.u(false).cols()){
        re.u_.conservativeResize(model.covariance.Q(),u_.cols());
        re.zu_.conservativeResize(model.covariance.Q(),u_.cols());
      }
      re.zu_ = model.covariance.ZLu(re.u_);
    }
    void set_trace(int trace_){
      optim.trace = trace_;
    }
};

template<>
class rtsModel<BitsNNGP> : public rtsModelBase<BitsNNGP>{
public:
  BitsNNGP model;
  rtsModel(const std::string& formula_,
           const ArrayXXd& data_,
           const strvec& colnames_,
           std::string family_, 
           std::string link_,
           int T, int m) : model(formula_,data_,colnames_,family_,link_,T, m), rtsModelBase<BitsNNGP>(model) {};
  void set_offset(const VectorXd& offset_){
    model.data.set_offset(offset_);
  }
  void set_weights(const ArrayXd& weights_){
    model.data.set_weights(weights_);
    if((weights_ != 1.0).any()){
      model.weighted = true;
    }
  }
  void set_y(const VectorXd& y_){
    model.data.update_y(y_);
  }
  void update_beta(const dblvec &beta_){
    model.linear_predictor.update_parameters(beta_);
  }
  void update_theta(const dblvec &theta_){
    model.covariance.update_parameters(theta_);
    re.zu_ = model.covariance.ZLu(re.u_);
  }
  void update_u(const MatrixXd &u_){
    if(u_.cols()!=re.u(false).cols()){
      re.u_.conservativeResize(model.covariance.Q(),u_.cols());
      re.zu_.conservativeResize(model.covariance.Q(),u_.cols());
    }
    re.zu_ = model.covariance.ZLu(re.u_);
  }
  void set_trace(int trace_){
    optim.trace = trace_;
  }
};

}

