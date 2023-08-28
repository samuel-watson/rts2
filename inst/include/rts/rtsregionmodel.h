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
#include "rtsmodelbits.h"
#include "rtsregionmodeloptim.h"

namespace rts {

using namespace Eigen;
using namespace glmmr;

template<typename modeltype>
class rtsRegionModel {
public:
  modeltype& model;
  glmmr::RandomEffects<modeltype> re;
  glmmr::ModelMatrix<modeltype> matrix;
  rts::rtsRegionModelOptim<modeltype> optim;
  
  rtsRegionModel(modeltype& model_) : model(model_), re(model), matrix(model,re), 
    optim(model,matrix,re,model.linear_predictor.region) { model.linear_predictor.u = &re.u_; };
  
  virtual void set_offset(const VectorXd& offset_);
  virtual void set_weights(const ArrayXd& weights_);
  virtual void set_y(const VectorXd& y_);
  virtual void update_beta(const dblvec &beta_);
  virtual void update_theta(const dblvec &theta_);
  virtual void update_u(const MatrixXd &u_);
  virtual void set_trace(int trace_);
  
};

}

template<typename modeltype>
inline void rts::rtsRegionModel<modeltype>::set_offset(const VectorXd& offset_){
  model.data.set_offset(offset_);
}

template<typename modeltype>
inline void rts::rtsRegionModel<modeltype>::set_weights(const ArrayXd& weights_){
  model.data.set_weights(weights_);
  if((weights_ != 1.0).any()){
    model.weighted = true;
  }
}

template<typename modeltype>
inline void rts::rtsRegionModel<modeltype>::set_y(const VectorXd& y_){
  model.data.update_y(y_);
}

template<typename modeltype>
inline void rts::rtsRegionModel<modeltype>::update_beta(const dblvec &beta_){
  model.linear_predictor.update_parameters(beta_);
}

template<typename modeltype>
inline void rts::rtsRegionModel<modeltype>::update_theta(const dblvec &theta_){
  model.covariance.update_parameters(theta_);
  re.zu_ = model.covariance.ZLu(re.u_);
}

template<typename modeltype>
inline void rts::rtsRegionModel<modeltype>::update_u(const MatrixXd &u_){
  if(u_.cols()!=re.u(false).cols()){
    re.u_.conservativeResize(model.covariance.Q(),u_.cols());
    re.zu_.conservativeResize(model.covariance.Q(),u_.cols());
  }
  re.u_ = u_;
  re.zu_ = model.covariance.ZLu(re.u_);
}

template<typename modeltype>
inline void rts::rtsRegionModel<modeltype>::set_trace(int trace_){
  optim.trace = trace_;
}
