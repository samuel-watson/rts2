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
class rtsModel {
public:
  rtsModel();
  ~rtsModel() = default;
};

template<>
class rtsModel<BitsAR> {
  public:
    BitsAR model;
    glmmr::RandomEffects<BitsAR> re;
    glmmr::ModelMatrix<BitsAR> matrix;
    rts::rtsModelOptim<BitsAR> optim;
    
    rtsModel(const std::string& formula_,
             const ArrayXXd& data_,
             const ArrayXXd& grid_data_,
             const strvec& colnames_,
             int T) : model(formula_,data_,colnames_,T,grid_data_), re(model), matrix(model,re,false,false), optim(model,matrix,re) {};
    
    rtsModel(const rts::rtsModel<BitsAR>& mod) : model(mod.model), re(mod.re), matrix(mod.matrix), optim(mod.optim) {};
    
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
    void update_rho(double rho_){
      model.covariance.update_rho(rho_);
      re.zu_ = model.covariance.ZLu(re.u_);
    }
    void update_u(const MatrixXd &u_){
      if(u_.cols()!=re.u(false).cols()){
        re.u_.conservativeResize(model.covariance.Q(),u_.cols());
        re.zu_.conservativeResize(model.covariance.Q(),u_.cols());
      }
      re.u_ = u_;
      re.zu_ = model.covariance.ZLu(re.u_);
    }
    void set_trace(int trace_){
      optim.trace = trace_;
    }
};

template<>
class rtsModel<BitsNNGP>{
public:
  BitsNNGP model;
  glmmr::RandomEffects<BitsNNGP> re;
  glmmr::ModelMatrix<BitsNNGP> matrix;
  rts::rtsModelOptim<BitsNNGP> optim;
  
  rtsModel(const std::string& formula_,
           const ArrayXXd& data_,
           const ArrayXXd& grid_data_,
           const strvec& colnames_,
           int T, int m,
           const rts::griddata& grid_) : model(formula_,data_,colnames_,T, m, grid_,grid_data_), re(model), matrix(model,re,false,false), optim(model,matrix,re) {};
  
  rtsModel(const rts::rtsModel<BitsNNGP>& mod) : model(mod.model), re(mod.re), matrix(mod.matrix), optim(mod.optim) {};
  
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
  void update_rho(double rho_){
    model.covariance.update_rho(rho_);
    re.zu_ = model.covariance.ZLu(re.u_);
  }
  void update_u(const MatrixXd &u_){
    if(u_.cols()!=re.u(false).cols()){
      re.u_.conservativeResize(model.covariance.Q(),u_.cols());
      re.zu_.conservativeResize(model.covariance.Q(),u_.cols());
    }
    re.u_ = u_;
    re.zu_ = model.covariance.ZLu(re.u_);
  }
  void set_trace(int trace_){
    optim.trace = trace_;
  }
};

template<>
class rtsModel<BitsHSGP>{
public:
  BitsHSGP model;
  glmmr::RandomEffects<BitsHSGP> re;
  glmmr::ModelMatrix<BitsHSGP> matrix;
  rts::rtsModelOptim<BitsHSGP> optim;
  
  rtsModel(const std::string& formula_,
           const ArrayXXd& data_,
           const ArrayXXd& grid_data_,
           const strvec& colnames_,
           int T, int m,
           const ArrayXd& L) : model(formula_,data_,colnames_,T, m, L,grid_data_), re(model), matrix(model,re,false,false), optim(model,matrix,re) {};
  
  rtsModel(const rts::rtsModel<BitsHSGP>& mod) : model(mod.model), re(mod.re), matrix(mod.matrix), optim(mod.optim) {};
  
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
  void update_rho(double rho_){
    model.covariance.update_rho(rho_);
    re.zu_ = model.covariance.ZLu(re.u_);
  }
  void update_u(const MatrixXd &u_){
    if(u_.cols()!=re.u(false).cols()){
      re.u_.conservativeResize(model.covariance.Q(),u_.cols());
      re.zu_.conservativeResize(model.covariance.Q(),u_.cols());
    }
    re.u_ = u_;
    re.zu_ = model.covariance.ZLu(re.u_);
  }
  void set_trace(int trace_){
    optim.trace = trace_;
  }
};

}

