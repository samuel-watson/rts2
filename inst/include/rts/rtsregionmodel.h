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
    rtsRegionModel();
    ~rtsRegionModel() = default;
};

template<>
class rtsRegionModel<BitsAR>{
public:
  BitsAR model;
  glmmr::RandomEffects<BitsAR> re;
  glmmr::ModelMatrix<BitsAR> matrix;
  rts::rtsRegionModelOptim<BitsAR> optim;
  
  rtsRegionModel(const std::string& formula_,
           const ArrayXXd& data_,
           const strvec& colnames_,
           std::string family_, 
           std::string link_,
           int T,
           rts::RegionData& region_) : model(formula_,data_,colnames_,family_,link_,T), re(model), matrix(model,re), 
            optim(model,matrix,re,region_) {};
  
  rtsRegionModel(const rts::rtsRegionModel<BitsAR>& mod) : model(mod.model), re(mod.re), matrix(mod.matrix), optim(mod.optim) {};
  
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
    re.u_ = u_;
    re.zu_ = model.covariance.ZLu(re.u_);
  }
  void set_trace(int trace_){
    optim.trace = trace_;
  }
};

template<>
class rtsRegionModel<BitsNNGP> {
public:
  BitsNNGP model;
  glmmr::RandomEffects<BitsNNGP> re;
  glmmr::ModelMatrix<BitsNNGP> matrix;
  rts::rtsRegionModelOptim<BitsNNGP> optim;
  
  rtsRegionModel(const std::string& formula_,
           const ArrayXXd& data_,
           const strvec& colnames_,
           std::string family_, 
           std::string link_,
           int T, int m,
           rts::RegionData& region_) : model(formula_,data_,colnames_,family_,link_,T, m), re(model), matrix(model,re), 
            optim(model,matrix,re,region_) {};
  
  rtsRegionModel(const rts::rtsRegionModel<BitsNNGP>& mod) : model(mod.model), re(mod.re), matrix(mod.matrix), optim(mod.optim) {};
  
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
    re.u_ = u_;
    re.zu_ = model.covariance.ZLu(re.u_);
  }
  void set_trace(int trace_){
    optim.trace = trace_;
  }
};

template<>
class rtsRegionModel<BitsARRegion> {
public:
  BitsARRegion model;
  glmmr::RandomEffects<BitsARRegion> re;
  glmmr::ModelMatrix<BitsARRegion> matrix;
  rts::rtsRegionModelOptim<BitsARRegion> optim;
  
  rtsRegionModel(const std::string& form_region,
                 const std::string& form_grid,
                 const Eigen::ArrayXXd &data_region,
                 const Eigen::ArrayXXd &data_grid,
                 const strvec& colnames_region,
                 const strvec& colnames_grid,
                 std::string family_, 
                 std::string link_,
                 int T,
                 rts::RegionData& region) : model(form_region,form_grid,data_region,data_grid,colnames_region,colnames_grid,family_,link_,T,region), 
                    re(model), matrix(model,re), 
                    optim(model,matrix,re,region)  {model.linear_predictor.u = &re.u_; };
  
  rtsRegionModel(const rts::rtsRegionModel<BitsARRegion>& mod) : model(mod.model), re(mod.re), matrix(mod.matrix), optim(mod.optim) {};
  
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
    re.u_ = u_;
    re.zu_ = model.covariance.ZLu(re.u_);
  }
  void set_trace(int trace_){
    optim.trace = trace_;
  }
};

template<>
class rtsRegionModel<BitsNNGPRegion> {
public:
  BitsNNGPRegion model;
  glmmr::RandomEffects<BitsNNGPRegion> re;
  glmmr::ModelMatrix<BitsNNGPRegion> matrix;
  rts::rtsRegionModelOptim<BitsNNGPRegion> optim;
  
  rtsRegionModel(const std::string& form_region,
                 const std::string& form_grid,
                 const Eigen::ArrayXXd &data_region,
                 const Eigen::ArrayXXd &data_grid,
                 const strvec& colnames_region,
                 const strvec& colnames_grid,
                 std::string family_, 
                 std::string link_,
                 rts::RegionData& region,
                 int T, int m) : model(form_region,form_grid,data_region,data_grid,colnames_region,colnames_grid,family_,link_,region,T,m), 
                    re(model), matrix(model,re), 
                    optim(model,matrix,re,region)  {model.linear_predictor.u = &re.u_; };
  
  rtsRegionModel(const rts::rtsRegionModel<BitsNNGPRegion>& mod) : model(mod.model), re(mod.re), matrix(mod.matrix), optim(mod.optim) {};
  
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
    re.u_ = u_;
    re.zu_ = model.covariance.ZLu(re.u_);
  }
  void set_trace(int trace_){
    optim.trace = trace_;
  }
};

}

