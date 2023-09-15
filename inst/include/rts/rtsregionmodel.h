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
  rts::RegionData region;
  BitsAR model;
  glmmr::RandomEffects<BitsAR> re;
  glmmr::ModelMatrix<BitsAR> matrix;
  rts::rtsRegionModelOptim<BitsAR> optim;
  
  rtsRegionModel(const std::string& formula_,
           const ArrayXXd& data_,
           const ArrayXXd& grid_data_,
           const strvec& colnames_,
           std::string family_, 
           std::string link_,
           int T,
           const rts::RegionData& region_) : region(region_), 
            model(formula_,data_,colnames_,family_,link_,T,grid_data_), 
            re(model,model.covariance.Q(),model.covariance.Q()), matrix(model,re), 
            optim(model,matrix,re,region) {};
  
  rtsRegionModel(const rts::rtsRegionModel<BitsAR>& mod) : region(mod.region), model(mod.model), re(mod.re), matrix(mod.matrix), optim(mod.optim) {};
  
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
  rts::RegionData region;
  BitsNNGP model;
  glmmr::RandomEffects<BitsNNGP> re;
  glmmr::ModelMatrix<BitsNNGP> matrix;
  rts::rtsRegionModelOptim<BitsNNGP> optim;
  
  rtsRegionModel(const std::string& formula_,
           const ArrayXXd& data_,
           const ArrayXXd& grid_data_,
           const strvec& colnames_,
           std::string family_, 
           std::string link_,
           int T, int m,
           const rts::RegionData& region_,
           const rts::griddata& grid_) : region(region_), 
            model(formula_,data_,colnames_,family_,link_,T, m, grid_, grid_data_), 
            re(model,model.covariance.Q(),model.covariance.Q()), matrix(model,re), 
            optim(model,matrix,re,region) {};
  
  rtsRegionModel(const rts::rtsRegionModel<BitsNNGP>& mod) : region(mod.region), model(mod.model), re(mod.re), matrix(mod.matrix), optim(mod.optim) {};
  
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
  rts::RegionData region;
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
                 rts::RegionData& region_) : region(region_), model(form_region,form_grid,data_region,data_grid,colnames_region,colnames_grid,family_,link_,T,region), 
                    re(model,model.covariance.Q(),model.covariance.Q()), matrix(model,re), 
                    optim(model,matrix,re,region)  {model.linear_predictor.u = &re.u_; };
  
  rtsRegionModel(const rts::rtsRegionModel<BitsARRegion>& mod) : region(mod.region), model(mod.model), re(mod.re), matrix(mod.matrix), optim(mod.optim) {};
  
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
  rts::RegionData region;
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
                 rts::RegionData& region_,
                 const rts::griddata& grid_,
                 int T, int m) : region(region_), model(form_region,form_grid,data_region,data_grid,colnames_region,colnames_grid,family_,link_,region,grid_,T,m), 
                    re(model,model.covariance.Q(),model.covariance.Q()), matrix(model,re), 
                    optim(model,matrix,re,region)  {model.linear_predictor.u = &re.u_; };
  
  rtsRegionModel(const rts::rtsRegionModel<BitsNNGPRegion>& mod) : region(mod.region), model(mod.model), re(mod.re), matrix(mod.matrix), optim(mod.optim) {};
  
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

