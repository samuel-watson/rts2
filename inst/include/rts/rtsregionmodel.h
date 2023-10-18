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

// this can be tidied up as lots of repeated code
// at the moment it works though!
// Problem is needing to specialise the constructor for the different specialisations

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
           int T,
           const rts::RegionData& region_) : region(region_), 
            model(formula_,data_,colnames_,T,grid_data_), 
            re(model,model.covariance.Q(),model.covariance.Q()), matrix(model,re,false,false), 
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
  MatrixXd intersection_infomat();
  sparse grid_to_region_multiplier_matrix();
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
           int T, int m,
           const rts::RegionData& region_,
           const rts::griddata& grid_) : region(region_), 
            model(formula_,data_,colnames_,T, m, grid_, grid_data_), 
            re(model,model.covariance.Q(),model.covariance.Q()), matrix(model,re,false,false), 
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
  MatrixXd intersection_infomat();
  sparse grid_to_region_multiplier_matrix();
};

template<>
class rtsRegionModel<BitsHSGP> {
public:
  rts::RegionData region;
  BitsHSGP model;
  glmmr::RandomEffects<BitsHSGP> re;
  glmmr::ModelMatrix<BitsHSGP> matrix;
  rts::rtsRegionModelOptim<BitsHSGP> optim;
  
  rtsRegionModel(const std::string& formula_,
                 const ArrayXXd& data_,
                 const ArrayXXd& grid_data_,
                 const strvec& colnames_,
                 int T, int m,
                 const ArrayXd& L,
                 const rts::RegionData& region_) : region(region_), 
                 model(formula_,data_,colnames_,T, m, L, grid_data_), 
                 re(model,grid_data_.rows()*T,model.covariance.Q()), matrix(model,re,false,false), 
                 optim(model,matrix,re,region) {};
  
  rtsRegionModel(const rts::rtsRegionModel<BitsHSGP>& mod) : region(mod.region), model(mod.model), re(mod.re), matrix(mod.matrix), optim(mod.optim) {};
  
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
  MatrixXd intersection_infomat();
  sparse grid_to_region_multiplier_matrix();
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
                 int T,
                 rts::RegionData& region_) : region(region_), model(form_region,form_grid,data_region,data_grid,colnames_region,colnames_grid,T,region), 
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
  MatrixXd intersection_infomat();
  sparse grid_to_region_multiplier_matrix();
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
                 rts::RegionData& region_,
                 const rts::griddata& grid_,
                 int T, int m) : region(region_), model(form_region,form_grid,data_region,data_grid,colnames_region,colnames_grid,region,grid_,T,m), 
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
  MatrixXd intersection_infomat();
  sparse grid_to_region_multiplier_matrix();
};

template<>
class rtsRegionModel<BitsHSGPRegion> {
public:
  rts::RegionData region;
  BitsHSGPRegion model;
  glmmr::RandomEffects<BitsHSGPRegion> re;
  glmmr::ModelMatrix<BitsHSGPRegion> matrix;
  rts::rtsRegionModelOptim<BitsHSGPRegion> optim;
  
  rtsRegionModel(const std::string& form_region,
                 const std::string& form_grid,
                 const Eigen::ArrayXXd &data_region,
                 const Eigen::ArrayXXd &data_grid,
                 const strvec& colnames_region,
                 const strvec& colnames_grid,
                 int T, int m,
                 const ArrayXd& L,
                 rts::RegionData& region_) : region(region_), model(form_region,form_grid,data_region,data_grid,colnames_region,colnames_grid,T,m, L,region), 
                 re(model,model.covariance.Q(),model.covariance.Q()), matrix(model,re), 
                 optim(model,matrix,re,region)  {model.linear_predictor.u = &re.u_; };
  
  rtsRegionModel(const rts::rtsRegionModel<BitsHSGPRegion>& mod) : region(mod.region), model(mod.model), re(mod.re), matrix(mod.matrix), optim(mod.optim) {};
  
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
  MatrixXd intersection_infomat();
  sparse grid_to_region_multiplier_matrix();
};

}

// most of this is repeated code
// need to find a way to specialise template and reuse code.

inline MatrixXd rts::rtsRegionModel<BitsAR>::intersection_infomat(){
  MatrixXd A = region.region_design_matrix();
  MatrixXd B = region.grid_design_matrix();
  VectorXd mu = A * model.linear_predictor.X() * model.linear_predictor.parameter_vector();
  VectorXd meanzu = re.zu_.rowwise().mean();
  mu += B * meanzu;
  mu += A * model.data.offset;
  for(int i = 0; i < region.gridT; i++){
    mu.segment(i*region.q_weights.size(),region.q_weights.size()) += region.q_weights.log().matrix();
  }
  mu = -1.0 * mu;
  VectorXd w = mu.array().exp().matrix();
  MatrixXd BDB = B * model.covariance.D(false,false) * B.transpose();
  BDB += w.asDiagonal();
  BDB = BDB.llt().solve(MatrixXd::Identity(BDB.rows(),BDB.cols()));
  MatrixXd AX = A * model.linear_predictor.X();
  MatrixXd M = AX.transpose() * BDB * AX;
  return M;
}

inline MatrixXd rts::rtsRegionModel<BitsNNGP>::intersection_infomat(){
  MatrixXd A = region.region_design_matrix();
  MatrixXd B = region.grid_design_matrix();
  if(A.cols()!=model.linear_predictor.X().rows())Rcpp::stop("A != X");
  VectorXd mu = A * model.linear_predictor.X() * model.linear_predictor.parameter_vector();
  if(B.rows() != mu.size())Rcpp::stop("B != mu");
  if(B.cols() != re.zu_.rows())Rcpp::stop("B != Zu");
  VectorXd meanzu = re.zu_.rowwise().mean();
  mu += B * meanzu;
  if(A.cols() != model.data.offset.size())Rcpp::stop("A != offset");
  mu += A * model.data.offset;
  for(int i = 0; i < region.gridT; i++){
    mu.segment(i*region.q_weights.size(),region.q_weights.size()) += region.q_weights.log().matrix();
  }
  mu = -1.0 * mu;
  VectorXd w = mu.array().exp().matrix();
  if(B.rows() != w.size())Rcpp::stop("B != w");
  if(B.cols() != model.covariance.D(false,false).rows())Rcpp::stop("B != D");
  MatrixXd BDB = B * model.covariance.D(false,false) * B.transpose();
  BDB += w.asDiagonal();
  BDB = BDB.llt().solve(MatrixXd::Identity(BDB.rows(),BDB.cols()));
  MatrixXd AX = A * model.linear_predictor.X();
  if(AX.rows() != BDB.rows())Rcpp::stop("AX != BDB");
  MatrixXd M = AX.transpose() * BDB * AX;
  return M;
}

inline MatrixXd rts::rtsRegionModel<BitsHSGP>::intersection_infomat(){
  MatrixXd A = region.region_design_matrix();
  MatrixXd B = region.grid_design_matrix();
  if(A.cols()!=model.linear_predictor.X().rows())Rcpp::stop("A != X");
  VectorXd mu = A * model.linear_predictor.X() * model.linear_predictor.parameter_vector();
  if(B.rows() != mu.size())Rcpp::stop("B != mu");
  VectorXd meanzu = re.zu_.rowwise().mean();
  if(B.cols() != meanzu.size())Rcpp::stop("B != Zu");
  mu += B * meanzu;
  if(A.cols() != model.data.offset.size())Rcpp::stop("A != offset");
  mu += A * model.data.offset;
  for(int i = 0; i < region.gridT; i++){
    mu.segment(i*region.q_weights.size(),region.q_weights.size()) += region.q_weights.log().matrix();
  }
  mu = -1.0 * mu;
  VectorXd w = mu.array().exp().matrix();
  if(B.rows() != w.size())Rcpp::stop("B != w");
  if(B.cols() != model.covariance.D(false,false).rows())Rcpp::stop("B != D");
  MatrixXd BDB = B * model.covariance.D(false,false) * B.transpose();
  BDB += w.asDiagonal();
  bool BDBsympd = glmmr::Eigen_ext::issympd(BDB);
  if(!BDBsympd){
    BDB = BDB.llt().solve(MatrixXd::Identity(BDB.rows(),BDB.cols()));
    MatrixXd AX = A * model.linear_predictor.X();
    if(AX.rows() != BDB.rows())Rcpp::stop("AX != BDB");
    MatrixXd M = AX.transpose() * BDB * AX;
    return M;
  } else {
    return MatrixXd::Identity(1,1);
  }
  
}

inline MatrixXd rts::rtsRegionModel<BitsARRegion>::intersection_infomat(){
  MatrixXd A = region.region_design_matrix();
  MatrixXd B = region.grid_design_matrix();
  VectorXd mu = A * model.linear_predictor.region_predictor.X() * model.linear_predictor.region_predictor.parameter_vector();
  mu += B * model.linear_predictor.grid_predictor.X() * model.linear_predictor.grid_predictor.parameter_vector();
  VectorXd meanzu = re.zu_.rowwise().mean();
  mu += B * meanzu;
  mu += A * model.data.offset;
  for(int i = 0; i < region.gridT; i++){
    mu.segment(i*region.q_weights.size(),region.q_weights.size()) += region.q_weights.log().matrix();
  }
  mu = -1.0 * mu;
  VectorXd w = mu.array().exp().matrix();
  MatrixXd BDB = B * model.covariance.D(false,false) * B.transpose();
  BDB += w.asDiagonal();
  BDB = BDB.llt().solve(MatrixXd::Identity(BDB.rows(),BDB.cols()));
  MatrixXd Xcombined(region.q_weights.size(),model.linear_predictor.P());
  Xcombined.leftCols(model.linear_predictor.region_predictor.P()) = A * model.linear_predictor.region_predictor.X();
  Xcombined.rightCols(model.linear_predictor.grid_predictor.P()) = B * model.linear_predictor.grid_predictor.X();
  MatrixXd M = Xcombined.transpose() * BDB * Xcombined;
  return M;
}

inline MatrixXd rts::rtsRegionModel<BitsNNGPRegion>::intersection_infomat(){
  MatrixXd A = region.region_design_matrix();
  MatrixXd B = region.grid_design_matrix();
  VectorXd mu = A * model.linear_predictor.region_predictor.X() * model.linear_predictor.region_predictor.parameter_vector();
  mu += B * model.linear_predictor.grid_predictor.X() * model.linear_predictor.grid_predictor.parameter_vector();
  VectorXd meanzu = re.zu_.rowwise().mean();
  mu += B * meanzu;
  mu += A * model.data.offset;
  for(int i = 0; i < region.gridT; i++){
    mu.segment(i*region.q_weights.size(),region.q_weights.size()) += region.q_weights.log().matrix();
  }
  mu = -1.0 * mu;
  VectorXd w = mu.array().exp().matrix();
  MatrixXd BDB = B * model.covariance.D(false,false) * B.transpose();
  BDB += w.asDiagonal();
  BDB = BDB.llt().solve(MatrixXd::Identity(BDB.rows(),BDB.cols()));
  MatrixXd Xcombined(region.q_weights.size(),model.linear_predictor.P());
  Xcombined.leftCols(model.linear_predictor.region_predictor.P()) = A * model.linear_predictor.region_predictor.X();
  Xcombined.rightCols(model.linear_predictor.grid_predictor.P()) = B * model.linear_predictor.grid_predictor.X();
  MatrixXd M = Xcombined.transpose() * BDB * Xcombined;
  return M;
}

inline MatrixXd rts::rtsRegionModel<BitsHSGPRegion>::intersection_infomat(){
  MatrixXd A = region.region_design_matrix();
  MatrixXd B = region.grid_design_matrix();
  VectorXd mu = A * model.linear_predictor.region_predictor.X() * model.linear_predictor.region_predictor.parameter_vector();
  mu += B * model.linear_predictor.grid_predictor.X() * model.linear_predictor.grid_predictor.parameter_vector();
  VectorXd meanzu = re.zu_.rowwise().mean();
  mu += B * meanzu;
  mu += A * model.data.offset;
  for(int i = 0; i < region.gridT; i++){
    mu.segment(i*region.q_weights.size(),region.q_weights.size()) += region.q_weights.log().matrix();
  }
  mu = -1.0 * mu;
  VectorXd w = mu.array().exp().matrix();
  MatrixXd BDB = B * model.covariance.D(false,false) * B.transpose();
  BDB += w.asDiagonal();
  bool BDBsympd = glmmr::Eigen_ext::issympd(BDB);
  if(!BDBsympd){
    BDB = BDB.llt().solve(MatrixXd::Identity(BDB.rows(),BDB.cols()));
    MatrixXd Xcombined(region.q_weights.size(),model.linear_predictor.P());
    Xcombined.leftCols(model.linear_predictor.region_predictor.P()) = A * model.linear_predictor.region_predictor.X();
    Xcombined.rightCols(model.linear_predictor.grid_predictor.P()) = B * model.linear_predictor.grid_predictor.X();
    MatrixXd M = Xcombined.transpose() * BDB * Xcombined;
    return M;
  } else {
    return MatrixXd::Identity(1,1);
  }
}

inline sparse rts::rtsRegionModel<BitsAR>::grid_to_region_multiplier_matrix(){
  sparse A = region.grid_to_region_matrix();
  ArrayXd xb = model.xb();
  xb = xb.exp();
  for(int i = 0; i < A.Ap.size()-1; i++){
    for(int j = A.Ap[i]; j < A.Ap[i+1]; j++){
      A.Ax[j] *= xb(i);
    }
  }
  return A;
}

inline sparse rts::rtsRegionModel<BitsNNGP>::grid_to_region_multiplier_matrix(){
  sparse A = region.grid_to_region_matrix();
  // ArrayXd xb = model.xb();
  // xb = xb.exp();
  // for(int i = 0; i < A.Ap.size()-1; i++){
  //   for(int j = A.Ap[i]; j < A.Ap[i+1]; j++){
  //     A.Ax[j] *= xb(i);
  //   }
  // }
  return A;
}

inline sparse rts::rtsRegionModel<BitsHSGP>::grid_to_region_multiplier_matrix(){
  sparse A = region.grid_to_region_matrix();
  ArrayXd xb = model.xb();
  xb = xb.exp();
  for(int i = 0; i < A.Ap.size()-1; i++){
    for(int j = A.Ap[i]; j < A.Ap[i+1]; j++){
      A.Ax[j] *= xb(i);
    }
  }
  return A;
}

inline sparse rts::rtsRegionModel<BitsARRegion>::grid_to_region_multiplier_matrix(){
  sparse A = region.grid_to_region_matrix();
  ArrayXd xb = model.linear_predictor.region_predictor.xb();
  xb += model.data.offset.array();
  xb = xb.exp();
  ArrayXd xbg = model.linear_predictor.grid_predictor.xb();
  xbg = xbg.exp();
  for(int i = 0; i < A.Ap.size()-1; i++){
    for(int j = A.Ap[i]; j < A.Ap[i+1]; j++){
      A.Ax[j] *= xb(i);
    }
  }
  // for(int i = 0; i < A.rows(); i++){
  //   A.row(i) *= xb(i);
  // }
  // for(int i = 0; i < A.cols(); i++){
  //   A.col(i) *= xbg(i);
  // }
  for(int i = 0; i < A.Ax.size(); i++){
    A.Ax[i] *= xbg(A.Ai[i]);
  }
  return A;
}

inline sparse rts::rtsRegionModel<BitsNNGPRegion>::grid_to_region_multiplier_matrix(){
  sparse A = region.grid_to_region_matrix();
  ArrayXd xb = model.linear_predictor.region_predictor.xb();
  xb += model.data.offset.array();
  xb = xb.exp();
  ArrayXd xbg = model.linear_predictor.grid_predictor.xb();
  xbg = xbg.exp();
  for(int i = 0; i < A.Ap.size()-1; i++){
    for(int j = A.Ap[i]; j < A.Ap[i+1]; j++){
      A.Ax[j] *= xb(i);
    }
  }
  for(int i = 0; i < A.Ax.size(); i++){
    A.Ax[i] *= xbg(A.Ai[i]);
  }
  return A;
}

inline sparse rts::rtsRegionModel<BitsHSGPRegion>::grid_to_region_multiplier_matrix(){
  sparse A = region.grid_to_region_matrix();
  ArrayXd xb = model.linear_predictor.region_predictor.xb();
  xb += model.data.offset.array();
  xb = xb.exp();
  ArrayXd xbg = model.linear_predictor.grid_predictor.xb();
  xbg = xbg.exp();
  for(int i = 0; i < A.Ap.size()-1; i++){
    for(int j = A.Ap[i]; j < A.Ap[i+1]; j++){
      A.Ax[j] *= xb(i);
    }
  }
  for(int i = 0; i < A.Ax.size(); i++){
    A.Ax[i] *= xbg(A.Ai[i]);
  }
  return A;
}