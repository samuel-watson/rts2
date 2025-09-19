#pragma once

#include <glmmr/modelbits.hpp>
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
  rts::RegionData                   region;
  BitsAR                            model;
  glmmr::RandomEffects<BitsAR>      re;
  glmmr::ModelMatrix<BitsAR>        matrix;
  rts::rtsRegionModelOptim<BitsAR>  optim;
  
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
  void update_u(const MatrixXd &u_, bool append = false){
      int newcolsize = u_.cols();
      int currcolsize = re.u_.cols();

      if(append){
        re.u_.conservativeResize(NoChange,currcolsize + newcolsize);
        re.zu_.conservativeResize(NoChange,currcolsize + newcolsize);
        re.u_.rightCols(newcolsize) = u_;
        optim.ll_current.resize(currcolsize + newcolsize,NoChange);
      } else {
        if(u_.cols()!=re.u_.cols()){
          re.u_.resize(NoChange,newcolsize);
          re.zu_.resize(NoChange,newcolsize);
        }
        re.u_ = u_;
        if(newcolsize != optim.ll_current.rows()) optim.ll_current.resize(newcolsize,NoChange);
      }
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
  rts::RegionData                     region;
  BitsNNGP                            model;
  glmmr::RandomEffects<BitsNNGP>      re;
  glmmr::ModelMatrix<BitsNNGP>        matrix;
  rts::rtsRegionModelOptim<BitsNNGP>  optim;
  
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
  void update_u(const MatrixXd &u_, bool append = false){
      int newcolsize = u_.cols();
      int currcolsize = re.u_.cols();

      if(append){
        re.u_.conservativeResize(NoChange,currcolsize + newcolsize);
        re.zu_.conservativeResize(NoChange,currcolsize + newcolsize);
        re.u_.rightCols(newcolsize) = u_;
        optim.ll_current.resize(currcolsize + newcolsize,NoChange);
      } else {
        if(u_.cols()!=re.u_.cols()){
          re.u_.resize(NoChange,newcolsize);
          re.zu_.resize(NoChange,newcolsize);
        }
        re.u_ = u_;
        if(newcolsize != optim.ll_current.rows()) optim.ll_current.resize(newcolsize,NoChange);
      }
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
  rts::RegionData                     region;
  BitsHSGP                            model;
  glmmr::RandomEffects<BitsHSGP>      re;
  glmmr::ModelMatrix<BitsHSGP>        matrix;
  rts::rtsRegionModelOptim<BitsHSGP>  optim;
  
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
  void update_u(const MatrixXd &u_, bool append = false){
      int newcolsize = u_.cols();
      int currcolsize = re.u_.cols();

      if(append){
        re.u_.conservativeResize(NoChange,currcolsize + newcolsize);
        re.zu_.conservativeResize(NoChange,currcolsize + newcolsize);
        re.u_.rightCols(newcolsize) = u_;
        optim.ll_current.resize(currcolsize + newcolsize,NoChange);
      } else {
        if(u_.cols()!=re.u_.cols()){
          re.u_.resize(NoChange,newcolsize);
          re.zu_.resize(NoChange,newcolsize);
        }
        re.u_ = u_;
        if(newcolsize != optim.ll_current.rows()) optim.ll_current.resize(newcolsize,NoChange);
      }
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
  rts::RegionData                         region;
  BitsARRegion                            model;
  glmmr::RandomEffects<BitsARRegion>      re;
  glmmr::ModelMatrix<BitsARRegion>        matrix;
  rts::rtsRegionModelOptim<BitsARRegion>  optim;
  
  rtsRegionModel(const std::string& form_region,
                 const std::string& form_grid,
                 const Eigen::ArrayXXd &data_region,
                 const Eigen::ArrayXXd &data_grid,
                 const strvec& colnames_region,
                 const strvec& colnames_grid,
                 int T,
                 rts::RegionData& region_) : region(region_), model(form_region,form_grid,data_region,data_grid,colnames_region,colnames_grid,T,region), 
                    re(model,model.covariance.Q(),model.covariance.Q()), matrix(model,re), 
                    optim(model,matrix,re,region)  {
    model.linear_predictor.u = &re.zu_; };
  
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
  void update_u(const MatrixXd &u_, bool append = false){
      int newcolsize = u_.cols();
      int currcolsize = re.u_.cols();

      if(append){
        re.u_.conservativeResize(NoChange,currcolsize + newcolsize);
        re.zu_.conservativeResize(NoChange,currcolsize + newcolsize);
        re.u_.rightCols(newcolsize) = u_;
        optim.ll_current.resize(currcolsize + newcolsize,NoChange);
      } else {
        if(u_.cols()!=re.u_.cols()){
          re.u_.resize(NoChange,newcolsize);
          re.zu_.resize(NoChange,newcolsize);
        }
        re.u_ = u_;
        if(newcolsize != optim.ll_current.rows()) optim.ll_current.resize(newcolsize,NoChange);
      }
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
  rts::RegionData                           region;
  BitsNNGPRegion                            model;
  glmmr::RandomEffects<BitsNNGPRegion>      re;
  glmmr::ModelMatrix<BitsNNGPRegion>        matrix;
  rts::rtsRegionModelOptim<BitsNNGPRegion>  optim;
  
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
                    optim(model,matrix,re,region)  {
    model.linear_predictor.u = &re.zu_; };
  
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
  void update_u(const MatrixXd &u_, bool append = false){
      int newcolsize = u_.cols();
      int currcolsize = re.u_.cols();

      if(append){
        re.u_.conservativeResize(NoChange,currcolsize + newcolsize);
        re.zu_.conservativeResize(NoChange,currcolsize + newcolsize);
        re.u_.rightCols(newcolsize) = u_;
        optim.ll_current.resize(currcolsize + newcolsize,NoChange);
      } else {
        if(u_.cols()!=re.u_.cols()){
          re.u_.resize(NoChange,newcolsize);
          re.zu_.resize(NoChange,newcolsize);
        }
        re.u_ = u_;
        if(newcolsize != optim.ll_current.rows()) optim.ll_current.resize(newcolsize,NoChange);
      }
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
  rts::RegionData                             region;
  BitsHSGPRegion                              model;
  glmmr::RandomEffects<BitsHSGPRegion>        re;
  glmmr::ModelMatrix<BitsHSGPRegion>          matrix;
  rts::rtsRegionModelOptim<BitsHSGPRegion>    optim;
  
  rtsRegionModel(const std::string& form_region,
                 const std::string& form_grid,
                 const Eigen::ArrayXXd &data_region,
                 const Eigen::ArrayXXd &data_grid,
                 const strvec& colnames_region,
                 const strvec& colnames_grid,
                 int T, int m,
                 const ArrayXd& L,
                 rts::RegionData& region_) : region(region_), model(form_region,form_grid,data_region,data_grid,colnames_region,colnames_grid,T,m, L,region), 
                 re(model,data_grid.rows()*T,model.covariance.Q()), matrix(model,re), 
                 optim(model,matrix,re,region)  {
    model.linear_predictor.u = &re.zu_; };
  
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
  void update_u(const MatrixXd &u_, bool append = false){
      int newcolsize = u_.cols();
      int currcolsize = re.u_.cols();

      if(append){
        re.u_.conservativeResize(NoChange,currcolsize + newcolsize);
        re.zu_.conservativeResize(NoChange,currcolsize + newcolsize);
        re.u_.rightCols(newcolsize) = u_;
        optim.ll_current.resize(currcolsize + newcolsize,NoChange);
      } else {
        if(u_.cols()!=re.u_.cols()){
          re.u_.resize(NoChange,newcolsize);
          re.zu_.resize(NoChange,newcolsize);
        }
        re.u_ = u_;
        if(newcolsize != optim.ll_current.rows()) optim.ll_current.resize(newcolsize,NoChange);
      }
      re.zu_ = model.covariance.ZLu(re.u_);
    }
  void set_trace(int trace_){
    optim.trace = trace_;
  }
  MatrixXd intersection_infomat();
  sparse grid_to_region_multiplier_matrix();
};

}


inline MatrixXd rts::rtsRegionModel<BitsAR>::intersection_infomat(){
  MatrixXd X = model.linear_predictor.X();
  VectorXd xb = model.linear_predictor.xb();
  xb += model.data.offset;
  xb = xb.array().exp().inverse().matrix();
  rts::ar1Covariance newcov(model.covariance.form_.formula_, model.covariance.grid.X, model.covariance.colnames_, 1);
  newcov.update_parameters(model.covariance.parameters_);
  MatrixXd D = newcov.D(false,false);
  MatrixXd AR = model.covariance.ar_matrix();
  MatrixXd ARD = rts::kronecker(AR,D); // on grid
  MatrixXd Sigma = region.grid_to_region(ARD,false);//sparse_matrix_mult(C,ARD,false);
  MatrixXd Sigma2 = region.grid_to_region(Sigma.transpose(),false);
  Sigma2 += xb.asDiagonal();
  Sigma2 = Sigma2.llt().solve(MatrixXd::Identity(Sigma2.rows(),Sigma2.cols()));
  if(X.rows() != Sigma2.cols()){
    Rcpp::Rcout << "\nX dim: " << X.rows() << " " << X.cols() << " Sigma2 dim " << Sigma2.rows() << " " << Sigma2.cols();
    throw std::runtime_error("X rows != Sigma cols ");
  }
  MatrixXd M = X.transpose() * Sigma2 * X;
  return M;
}

inline MatrixXd rts::rtsRegionModel<BitsNNGP>::intersection_infomat(){
  MatrixXd X = model.linear_predictor.X();
  VectorXd xb = model.linear_predictor.xb();
  xb += model.data.offset;
  xb = xb.array().exp().inverse().matrix();
  rts::ar1Covariance newcov(model.covariance.form_.formula_, model.covariance.grid.X, model.covariance.colnames_, 1);
  newcov.update_parameters(model.covariance.parameters_);
  MatrixXd D = newcov.D(false,false);
  MatrixXd AR = model.covariance.ar_matrix();
  MatrixXd ARD = rts::kronecker(AR,D); // on grid
  MatrixXd Sigma = region.grid_to_region(ARD,false);//sparse_matrix_mult(C,ARD,false);
  MatrixXd Sigma2 = region.grid_to_region(Sigma.transpose(),false);
  Sigma2 += xb.asDiagonal();
  Sigma2 = Sigma2.llt().solve(MatrixXd::Identity(Sigma2.rows(),Sigma2.cols()));
  if(X.rows() != Sigma2.cols()){
    Rcpp::Rcout << "\nX dim: " << X.rows() << " " << X.cols() << " Sigma2 dim " << Sigma2.rows() << " " << Sigma2.cols();
    throw std::runtime_error("X rows != Sigma cols ");
  }
  MatrixXd M = X.transpose() * Sigma2 * X;
  return M;
}

inline MatrixXd rts::rtsRegionModel<BitsHSGP>::intersection_infomat(){
  MatrixXd X = model.linear_predictor.X();
  VectorXd xb = model.linear_predictor.xb();
  xb += model.data.offset;
  xb = xb.array().exp().inverse().matrix();
  rts::ar1Covariance newcov(model.covariance.form_.formula_, model.covariance.grid.X, model.covariance.colnames_, 1);
  newcov.update_parameters(model.covariance.parameters_);
  MatrixXd D = newcov.D(false,false);
  MatrixXd AR = model.covariance.ar_matrix();
  MatrixXd ARD = rts::kronecker(AR,D); // on grid
  MatrixXd Sigma = region.grid_to_region(ARD,false);//sparse_matrix_mult(C,ARD,false);
  MatrixXd Sigma2 = region.grid_to_region(Sigma.transpose(),false);
  Sigma2 += xb.asDiagonal();
  Sigma2 = Sigma2.llt().solve(MatrixXd::Identity(Sigma2.rows(),Sigma2.cols()));
  if(X.rows() != Sigma2.cols()){
    Rcpp::Rcout << "\nX dim: " << X.rows() << " " << X.cols() << " Sigma2 dim " << Sigma2.rows() << " " << Sigma2.cols();
    throw std::runtime_error("X rows != Sigma cols ");
  }
  MatrixXd M = X.transpose() * Sigma2 * X;

    // sparse A = region.region_design_matrix();
  // sparse B = region.grid_design_matrix();
  // int N = model.covariance.grid.N;
  // int T = model.covariance.grid.T;
  // int Q = region.q_weights.size();
  // int R = region.nRegion;
  // sparse Bt = B;
  // Bt.transpose();
  // MatrixXd X = model.linear_predictor.X();
  // MatrixXd AX = sparse_matrix_t_mult(A,X,T); 
  // VectorXd mu = AX * model.linear_predictor.parameter_vector();
  // VectorXd meanzu = re.zu_.rowwise().mean();
  // mu += sparse_vector_t_mult(B,meanzu,T);
  // mu += sparse_vector_t_mult(A,model.data.offset,T);
  // for(int t = 0; t < T; t++){
  //   mu.segment(t*Q,Q) += region.q_weights.log().matrix();
  // }
  // mu = mu.array().exp().matrix();
  // MatrixXd D = model.covariance.D(false,false);
  // D = D.llt().solve(MatrixXd::Identity(D.rows(),D.cols()));
  // MatrixXd AR = model.covariance.ar_matrix();
  // AR = AR.llt().solve(MatrixXd::Identity(AR.rows(),AR.cols()));
  // MatrixXd ARD = rts::kronecker(AR,D);
  // VectorXd WD = sparse_vector_t_mult(Bt,mu,T);
  // ARD += WD.asDiagonal();
  // ARD = ARD.llt().solve(MatrixXd::Identity(ARD.rows(),ARD.cols()));
  // MatrixXd S(AX.rows(),AX.rows());
  // for(int t = 0; t < T; t++){
  //   for(int s = t; s < T; s++){
  //     if(t==s){
  //       MatrixXd Ssub = sparse_matrix_mult(B,ARD.block(s*N,t*N,N,N));
  //       MatrixXd Ssub2 = sparse_matrix_mult(B, Ssub,true);
  //       S.block(s*Q,t*Q,Q,Q) = mu.segment(t*Q,Q).asDiagonal();
  //       S.block(s*Q,t*Q,Q,Q) -= mu.segment(s*Q,Q).asDiagonal()*Ssub2*mu.segment(t*Q,Q).asDiagonal();
  //     } else {
  //       MatrixXd Ssub = sparse_matrix_mult(B,ARD.block(s*N,t*N,N,N));
  //       MatrixXd Ssub2 = sparse_matrix_mult(B, Ssub, true);
  //       S.block(s*Q,t*Q,Q,Q) = -1.0 * mu.segment(s*Q,Q).asDiagonal()*Ssub2*mu.segment(t*Q,Q).asDiagonal();
  //       S.block(t*Q,s*Q,Q,Q) = S.block(s*Q,t*Q,Q,Q).transpose();
  //     }
  //   }
  // }
  // MatrixXd M = AX.transpose() * S * AX;
  return M;
}

inline MatrixXd rts::rtsRegionModel<BitsARRegion>::intersection_infomat(){
  MatrixXd X = model.linear_predictor.X();
  MatrixXd Z = model.linear_predictor.Z();
  VectorXd xb = model.linear_predictor.xb();
  xb += model.data.offset;
  xb = xb.array().exp().inverse().matrix();
  rts::ar1Covariance newcov(model.covariance.form_.formula_, model.covariance.data_, model.covariance.colnames_, 1);
  newcov.update_parameters(model.covariance.parameters_);
  MatrixXd D = newcov.D(false,false);
  MatrixXd AR = model.covariance.ar_matrix();
  MatrixXd ARD = rts::kronecker(AR,D); // on grid
  // MatrixXd Sigma = region.grid_to_region(ARD,false);//sparse_matrix_mult(C,ARD,false);
  // MatrixXd Sigma2 = region.grid_to_region(Sigma.transpose(),false);
  MatrixXd Sigma2 = Z * ARD * Z.transpose();
  Sigma2 += xb.asDiagonal();
  Sigma2 = Sigma2.llt().solve(MatrixXd::Identity(Sigma2.rows(),Sigma2.cols()));
  if(X.rows() != Sigma2.cols()){
    Rcpp::Rcout << "\nX dim: " << X.rows() << " " << X.cols() << " Sigma2 dim " << Sigma2.rows() << " " << Sigma2.cols();
    throw std::runtime_error("X rows != Sigma cols ");
  }
  MatrixXd M = X.transpose() * Sigma2 * X;
  return M;
}

inline MatrixXd rts::rtsRegionModel<BitsNNGPRegion>::intersection_infomat(){
  MatrixXd X = model.linear_predictor.X();
  MatrixXd Z = model.linear_predictor.Z();
  VectorXd xb = model.linear_predictor.xb();
  xb += model.data.offset;
  xb = xb.array().exp().inverse().matrix();
  rts::ar1Covariance newcov(model.covariance.form_.formula_, model.covariance.data_, model.covariance.colnames_, 1);
  newcov.update_parameters(model.covariance.parameters_);
  MatrixXd D = newcov.D(false,false);
  MatrixXd AR = model.covariance.ar_matrix();
  MatrixXd ARD = rts::kronecker(AR,D); // on grid
  // MatrixXd Sigma = region.grid_to_region(ARD,false);//sparse_matrix_mult(C,ARD,false);
  // MatrixXd Sigma2 = region.grid_to_region(Sigma.transpose(),false);
  MatrixXd Sigma2 = Z * ARD * Z.transpose();
  Sigma2 += xb.asDiagonal();
  Sigma2 = Sigma2.llt().solve(MatrixXd::Identity(Sigma2.rows(),Sigma2.cols()));
  if(X.rows() != Sigma2.cols()){
    Rcpp::Rcout << "\nX dim: " << X.rows() << " " << X.cols() << " Sigma2 dim " << Sigma2.rows() << " " << Sigma2.cols();
    throw std::runtime_error("X rows != Sigma cols ");
  }
  MatrixXd M = X.transpose() * Sigma2 * X;
  return M;
}

inline MatrixXd rts::rtsRegionModel<BitsHSGPRegion>::intersection_infomat(){
  VectorXd xbg = model.linear_predictor.grid_predictor.xb();
  VectorXd xbr = model.linear_predictor.region_predictor.xb();
  VectorXd Xg = model.linear_predictor.grid_predictor.X();
  VectorXd Xr = model.linear_predictor.region_predictor.X();
  xbg = xbg.array().exp().matrix();
  xbr = xbr.array().exp().matrix();
  xbr = (xbr.array() * model.data.offset.array()).matrix();
  sparse A = model.linear_predictor.region.grid_to_region_matrix();
  sparse A_mug = sparse_times_diagonal_l(A,xbg);
  sparse A_mug_t(A_mug);
  A_mug_t.transpose();
  VectorXd A_mug_vec = A * xbg;
  ArrayXd ydiva = model.data.y.array() * A_mug_vec.array().inverse();
  sparse A_mug_y = sparse_times_diagonal_l(A_mug_t,ydiva.matrix());
  sparse A_mug_rg = sparse_times_diagonal_l(A,xbr);
  MatrixXd Ag = sparse_to_dense(A_mug_y * A_mug);
  VectorXd h(xbg.size());
  h.setZero();
  VectorXd htmp(h);
  double doth = 0;
  for(int i = 0; i < xbr.size(); i++){
    htmp = sparse_row_hademard_col(A,xbg,i);
    doth = sparse_row_dot_col(A,xbg,i);
    h += htmp * (model.data.y(i)/doth - xbr(i));
  }
  // try it with the existing covariance matrix
  sparse hdiag = make_sparse_diagonal(h);
  A_mug_y *= A_mug;
  //negate
  for(int i = 0; i < A_mug_y.Ax.size(); i++)A_mug_y.Ax[i] *= -1.0;
  hdiag += A_mug_y;
  
  //MatrixXd Wg = h.asDiagonal() - sparse_to_dense(A_mug_y * A_mug);
  MatrixXd D = model.covariance.D();
  D.noalias() = D + sparse_to_dense(hdiag);
  D = D.llt().solve(MatrixXd::Identity(D.rows(), D.cols()));
  VectorXd mu = (xbr.array() * A_mug_vec.array()).matrix();
  sparse mudiag = make_sparse_diagonal(mu);
  sparse wrg = sparse_times_diagonal_l(A_mug_t,xbr);
  wrg.transpose();
  
  MatrixXd WrgD = wrg * D;
  MatrixXd WgD = hdiag * D;
  
  MatrixXd M(Xr.cols() + Xg.cols(), Xr.cols() + Xg.cols());
  MatrixXd Wr = sparse_to_dense(mudiag);
  MatrixXd Wg = sparse_to_dense(hdiag);
  MatrixXd Wrg = sparse_to_dense(wrg);
  wrg.transpose();
  M.block(0,0,Xr.cols(),Xr.cols()) = Xr.transpose() * (Wr - WrgD*wrg) * Xr;
  M.block(Xr.cols(),Xr.cols(),Xg.cols(),Xg.cols()) = Xg.transpose() * (Wg - WgD*hdiag) * Xg;
  M.block(0,Xr.cols(),Xr.cols(),Xg.cols()) = Xr.transpose() * (Wrg - WrgD*hdiag) * Xg;
  M.block(Xr.cols(),0,Xg.cols(),Xr.cols()) = M.block(0,Xr.cols(),Xr.cols(),Xg.cols()).transpose();
  
  return M;
}

inline sparse rts::rtsRegionModel<BitsAR>::grid_to_region_multiplier_matrix(){
  sparse A = region.grid_to_region_matrix();
  return A;
}

inline sparse rts::rtsRegionModel<BitsNNGP>::grid_to_region_multiplier_matrix(){
  sparse A = region.grid_to_region_matrix();
  return A;
}

inline sparse rts::rtsRegionModel<BitsHSGP>::grid_to_region_multiplier_matrix(){
  sparse A = region.grid_to_region_matrix();
  return A;
}

inline sparse rts::rtsRegionModel<BitsARRegion>::grid_to_region_multiplier_matrix(){
  sparse A = region.grid_to_region_matrix();
  return A;
}

inline sparse rts::rtsRegionModel<BitsNNGPRegion>::grid_to_region_multiplier_matrix(){
  sparse A = region.grid_to_region_matrix();
  return A;
}

inline sparse rts::rtsRegionModel<BitsHSGPRegion>::grid_to_region_multiplier_matrix(){
  sparse A = region.grid_to_region_matrix();
  return A;
}