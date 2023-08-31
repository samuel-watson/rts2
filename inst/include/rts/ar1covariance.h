#pragma once

#include <glmmr/covariance.hpp>
#include "griddata.h"

namespace rts {

using namespace Eigen;
using namespace glmmr;

class ar1Covariance : public Covariance {
public:
  double rho = 0.0;
  rts::griddata grid;
  
  ar1Covariance(const str& formula,
                const ArrayXXd &data,
                const strvec& colnames, int T) : Covariance(formula, data, colnames), grid(data, T) { isSparse = false; };
  
  ar1Covariance(const rts::ar1Covariance& cov) : Covariance(cov.form_, cov.data_, cov.colnames_), grid(cov.grid) {isSparse = false; };
  
  MatrixXd ZL() override;
  MatrixXd LZWZL(const VectorXd& w) override;
  MatrixXd ZLu(const MatrixXd& u) override;
  MatrixXd Lu(const MatrixXd& u) override;
  sparse ZL_sparse() override;
  int Q() override;
  double log_likelihood(const VectorXd &u) override;
  double log_determinant() override;
  void update_rho(const double rho_);
  void update_grid(int T);
  void update_parameters(const dblvec& parameters) override;
  void update_parameters(const ArrayXd& parameters) override;
  void update_parameters_extern(const dblvec& parameters) override;
  
};

}

inline void rts::ar1Covariance::update_grid(int T){
  grid.setup(data_,T);
}

inline void rts::ar1Covariance::update_parameters_extern(const dblvec& parameters){
  parameters_ = parameters;
  update_parameters_in_calculators();
};

inline void rts::ar1Covariance::update_parameters(const dblvec& parameters){
  parameters_ = parameters;
  update_parameters_in_calculators();
};

inline void rts::ar1Covariance::update_parameters(const ArrayXd& parameters){
  if(parameters_.size()==0){
    for(unsigned int i = 0; i < parameters.size(); i++){
      parameters_.push_back(parameters(i));
    }
    update_parameters_in_calculators();
  } else {
    for(unsigned int i = 0; i < parameters.size(); i++){
      parameters_[i] = parameters(i);
    }
    update_parameters_in_calculators();
  } 
};

inline MatrixXd rts::ar1Covariance::ZL(){
  MatrixXd ZL = MatrixXd::Zero(grid.N*grid.T,grid.N*grid.T);
  MatrixXd L = Covariance::D(true,false);//glmmr::sparse_to_dense(matL,false);
  
  if(grid.T==1){
    ZL = L;
  } else {
    for(int t=0; t<grid.T;t++){
      for(int s=t; s<grid.T;s++){
        if(t==0){
          ZL.block(s*grid.N,t*grid.N,grid.N,grid.N) = (pow(rho,s)/(1-rho*rho))*L;
        } else {
          ZL.block(s*grid.N,t*grid.N,grid.N,grid.N) = pow(rho,s-t)*L;
        }
      }
    }
  }
  
  return ZL;
}

inline MatrixXd rts::ar1Covariance::LZWZL(const VectorXd& w){
  MatrixXd ZL = rts::ar1Covariance::ZL();
  MatrixXd LZWZL = ZL.transpose() * w.asDiagonal() * ZL;
  LZWZL += MatrixXd::Identity(LZWZL.rows(), LZWZL.cols());
  return LZWZL;
}

inline MatrixXd rts::ar1Covariance::ZLu(const MatrixXd& u){
  MatrixXd ZLu = rts::ar1Covariance::ZL() * u;
  return ZLu;
}

inline MatrixXd rts::ar1Covariance::Lu(const MatrixXd& u){
  MatrixXd LU(u.rows(), u.cols());
  MatrixXd L = Covariance::D(true,false);//glmmr::sparse_to_dense(matL,false);
  for(int t = 0; t<grid.T; t++)LU.block(t*grid.N, 0, grid.N, u.cols()) = L*u;
  return LU;
}

inline sparse rts::ar1Covariance::ZL_sparse(){
  sparse dummy;
  return dummy;
}

inline int rts::ar1Covariance::Q(){
  return grid.N * grid.T;
}

inline double rts::ar1Covariance::log_likelihood(const VectorXd &u){
  double ll = 0;
  if(grid.T==1){
    ll = Covariance::log_likelihood(u);
  } else {
    VectorXd usub = u.head(grid.N);
    ll += Covariance::log_likelihood(usub);
    for(int t=1; t<grid.T; t++){
      usub = u.segment(t*grid.N,grid.N);
      ll += Covariance::log_likelihood(usub);
    }
  }
  return ll;
}

inline double rts::ar1Covariance::log_determinant(){
  double logdet = Covariance::log_determinant();
  return logdet * grid.T;
}

inline void rts::ar1Covariance::update_rho(const double rho_){
  rho = rho_;
}



