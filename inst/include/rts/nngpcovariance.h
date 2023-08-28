#pragma once

#include <glmmr/covariance.hpp>
#include "griddata.h"

namespace rts {

using namespace Eigen;
using namespace glmmr;

class nngpCovariance : public Covariance {
public:
  double rho = 0.0;
  rts::griddata grid;
  MatrixXd A;
  VectorXd Dvec;
  
  nngpCovariance(const str& formula,
                 const ArrayXXd &data,
                 const strvec& colnames, 
                 int T, int m) : Covariance(formula, data, colnames),  
                   grid(data, T, m),
                   A(m,grid.N), Dvec(grid.N), 
                   S(m,m), Sv(m), ll1(0.0), ja(1) {
    isSparse = false;
    ja(0) = 0;
    gen_AD();
  }
  
  using Covariance::ZL;
  MatrixXd D(bool chol = true, bool upper = false) override;
  MatrixXd ZL() override;
  MatrixXd LZWZL(const VectorXd& w) override;
  MatrixXd ZLu(const MatrixXd& u) override;
  MatrixXd Lu(const MatrixXd& u) override;
  sparse ZL_sparse() override;
  int Q() override;
  double log_likelihood(const VectorXd &u) override;
  double log_determinant() override;
  void update_rho(const double rho_);
  void gen_AD();
  void update_parameters(const dblvec& parameters) override;
  void update_parameters(const ArrayXd& parameters) override;
  void update_parameters_extern(const dblvec& parameters) override;
  
private:
  MatrixXd S;
  VectorXd Sv;
  double ll1;
  ArrayXi ja;
};

}

inline MatrixXd rts::nngpCovariance::D(bool chol, bool upper){
  if(chol){
    return rts::inv_ldlt_AD(A,Dvec,grid.NN);
  } else {
    Rcpp::stop("Only chol for now.");
  }
}

inline MatrixXd rts::nngpCovariance::ZL(){
  MatrixXd ZL = MatrixXd::Zero(grid.N*grid.T,grid.N*grid.T);
  MatrixXd L = D();
  
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

inline MatrixXd rts::nngpCovariance::LZWZL(const VectorXd& w){
  MatrixXd ZL = rts::nngpCovariance::ZL();
  MatrixXd LZWZL = ZL.transpose() * w.asDiagonal() * ZL;
  LZWZL += MatrixXd::Identity(LZWZL.rows(), LZWZL.cols());
  return LZWZL;
}

inline MatrixXd rts::nngpCovariance::ZLu(const MatrixXd& u){
  MatrixXd ZLu = rts::nngpCovariance::ZL() * u;
  return ZLu;
}

inline MatrixXd rts::nngpCovariance::Lu(const MatrixXd& u){
  MatrixXd LU(u.rows(), u.cols());
  MatrixXd L = D();
  for(int t = 0; t<grid.T; t++)LU.block(t*grid.N, 0, grid.N, u.cols()) = L*u;
  return LU;
}

inline sparse rts::nngpCovariance::ZL_sparse(){
  sparse dummy;
  return dummy;
}

inline int rts::nngpCovariance::Q(){
  return grid.N * grid.T;
}

inline double rts::nngpCovariance::log_likelihood(const VectorXd &u){
  ll1 = 0.0;
  double qf = 0.0;
  double logdet = 0.0;
  
  for(int t = 0; t<grid.T; t++){
    int idxlim;
    for(int i = 0; i < grid.N; i++){
      qf = 0.0;
      logdet = 0.0;
      idxlim = i <= grid.m ? i - 1 : grid.m;
      VectorXd usec(idxlim);
      for(int j = 0; j < idxlim; j++) usec(j) = u(grid.NN(j,i));
      logdet += log(Dvec(i));
      if(i == 0){
        qf += u(t*grid.N)*u(t*grid.N)/Dvec(0); //qf
      } else {
        double au = u(t*grid.N + i) - (A.col(i).segment(0,idxlim).transpose() * usec);
        qf += au*au/Dvec(i);//qf
      }
    }
    ll1 += -0.5*qf - 0.5*logdet - 0.5*grid.N*3.141593; 
  }
  return ll1;
}

inline double rts::nngpCovariance::log_determinant(){
  double logdet = Covariance::log_determinant();
  return logdet * grid.T;
}

inline void rts::nngpCovariance::update_rho(const double rho_){
  rho = rho_;
}

inline void rts::nngpCovariance::gen_AD(){
  A.setZero();
  Dvec.setZero();
  S.setZero();
  Sv.setZero();
  int idxlim;
  double val = get_val(0,0,0);
  Dvec(0) = val;
  
  for(int i = 1; i < grid.N; i++){
    idxlim = i <= grid.m ? i : grid.m;
    S.setZero();
    Sv.setZero();
    for(int j = 0; j<idxlim; j++){
      S(j,j) = val;
    }
    if(idxlim > 1){
      for(int j = 0; j<(idxlim-1); j++){
        for(int k = j+1; k<idxlim; k++){
          S(j,k) = get_val(0,grid.NN(j,i),grid.NN(k,i));
          S(k,j) = S(j,k);
        }
      }
    }
    for(int j = 0; j<idxlim; j++){
      Sv(j) = get_val(0,i,grid.NN(j,i));
    }
    A.block(0,i,idxlim,1) = S.block(0,0,idxlim,idxlim).ldlt().solve(Sv.segment(0,idxlim));
    Dvec(i) = val - (A.col(i).segment(0,idxlim).transpose() * Sv.segment(0,idxlim))(0);
  }
}

inline void rts::nngpCovariance::update_parameters(const dblvec& parameters){
  if(parameters_.size()==0){
    parameters_ = parameters;
    update_parameters_in_calculators();
  } else {
    parameters_ = parameters;
    update_parameters_in_calculators();
  }
}

inline void rts::nngpCovariance::update_parameters_extern(const dblvec& parameters){
  if(parameters_.size()==0){
    parameters_ = parameters;
    update_parameters_in_calculators();
  } else {
    parameters_ = parameters;
    update_parameters_in_calculators();
  }
}

inline void rts::nngpCovariance::update_parameters(const ArrayXd& parameters){
  if(parameters_.size()==0){
    for(unsigned int i = 0; i < parameters.size(); i++){
      parameters_.push_back(parameters(i));
    }
    update_parameters_in_calculators();
  } else if(parameters_.size() == parameters.size()){
    for(unsigned int i = 0; i < parameters.size(); i++){
      parameters_[i] = parameters(i);
    }
    update_parameters_in_calculators();
  } 
};
