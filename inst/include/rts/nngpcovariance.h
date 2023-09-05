#pragma once

#define _USE_MATH_DEFINES

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
  int m;
  
  nngpCovariance(const str& formula,
                 const ArrayXXd &data,
                 const strvec& colnames,
                 const dblvec& parameters,
                 int T, int m_,
                 const rts::griddata& grid_) : Covariance(formula, data, colnames, parameters),  
                   grid(grid_),
                   A(m_,data.rows()), Dvec(data.rows()), m(m_) {
    isSparse = false;
    grid.genNN(m);
    gen_AD();
  }
  
  nngpCovariance(const str& formula,
                 const ArrayXXd &data,
                 const strvec& colnames,
                 int T, int m_,
                 const rts::griddata& grid_) : Covariance(formula, data, colnames),  
                 grid(grid_),
                 A(m_,data.rows()), Dvec(data.rows()), m(m_) {
    isSparse = false;
    grid.genNN(m);
  }
  
  nngpCovariance(const rts::nngpCovariance& cov) : Covariance(cov.form_, cov.data_, cov.colnames_, cov.parameters_),  
    grid(cov.grid), A(grid.m,grid.N), Dvec(grid.N), m(cov.m) {
      isSparse = false;
      grid.genNN(m);
      gen_AD();
    }
  
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
  vector_matrix submatrix(int i);
  
// private:
//   //MatrixXd S;
//   //VectorXd Sv;
//   double ll1;
};

}

inline MatrixXd rts::nngpCovariance::D(bool chol, bool upper){
  MatrixXd As = rts::inv_ldlt_AD(A,Dvec,grid.NN);
  if(chol){
    if(upper){
      return As.transpose();
    } else {
      return As;
    }
  } else {
    return As * As.transpose();
  }
}

inline MatrixXd rts::nngpCovariance::ZL(){
  MatrixXd ZL = MatrixXd::Zero(grid.N*grid.T,grid.N*grid.T);
  MatrixXd L = D(true,false);
  
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
  double ll1 = 0.0;
  double logdet = log_determinant();
  
  for(int t = 0; t<grid.T; t++){
    double qf = u(t*grid.N)*u(t*grid.N)/Dvec(0);
#pragma omp parallel for reduction (+:qf)
    for(int i = 1; i < grid.N; i++){
      int idxlim = i <= m ? i : m;
      VectorXd usec(idxlim);
      for(int j = 0; j < idxlim; j++) usec(j) = u(grid.NN(j,i));
      double au = u(t*grid.N + i) - (A.col(i).segment(0,idxlim).transpose() * usec)(0);
      qf += au*au/Dvec(i);
    }
    ll1 -= 0.5*qf + 0.5*grid.N*M_PI; 
  }
  ll1 -= 0.5*logdet;
  return ll1;
}

inline double rts::nngpCovariance::log_determinant(){
  double logdet = Dvec.array().log().sum();
  return logdet * grid.T;
}

inline void rts::nngpCovariance::update_rho(const double rho_){
  rho = rho_;
}

inline void rts::nngpCovariance::gen_AD(){
  A.setZero();
  Dvec.setZero();
  // S.setZero();
  // Sv.setZero();
  
  int idxlim;
  double val = Covariance::get_val(0,0,0);
  Dvec(0) = val;
  
#pragma omp parallel for
  for(int i = 1; i < grid.N; i++){
    idxlim = i <= m ? i : m;
    //S.setZero();
    //Sv.setZero();
    MatrixXd S(idxlim,idxlim);
    VectorXd Sv(idxlim);
    for(int j = 0; j<idxlim; j++){
      S(j,j) = val;
    }
    if(idxlim > 1){
      for(int j = 0; j<(idxlim-1); j++){
        for(int k = j+1; k<idxlim; k++){
          S(j,k) = Covariance::get_val(0,grid.NN(j,i),grid.NN(k,i));
          S(k,j) = S(j,k);
        }
      }
    }
    for(int j = 0; j<idxlim; j++){
      Sv(j) = Covariance::get_val(0,i,grid.NN(j,i));
    }
    // A.block(0,i,idxlim,1) = S.block(0,0,idxlim,idxlim).ldlt().solve(Sv.segment(0,idxlim));
    // Dvec(i) = val - (A.col(i).segment(0,idxlim).transpose() * Sv.segment(0,idxlim))(0);
    A.block(0,i,idxlim,1) = S.ldlt().solve(Sv);
    Dvec(i) = val - (A.col(i).segment(0,idxlim).transpose() * Sv)(0);
  }
}

inline vector_matrix rts::nngpCovariance::submatrix(int i){
  // S.setZero();
  // Sv.setZero();
  int idxlim = i <= m ? i : m;
  double val = Covariance::get_val(0,0,0);
  Dvec(0) = val;
  MatrixXd S(idxlim,idxlim);
  VectorXd Sv(idxlim);
  for(int j = 0; j<idxlim; j++){
    S(j,j) = val;
  }
  if(idxlim > 1){
    for(int j = 0; j<(idxlim-1); j++){
      for(int k = j+1; k<idxlim; k++){
        S(j,k) = Covariance::get_val(0,grid.NN(j,i),grid.NN(k,i));
        S(k,j) = S(j,k);
      }
    }
  }
  for(int j = 0; j<idxlim; j++){
    Sv(j) = Covariance::get_val(0,i,grid.NN(j,i));
  }
  vector_matrix result(idxlim);
  // result.vec = Sv.segment(0,idxlim);
  // result.mat = S.block(0,0,idxlim,idxlim);
  result.vec = Sv;
  result.mat = S;
  return result;
}

inline void rts::nngpCovariance::update_parameters(const dblvec& parameters){
  parameters_ = parameters;
  update_parameters_in_calculators();
  gen_AD();
}

inline void rts::nngpCovariance::update_parameters_extern(const dblvec& parameters){
  parameters_ = parameters;
  update_parameters_in_calculators();
  gen_AD();
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
  gen_AD();
};
