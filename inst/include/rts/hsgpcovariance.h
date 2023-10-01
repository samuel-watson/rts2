#pragma once

#include <glmmr/covariance.hpp>
#include "griddata.h"
#include "rtsmaths.h"

namespace rts {

using namespace Eigen;
using namespace glmmr;

class hsgpCovariance : public Covariance {
public:
  double rho;
  rts::griddata grid;
  int m;
  ArrayXd L_boundary;

  hsgpCovariance(const str& formula,
                const ArrayXXd &data,
                const strvec& colnames,
                int T,
                int m_,
                const ArrayXd& L_boundary_) : Covariance(formula, data, colnames),
                grid(data, T), m(m_), L_boundary(L_boundary_),
                L(m*m,m*m), Lambda(m*m),
                ar_factor(T,T), ar_factor_chol(T,T),
                indices(m*m,2), Phi(grid.N,m*m), PhiT(m*m,m*m) {
    gen_indices();
    gen_phi_prod();
    isSparse = false;
    update_rho(0.1);
  };

  hsgpCovariance(const rts::hsgpCovariance& cov) : Covariance(cov.form_, cov.data_, cov.colnames_), grid(cov.grid),
  L_boundary(cov.L_boundary), L(cov.L), Lambda(cov.Lambda), ar_factor(cov.ar_factor), ar_factor_chol(cov.ar_factor_chol),
  indices(cov.indices), Phi(cov.Phi), PhiT(cov.PhiT) {
    isSparse = false;
    update_rho(cov.rho);
  };

  Vector2d lambda_nD(int i);
  double spd_nD(Vector2d w);
  ArrayXd phi_nD(int i);
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
  void set_function(bool squared_exp);
  MatrixXd PhiSPD();

protected:
  MatrixXd L; // cholesky decomposition of Lambda + PhiTPhi m^2 * m^2
  ArrayXd Lambda;
  MatrixXd ar_factor;
  MatrixXd ar_factor_chol;
  ArrayXXi indices;
  MatrixXd Phi;
  MatrixXd PhiT;
  bool sq_exp = false;
  void gen_indices();
  void gen_phi_prod();
  void Z_chol();
  void update_lambda();
};

}

inline Vector2d rts::hsgpCovariance::lambda_nD(int i){
  Vector2d lambda;
  lambda(0) = (indices(i,0)*M_PI)/(2*L_boundary(0));
  lambda(1) = (indices(i,1)*M_PI)/(2*L_boundary(1));
  lambda(0) *= lambda(0);
  lambda(1) *= lambda(1);
  return lambda;
}

inline double rts::hsgpCovariance::spd_nD(Vector2d w){
  double S;
  Vector2d phisq;
  phisq(0) = this->parameters_[1] * this->parameters_[1];
  phisq(1) = phisq(0);
  Vector2d wsq = (w.array() * w.array()).matrix();
  if(sq_exp){
    S = this->parameters_[0] * 2 * M_PI * phisq(0) * exp(-0.5 * (2*phisq(0)) * (w(0) * w(0) + w(1) * w(1)));
  } else {
    double S1 = this->parameters_[0] * 4 * tgamma(1.5)/tgamma(0.5);
    S = S1 * phisq(0) * 1 /((1 + phisq.transpose() * wsq)*(1 + phisq.transpose() * wsq));
  }
  return S;
}

inline ArrayXd rts::hsgpCovariance::phi_nD(int i){
  ArrayXd fi1(grid.N);
  ArrayXd fi2(grid.N);
  fi1 = (1/sqrt(L_boundary(0))) * sin(indices(i,0)*M_PI*(grid.X.col(0)+L_boundary(0))/(2*L_boundary(0)));
  fi2 = (1/sqrt(L_boundary(1))) * sin(indices(i,1)*M_PI*(grid.X.col(1)+L_boundary(1))/(2*L_boundary(1)));
  fi1 *= fi2;
  return fi1;
}

inline void rts::hsgpCovariance::update_grid(int T){
  grid.setup(data_,T);
}

inline void rts::hsgpCovariance::update_parameters_extern(const dblvec& parameters){
  parameters_ = parameters;
  update_lambda();
};

inline void rts::hsgpCovariance::update_parameters(const dblvec& parameters){
  parameters_ = parameters;
  update_lambda();
};

inline void rts::hsgpCovariance::update_parameters(const ArrayXd& parameters){
  if(parameters_.size()==0){
    for(unsigned int i = 0; i < parameters.size(); i++){
      parameters_.push_back(parameters(i));
    }
    //update_parameters_in_calculators();
  } else {
    for(unsigned int i = 0; i < parameters.size(); i++){
      parameters_[i] = parameters(i);
    }
    //update_parameters_in_calculators();
  }
  update_lambda();
};

inline MatrixXd rts::hsgpCovariance::ZL(){
  MatrixXd ZL = rts::kronecker(ar_factor_chol, PhiSPD());
  return ZL;
}

inline MatrixXd rts::hsgpCovariance::LZWZL(const VectorXd& w){
  MatrixXd ZL = rts::hsgpCovariance::ZL();
  MatrixXd LZWZL = ZL.transpose() * w.asDiagonal() * ZL;
  LZWZL += MatrixXd::Identity(LZWZL.rows(), LZWZL.cols());
  return LZWZL;
}

inline MatrixXd rts::hsgpCovariance::ZLu(const MatrixXd& u){
  MatrixXd ZLu = rts::hsgpCovariance::ZL() * u;
  return ZLu;
}

inline MatrixXd rts::hsgpCovariance::Lu(const MatrixXd& u){
  MatrixXd ZLu = rts::hsgpCovariance::ZL() * u;
  return ZLu;
}

inline sparse rts::hsgpCovariance::ZL_sparse(){
  sparse dummy;
  return dummy;
}

inline int rts::hsgpCovariance::Q(){
  return m * m * grid.T;
}

inline double rts::hsgpCovariance::log_likelihood(const VectorXd &u){
  double ll = 0;
  Z_chol(); //ADD PARAMETERS
  MatrixXd ar_factor_inverse = ar_factor.llt().solve(MatrixXd::Identity(grid.T,grid.T));
  MatrixXd umat(grid.N,grid.T);
  for(int t = 0; t< grid.T; t++){
    umat.col(t) = u.segment(t*grid.N,grid.N);
  }
  MatrixXd vmat = umat * ar_factor_inverse;
  double logdet = log_determinant();
  VectorXd uquad(m*m);
  VectorXd vquad(m*m);
  VectorXd uphi(m*m);
  VectorXd vphi(m*m);
  for(int t = 0; t< grid.T; t++){
    uphi = Phi.transpose() * umat.col(t);
    vphi = Phi.transpose() * vmat.col(t);
    uquad = glmmr::algo::forward_sub(L,uphi,m*m);
    vquad = glmmr::algo::forward_sub(L,vphi,m*m);
    ll += (-0.5*grid.N * log(2*M_PI) - 0.5*(umat.col(t).transpose() * vmat.col(t) - uquad.transpose()*vquad)(0));
  }
  ll += 0.5*logdet;
  return -1.0*ll;
}

inline double rts::hsgpCovariance::log_determinant(){
  double logdet = (grid.N - (m*m))*log(this->parameters_[0]);
  for(int i = 0; i < (m*m); i++){
    logdet += 2*log(L(i,i));
    logdet += log(Lambda(i));
  }
  
  logdet *= grid.T;
  double logdet_ar = 0;
  if(grid.T > 1){
    for(int t = 0; t < grid.T; t++) logdet_ar += 2*log(ar_factor_chol(t,t));
    logdet_ar *= grid.N;
  }
  return logdet + logdet_ar;
}

inline void rts::hsgpCovariance::update_rho(const double rho_){
  rho = rho_;
  ar_factor.setConstant(1.0);
  if(grid.T > 1){
    for(int t = 0; t < (grid.T)-1; t++){
      for(int s = t+1; s < grid.T; s++){
        ar_factor(t,s) = pow(rho,s);
        ar_factor(s,t) = ar_factor(t,s);
      }
    }
  }
  rts::cholesky(ar_factor_chol, ar_factor);
}

inline void rts::hsgpCovariance::set_function(bool squared_exp){
  sq_exp = squared_exp;
}

inline void rts::hsgpCovariance::gen_indices(){
  int counter = 0;
  for(int i = 0; i < m; i++){
    for(int j = 0; j< m; j++){
      indices(counter,0) = i;
      indices(counter,1) = j;
      counter++;
    }
  }
}

inline void rts::hsgpCovariance::gen_phi_prod(){
  for(int i = 0; i < (m*m); i++){
    ArrayXd phi = phi_nD(i);
    Phi.col(i) = phi.matrix();
  }
  PhiT = Phi.transpose() * Phi;
}

inline void rts::hsgpCovariance::Z_chol(){
  MatrixXd pnew = PhiT;
  for(int i = 0; i < (m*m); i++){
    pnew(i,i) += sqrt(Lambda(i));
  }
  rts::cholesky(L, pnew);
}

inline MatrixXd rts::hsgpCovariance::PhiSPD(){
  ArrayXXd pnew = Phi.array();
  for(int i = 0; i < (m*m); i++){
    pnew.col(i) *= sqrt(Lambda(i));
  }
  return pnew.matrix();
}

inline void rts::hsgpCovariance::update_lambda(){
  for(int i = 0; i < (m*m); i++){
    Lambda(i) = spd_nD((lambda_nD(i)).array().sqrt().matrix());
  }
}
