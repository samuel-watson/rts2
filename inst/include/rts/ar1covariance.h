#pragma once

#include <glmmr/covariance.hpp>
#include "griddata.h"
#include "rtsmaths.h"

namespace rts {

using namespace Eigen;
using namespace glmmr;

class ar1Covariance : public Covariance {
public:
  double rho;
  rts::griddata grid;
  
  ar1Covariance(const str& formula,const ArrayXXd &data,const strvec& colnames, int T);  
  ar1Covariance(const rts::ar1Covariance& cov);
  
  MatrixXd  ZL() override;
  MatrixXd  LZWZL(const VectorXd& w) override;
  MatrixXd  ZLu(const MatrixXd& u) override;
  MatrixXd  Lu(const MatrixXd& u) override;
  sparse    ZL_sparse() override;
  int       Q() const override;
  double    log_likelihood(const VectorXd &u) override;
  double    log_determinant() override;
  void      update_rho(const double rho_);
  void      update_grid(int T);
  void      update_parameters(const dblvec& parameters) override;
  void      update_parameters(const ArrayXd& parameters) override;
  void      update_parameters_extern(const dblvec& parameters) override;
  MatrixXd  ar_matrix(bool chol = false);
  VectorXd  log_gradient(const MatrixXd& u, double& logl);
  VectorXd  log_gradient_rho(const MatrixXd& u);
  
protected:
  MatrixXd L;
  MatrixXd ar_factor;    
  MatrixXd ar_factor_chol;
  MatrixXd ar_factor_deriv;  
  MatrixXd ar_factor_inverse;  
};

}


inline rts::ar1Covariance::ar1Covariance(const str& formula,
                const ArrayXXd &data,
                const strvec& colnames, int T) : Covariance(formula, data, colnames), 
                  grid(data, T), L(data.rows(),data.rows()),
                  ar_factor(T,T), ar_factor_chol(T,T), ar_factor_deriv(T,T), ar_factor_inverse(T,T) { 
      isSparse = false;
      update_rho(0.1);
    };
  
  inline rts::ar1Covariance::ar1Covariance(const rts::ar1Covariance& cov) : Covariance(cov.form_, cov.data_, cov.colnames_), grid(cov.grid), 
    L(cov.L), ar_factor(cov.ar_factor), ar_factor_chol(cov.ar_factor_chol), ar_factor_deriv(cov.ar_factor_deriv), ar_factor_inverse(cov.ar_factor_inverse) {
      isSparse = false; 
      update_rho(cov.rho);
  };

inline MatrixXd rts::ar1Covariance::ar_matrix(bool chol)
{
  if(chol){
    return ar_factor_chol;
  } else {
    return ar_factor;
  }
}

inline void rts::ar1Covariance::update_grid(int T)
{
  grid.setup(data_,T);
}

inline void rts::ar1Covariance::update_parameters_extern(const dblvec& parameters)
{
  parameters_ = parameters;
  update_parameters_in_calculators();
  L = Covariance::D(true,false);
};

inline void rts::ar1Covariance::update_parameters(const dblvec& parameters)
{
  parameters_ = parameters;
  update_parameters_in_calculators();
  L = Covariance::D(true,false);
};

inline void rts::ar1Covariance::update_parameters(const ArrayXd& parameters)
{
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
  L = Covariance::D(true,false);
};

inline MatrixXd rts::ar1Covariance::ZL()
{
  MatrixXd ZL = rts::kronecker(ar_factor_chol, L);
  return ZL;
}

inline MatrixXd rts::ar1Covariance::LZWZL(const VectorXd& w)
{
  MatrixXd ZL = rts::ar1Covariance::ZL();
  MatrixXd LZWZL = ZL.transpose() * w.asDiagonal() * ZL;
  LZWZL += MatrixXd::Identity(LZWZL.rows(), LZWZL.cols());
  return LZWZL;
}

inline MatrixXd rts::ar1Covariance::ZLu(const MatrixXd& u)
{
  MatrixXd ZLu = rts::ar1Covariance::ZL() * u;
  return ZLu;
}

inline MatrixXd rts::ar1Covariance::Lu(const MatrixXd& u)
{
  MatrixXd ZLu = rts::ar1Covariance::ZL() * u;
  return ZLu;
}

inline sparse rts::ar1Covariance::ZL_sparse()
{
  sparse dummy;
  return dummy;
}

inline int rts::ar1Covariance::Q() const 
{
  return grid.N * grid.T;
}

inline double rts::ar1Covariance::log_likelihood(const VectorXd &u)
{
  // need to collapse u to v
  double ll = 0;
  MatrixXd ar_factor_inverse = ar_factor.llt().solve(MatrixXd::Identity(grid.T,grid.T));
  MatrixXd umat(grid.N,grid.T);
  for(int t = 0; t< grid.T; t++)
  {
    umat.col(t) = u.segment(t*grid.N,grid.N);
  }
  MatrixXd vmat = umat * ar_factor_inverse;
  double logdet = log_determinant();
  VectorXd uquad(grid.N);
  VectorXd vquad(grid.N);
  for(int t = 0; t< grid.T; t++)
  {
    uquad = glmmr::algo::forward_sub(L,umat.col(t),grid.N);
    vquad = glmmr::algo::forward_sub(L,vmat.col(t),grid.N);
    ll += (-0.5*grid.N * log(2*M_PI) - 0.5*uquad.transpose()*vquad);
  }
  ll += 0.5*logdet;
  return -1.0*ll;
}

// this needs to be updated to the version in glmmrBase which is faster
inline VectorXd rts::ar1Covariance::log_gradient(const MatrixXd& u, double& logl)
{
  std::vector<MatrixXd> derivs;
  derivatives(derivs,1);
  VectorXd grad(derivs.size()-1);
  grad.setZero();
  MatrixXd Sinv = D();
  Sinv = Sinv.llt().solve(MatrixXd::Identity(Sinv.rows(),Sinv.cols()));
  int niter = u.cols();
  MatrixXd ar_factor_inverse = ar_factor.llt().solve(MatrixXd::Identity(grid.T,grid.T));
  logl = 0;
#pragma omp parallel for reduction(+:logl)
  for(int j = 0; j < niter; j++) logl += log_likelihood(u.col(j));
  logl *= 1.0 / (double) niter;
  for(int i = 1; i < derivs.size(); i++)
  {
    MatrixXd Sdi = Sinv * derivs[i];
    double dlogdet = Sdi.trace();
    grad(i-1) = -grid.T*dlogdet;
    MatrixXd dinv = Sdi * Sinv;
    double val = 0;
  #pragma omp parallel for reduction(+:val)
    for(int j = 0; j < niter; j++)
    {
      MatrixXd umat(grid.N,grid.T);
      MatrixXd vmat(grid.N,grid.T);
      for(int t = 0; t < grid.T; t++)
      {
        umat.col(t) = u.col(j).segment(t*grid.N, grid.N);
      }
      vmat = umat * ar_factor_inverse;
      for(int t = 0; t < grid.T; t++)
      {
        val += ((umat.col(t).transpose()*dinv)*vmat.col(t))(0);
      }           
    }      
    grad(i-1) += val/((double)niter);
    grad(i-1) *= 0.5;
  }
  return grad;
}

inline VectorXd rts::ar1Covariance::log_gradient_rho(const MatrixXd& u)
{
  VectorXd grad(1);
  grad.setZero();
  MatrixXd Sinv = D();
  Sinv = Sinv.llt().solve(MatrixXd::Identity(Sinv.rows(),Sinv.cols()));
  int niter = u.cols();
  MatrixXd ar_factor_inverse = ar_factor.llt().solve(MatrixXd::Identity(grid.T,grid.T));
  MatrixXd umat(grid.N,grid.T);
  MatrixXd vmat(grid.N,grid.T);
  MatrixXd Sdi = ar_factor_inverse * ar_factor_deriv;
  double dlogdet = Sdi.trace();
  grad(0) = -grid.N*dlogdet;
  MatrixXd dinv = Sdi * ar_factor_inverse;
  double val = 0;

#pragma omp parallel for reduction(+:val) private(umat,vmat)
  for(int j = 0; j < niter; j++)
  {
    for(int t = 0; t< grid.T; t++)
    {
      umat.col(t) = u.col(j).segment(t*grid.N,grid.N);
    }
    vmat = umat * dinv;
    for(int t = 0; t < grid.T; t++)
    {
      val += ((umat.col(t).transpose()*Sinv)*vmat.col(t))(0);
    }           
  }      
  grad(0) += val/((double)niter);
  grad(0) *= 0.5;
  return grad;
}

inline double rts::ar1Covariance::log_determinant()
{
  double logdet = Covariance::log_determinant();
  logdet *= grid.T;
  double logdet_ar = 0;
  if(grid.T > 1){
    for(int t = 0; t < grid.T; t++) logdet_ar += 2*log(ar_factor_chol(t,t));
    logdet_ar *= grid.N;
  }
  return logdet + logdet_ar;
}

inline void rts::ar1Covariance::update_rho(const double rho_)
{
  rho = rho_;
  ar_factor.setConstant(1.0);
  if(grid.T > 1){
    for(int t = 0; t < (grid.T)-1; t++){
      for(int s = t+1; s < grid.T; s++){
        ar_factor(t,s) = pow(rho,s);
        ar_factor_deriv(t,s) = s*pow(rho,s-1);
        ar_factor(s,t) = ar_factor(t,s);
        ar_factor_deriv(s,t) = ar_factor_deriv(t,s);
      }
    }
  }
  //ar_factor_chol = rts::cholesky(ar_factor);
  ar_factor_chol = MatrixXd(ar_factor.llt().matrixL());
  ar_factor_inverse = ar_factor.llt().solve(MatrixXd::Identity(grid.T, grid.T));
}


