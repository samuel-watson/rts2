#pragma once

#include <memory>
#include <glmmr/covariance.hpp>
#include "griddata.h"
#include "rtsmaths.h"

namespace rts {

using namespace Eigen;
using namespace glmmr;

class nngpCovariance : public Covariance {
public:
  double          rho;
  rts::griddata   grid;
  MatrixXd        A;
  VectorXd        Dvec;
  int             m;
  
  nngpCovariance(const str& formula,const ArrayXXd &data,const strvec& colnames,const dblvec& parameters,int T, int m_,const rts::griddata& grid_);  
  nngpCovariance(const str& formula,const ArrayXXd &data,const strvec& colnames,int T, int m_,const rts::griddata& grid_);  
  nngpCovariance(const rts::nngpCovariance& cov);
  
  MatrixXd        inv_ldlt_AD(const MatrixXd &A, const VectorXd &D, const ArrayXXi &NN);
  MatrixXd        D(bool chol = true, bool upper = false) override;
  MatrixXd        ZL() override;
  MatrixXd        LZWZL(const VectorXd& w) override;
  MatrixXd        ZLu(const MatrixXd& u) override;
  MatrixXd        Lu(const MatrixXd& u) override;
  sparse          ZL_sparse() override;
  int             Q() const override;
  double          log_likelihood(const VectorXd &u) override;
  double          log_determinant() override;
  void            update_rho(const double rho_);
  void            gen_AD();
  void            update_parameters(const dblvec& parameters) override;
  void            update_parameters(const ArrayXd& parameters) override;
  void            update_parameters_d(const ArrayXd& parameters);
  void            update_parameters_extern(const dblvec& parameters) override;
  VectorMatrix    submatrix(int i);
  void            set_function(bool squared_exp);
  MatrixXd        ar_matrix(bool chol = false);
  // generates AD and the derivatives
  void            gen_AD_derivatives(VectorXd& D1, VectorXd& D2, MatrixXd& A1, MatrixXd& A2); 
  VectorXd        log_gradient(const MatrixXd& u, double& ll);
  VectorXd        log_gradient_rho(const MatrixXd& u);
  
protected:
  MatrixXd        ar_factor;    
  MatrixXd        ar_factor_chol;
  MatrixXd        ar_factor_inverse;
  bool            sq_exp = false;
  MatrixXd        ar_factor_deriv;
};
}


inline rts::nngpCovariance::nngpCovariance(const str& formula,
                 const ArrayXXd &data,
                 const strvec& colnames,
                 const dblvec& parameters,
                 int T, int m_,
                 const rts::griddata& grid_) : Covariance(formula, data, colnames, parameters),  
                   grid(grid_),
                   A(m_,data.rows()), Dvec(data.rows()), m(m_), ar_factor(T,T), ar_factor_chol(T,T), ar_factor_deriv(T,T) {
    isSparse = false;
    grid.genNN(m);
    gen_AD();
    update_rho(0.1);
  }
  
  inline rts::nngpCovariance::nngpCovariance(const str& formula,
                 const ArrayXXd &data,
                 const strvec& colnames,
                 int T, int m_,
                 const rts::griddata& grid_) : Covariance(formula, data, colnames),  
                 grid(grid_),
                 A(m_,data.rows()), Dvec(data.rows()), m(m_), ar_factor(T,T), ar_factor_chol(T,T), ar_factor_deriv(T,T) {
    isSparse = false;
    grid.genNN(m);
    update_rho(0.1);
  }
  
  inline rts::nngpCovariance::nngpCovariance(const rts::nngpCovariance& cov) : Covariance(cov.form_, cov.data_, cov.colnames_, cov.parameters_),  
    grid(cov.grid), A(grid.m,grid.N), Dvec(grid.N), m(cov.m), ar_factor(grid.T,grid.T), ar_factor_chol(grid.T,grid.T),
    ar_factor_deriv(grid.T,grid.T) {
      isSparse = false;
      grid.genNN(m);
      gen_AD();
      update_rho(cov.rho);
    }

inline MatrixXd rts::nngpCovariance::ar_matrix(bool chol)
{
  if(chol){
    return ar_factor_chol;
  } else {
    return ar_factor;
  }
}

inline void rts::nngpCovariance::set_function(bool squared_exp)
{
  sq_exp = squared_exp;
}

inline MatrixXd rts::nngpCovariance::inv_ldlt_AD(const MatrixXd &A, 
                            const VectorXd &D,
                            const ArrayXXi &NN)
{
  int n = A.cols();
  int m = A.rows();
  MatrixXd y = MatrixXd::Zero(n,n);
  ArrayXd dsqrt = Dvec.array().sqrt();
#pragma omp parallel for  
  for(int k=0; k<n; k++){
    int idxlim;
    for (int i = k; i < n; i++) {
      idxlim = i<=m ? i : m;
      double lsum = 0;
      for (int j = 0; j < idxlim; j++) {
        lsum += A(j,i) * y(NN(j,i),k);
      }
      y(i,k) = i==k ? (1+lsum)  : lsum ;
    }
  }
  
  return y * dsqrt.matrix().asDiagonal();
}

inline MatrixXd rts::nngpCovariance::D(bool chol, bool upper)
{
  MatrixXd As = inv_ldlt_AD(A,Dvec,grid.NN);
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

inline MatrixXd rts::nngpCovariance::ZL()
{
  MatrixXd L = D(true,false);
  MatrixXd ZL = rts::kronecker(ar_factor_chol, L);
  return ZL;
}

inline MatrixXd rts::nngpCovariance::LZWZL(const VectorXd& w)
{
  MatrixXd ZL = rts::nngpCovariance::ZL();
  MatrixXd LZWZL = ZL.transpose() * w.asDiagonal() * ZL;
  LZWZL += MatrixXd::Identity(LZWZL.rows(), LZWZL.cols());
  return LZWZL;
}

inline MatrixXd rts::nngpCovariance::ZLu(const MatrixXd& u)
{
  MatrixXd L = D(true,false);
  MatrixXd ZL = rts::kronecker(ar_factor_chol, L);
  return ZL * u;
}

inline MatrixXd rts::nngpCovariance::Lu(const MatrixXd& u)
{
  return rts::nngpCovariance::ZLu(u);
}

inline sparse rts::nngpCovariance::ZL_sparse()
{
  sparse dummy;
  return dummy;
}

inline int rts::nngpCovariance::Q() const 
{
  return grid.N * grid.T;
}

inline double rts::nngpCovariance::log_likelihood(const VectorXd &u)
{
  double ll1 = 0.0;
  double logdet = log_determinant();
  int idxlim;
  double au,av; 
  if(grid.T > 1)
  {
    // requires copying a lot of data, is there a better way of doing this? 
    MatrixXd umat(grid.N,grid.T);
    for(int t = 0; t< grid.T; t++){
      umat.col(t) = u.segment(t*grid.N,grid.N);
    }
    MatrixXd vmat = umat * ar_factor_inverse;
    for(int t = 0; t < grid.T; t++)
    {
      double qf = umat(0,t)*vmat(0,t)/Dvec(0);
      for(int i = 1; i < grid.N; i++)
      {
        idxlim = i <= m ? i : m;
        VectorXd usec(idxlim);
        VectorXd vsec(idxlim);
        for(int j = 0; j < idxlim; j++) 
        {
          usec(j) = umat(grid.NN(j,i),t);
          vsec(j) = vmat(grid.NN(j,i),t);
        }
        au = umat(i,t) - (A.col(i).segment(0,idxlim).transpose() * usec)(0);
        av = vmat(i,t) - (A.col(i).segment(0,idxlim).transpose() * vsec)(0);
        qf += au*av/Dvec(i);
      }
      ll1 -= 0.5*qf; 
    }
  } else {
    double qf = u(0)*u(0)/Dvec(0);
    for(int i = 1; i < grid.N; i++)
    {
      idxlim = i <= m ? i : m;
      VectorXd usec(idxlim);
      for(int j = 0; j < idxlim; j++) 
      {
        usec(j) = u(grid.NN(j,i));
      }
      au = u(i) - (A.col(i).segment(0,idxlim).transpose() * usec)(0);
      qf += au*au/Dvec(i);
    }
    ll1 -= 0.5*qf; 
  }
  ll1 -= 0.5*logdet  + 0.5*grid.N*grid.T*log(2*M_PI);
  return ll1;
}


inline double rts::nngpCovariance::log_determinant()
{
  double logdet = Dvec.array().log().sum();
  logdet *= grid.T;
  double logdet_ar = 0;
  if(grid.T > 1){
    for(int t = 0; t < grid.T; t++) logdet_ar += 2*log(ar_factor_chol(t,t));
    logdet_ar *= grid.N;
  }
  return logdet + logdet_ar;
}

inline void rts::nngpCovariance::update_rho(const double rho_)
{
  rho = rho_;
  ar_factor.setConstant(1.0);
  if(grid.T > 1){
    for(int t = 0; t < (grid.T)-1; t++){
      for(int s = t+1; s < grid.T; s++){
        ar_factor(t,s) = pow(rho,s);
        ar_factor_deriv(t,s) = s * pow(rho,s-1);
        ar_factor(s,t) = ar_factor(t,s);
        ar_factor_deriv(s,t) = ar_factor_deriv(t,s);
      }
    }
  }
  ar_factor_chol = MatrixXd(ar_factor.llt().matrixL());
  ar_factor_inverse = ar_factor.llt().solve(MatrixXd::Identity(grid.T,grid.T));
}

inline void rts::nngpCovariance::gen_AD()
{
  A.setZero();
  Dvec.setZero();
  int idxlim;
  double val = this->parameters_[0];
  Dvec(0) = val;
  double tmp;
  
#pragma omp parallel for
  for(int i = 1; i < grid.N; i++)
  {
    idxlim = i <= m ? i : m;
    MatrixXd S(idxlim,idxlim);
    VectorXd Sv(idxlim);
    for(int j = 0; j<idxlim; j++)
    {
      S(j,j) = val;
    }
    if(idxlim > 1)
    {
      for(int j = 0; j<(idxlim-1); j++)
      {
        for(int k = j+1; k<idxlim; k++)
        {
          tmp = this->calc_[0].get_covariance_data(grid.NN(j,i), grid.NN(k,i), 0);
          if(sq_exp){
            S(j,k) = this->parameters_[0] * exp(-1 * tmp * tmp / (this->parameters_[1] * this->parameters_[1]));//Covariance::get_val(0,grid.NN(j,i),grid.NN(k,i));
          } else {
            S(j,k) = this->parameters_[0] * exp(-1 * tmp / this->parameters_[1]);
          }
          S(k,j) = S(j,k);
        }
      }
    }
    for(int j = 0; j<idxlim; j++)
    {
      tmp = this->calc_[0].get_covariance_data(i, grid.NN(j,i), 0); // this->dists[0]((grid.N-1)*grid.NN(j,i) - ((grid.NN(j,i)-1)*grid.NN(j,i)/2) + (i-grid.NN(j,i)-1),0)/ this->parameters_[1];
      Sv(j) = sq_exp ? this->parameters_[0] * exp(-1 * tmp * tmp / (this->parameters_[1] * this->parameters_[1])) : this->parameters_[0] * exp(-1 * tmp / this->parameters_[1]);//Covariance::get_val(0,i,grid.NN(j,i));
    }
    A.block(0,i,idxlim,1) = S.ldlt().solve(Sv);
    Dvec(i) = val - (A.col(i).segment(0,idxlim).transpose() * Sv)(0);
  }
}

inline void rts::nngpCovariance::gen_AD_derivatives(VectorXd& D1, VectorXd& D2, MatrixXd& A1, MatrixXd& A2)
{
  A.setZero();
  Dvec.setZero();
  A1.setZero();
  A2.setZero();
  D1.setZero();
  D2.setZero();
  int idxlim;
  double tmp;
  double val = this->parameters_[0];
  Dvec(0) = val;
  D1(0) = 1;
  D2(0) = 0;

  #pragma omp parallel for private(idxlim, tmp)
  for(int i = 1; i < grid.N; i++)
  {
    idxlim = i <= m ? i : m;
    MatrixXd S(idxlim,idxlim);
    MatrixXd Sinv(idxlim,idxlim);
    VectorXd Sv(idxlim);
    MatrixXd dS1(idxlim,idxlim);
    MatrixXd dS2(idxlim,idxlim);
    VectorXd dSv1(idxlim);
    VectorXd dSv2(idxlim);
    double p3 = this->parameters_[1] * this->parameters_[1] * this->parameters_[1];
    for(int j = 0; j<idxlim; j++)
    {
      S(j,j) = val;
      dS1(j,j) = 1;
      dS2(j,j) = 0;
    }
    if(idxlim > 1)
    {
      for(int j = 0; j<(idxlim-1); j++)
      {
        for(int k = j+1; k<idxlim; k++)
        {
          tmp = this->calc_[0].get_covariance_data(grid.NN(j,i), grid.NN(k,i), 0);
          if(sq_exp){
            S(j,k) = this->parameters_[0] * exp(-1 * tmp * tmp / (this->parameters_[1] * this->parameters_[1]));
            dS1(j,k) = exp(-1 * tmp * tmp / (this->parameters_[1] * this->parameters_[1]));            
            dS2(j,k) = 2 * this->parameters_[0] * tmp * tmp * exp(-1 * tmp * tmp / (this->parameters_[1] * this->parameters_[1])) / p3;
          } else {
            S(j,k) = this->parameters_[0] * exp(-1 * tmp / this->parameters_[1]);
            dS1(j,k) = exp(-1 * tmp / this->parameters_[1]);
            dS2(j,k) = this->parameters_[0] * tmp * exp(-1 * tmp / this->parameters_[1]) / (this->parameters_[1] * this->parameters_[1]) ;
          }
          S(k,j) = S(j,k);
          dS1(k,j) = dS1(j,k);
          dS2(k,j) = dS2(j,k);
        }
      }
    }
    for(int j = 0; j<idxlim; j++)
    {
      tmp = this->calc_[0].get_covariance_data(i, grid.NN(j,i), 0); 
      Sv(j) = sq_exp ? this->parameters_[0] * exp(-1 * tmp * tmp / (this->parameters_[1] * this->parameters_[1])) : this->parameters_[0] * exp(-1 * tmp / this->parameters_[1]);
      dSv1(j) = sq_exp ? exp(-1 * tmp * tmp / (this->parameters_[1] * this->parameters_[1])) : exp(-1 * tmp / this->parameters_[1]);
      dSv2(j) = sq_exp ? 2 * this->parameters_[0] * tmp * tmp * exp(-1 * tmp * tmp / (this->parameters_[1] * this->parameters_[1])) / p3 : this->parameters_[0] * tmp * exp(-1 * tmp / this->parameters_[1]) / (this->parameters_[1] * this->parameters_[1]);
    }
    A.block(0,i,idxlim,1) = S.ldlt().solve(Sv);
    A1.block(0,i,idxlim,1) = S.ldlt().solve(dSv1 - dS1 * A.col(i).segment(0,idxlim)); 
    A2.block(0,i,idxlim,1) = S.ldlt().solve(dSv2 - dS2 * A.col(i).segment(0,idxlim)); 
    Dvec(i) = val - (A.col(i).segment(0,idxlim).transpose() * Sv)(0);
    D1(i) = 1 - (dSv1.transpose()*A.col(i).segment(0,idxlim))(0) - (Sv.transpose() * A1.col(i).segment(0,idxlim))(0);
    D2(i) =  -1.0 * (dSv2.transpose() * A.col(i).segment(0,idxlim))(0) - (Sv.transpose() * A2.col(i).segment(0,idxlim))(0);
  }
}

inline VectorXd rts::nngpCovariance::log_gradient(const MatrixXd& umat, double& ll)
{
  VectorXd grad(2);
  grad.setZero();
  MatrixXd A1(m,grid.N);
  MatrixXd A2(m,grid.N);
  VectorXd D1(grid.N);
  VectorXd D2(grid.N);
  gen_AD_derivatives(D1,D2,A1,A2);

  //log determinant derivatives
  double dlogdet1 = 0;
  double dlogdet2 = 0;
  double logdet = log_determinant();
  
  for(int i = 0; i < grid.N; i++)
  {
    dlogdet1 += D1(i)/Dvec(i);
    dlogdet2 += D2(i)/Dvec(i);
    D1(i) *= (1 / (Dvec(i) * Dvec(i)));
    D2(i) *= (1 / (Dvec(i) * Dvec(i)));
  } 
  dlogdet1 *= grid.T;
  dlogdet2 *= grid.T;

  double au, av, dau1, dav1, dau2, dav2, qf1, qf2, qf0;
  double ll1 = 0.0;
  double ll2 = 0.0;
  int niter = umat.cols();
  int idxlim;
  ll = 0;

//#pragma omp parallel for reduction(+:ll1,ll2) private(au, av, dau1, dav1, dau2, dav2, qf1, qf2, idxlim) if(niter > 50)
  for(int k = 0; k < niter; k++)
  {
    if(grid.T > 1)
    {
      VectorXd v(grid.N);    
      for(int t = 0; t < grid.T; t++)
      {
        v = ar_factor_inverse * umat.col(k).segment(t*grid.N,grid.N);
        qf0 = umat(t*grid.N,k)*v(t*grid.N)/Dvec(0);
        qf1 = -1.0 * umat(t*grid.N,k)*v(t*grid.N)*D1(0);
        qf2 = -1.0 * umat(t*grid.N,k)*v(t*grid.N)*D2(0);
        for(int i = 1; i < grid.N; i++)
        {
          idxlim = i <= m ? i : m;
          VectorXd usec(idxlim);
          VectorXd vsec(idxlim);
          for(int j = 0; j < idxlim; j++) 
          {
            usec(j) = umat(grid.NN(j,i)+t*grid.N,k);
            vsec(j) = v(grid.NN(j,i)+t*grid.N);
          }
          au = umat(i+t*grid.N,k) - (A.col(i).segment(0,idxlim).transpose() * usec)(0);
          dau1 = -1.0 * (A1.col(i).segment(0,idxlim).transpose() * usec)(0);
          dau2 = -1.0 * (A2.col(i).segment(0,idxlim).transpose() * usec)(0);
          av = v(i+t*grid.N) - (A.col(i).segment(0,idxlim).transpose() * vsec)(0);
          dav1 = -1.0 * (A1.col(i).segment(0,idxlim).transpose() * vsec)(0);
          dav2 = -1.0 * (A2.col(i).segment(0,idxlim).transpose() * vsec)(0);
          qf1 += dau1*av/Dvec(i) - au*av*D1(i) + au*dav1/Dvec(i);
          qf2 += dau2*av/Dvec(i) - au*av*D2(i) + au*dav2/Dvec(i);
          qf0 += au*av/Dvec(i);
        }
        ll1 -= 0.5*qf1;
        ll2 -= 0.5*qf2; 
        ll -= 0.5*qf0;
      }
    } else {
      qf0 = umat(0)*umat(0)/Dvec(0);
      qf1 = -1.0 * umat(0,k)*umat(0,k)*D1(0);
      qf2 = -1.0 * umat(0,k)*umat(0,k)*D2(0);
      for(int i = 1; i < grid.N; i++)
      {
        idxlim = i <= m ? i : m;
        VectorXd usec(idxlim);
        for(int j = 0; j < idxlim; j++) 
        {
          usec(j) = umat(grid.NN(j,i),k);
        }
        au = umat(i,k) - (A.col(i).segment(0,idxlim).transpose() * usec)(0);
        dau1 = -1.0 * (A1.col(i).segment(0,idxlim).transpose() * usec)(0);
        dau2 = -1.0 * (A2.col(i).segment(0,idxlim).transpose() * usec)(0);
        qf1 += 2*dau1*au/Dvec(i) - au*au*D1(i);
        qf2 += 2*dau2*au/Dvec(i) - au*au*D2(i);
        qf0 += au*au/Dvec(i);
      }
      ll1 -= 0.5*qf1; 
      ll2 -= 0.5*qf2; 
      ll -= 0.5*qf0;
    }
  }  

  grad(0) = -0.5 * dlogdet1 + ll1 / (double)niter;
  grad(1) = -0.5 * dlogdet2 + ll2 / (double)niter;
  ll *= 1.0/(double)niter;
  ll -= 0.5*logdet + 0.5*grid.N*grid.T*log(2*M_PI);

  return grad;
}

inline VectorXd rts::nngpCovariance::log_gradient_rho(const MatrixXd& u)
{
  VectorXd grad(1);
  grad.setZero();

  //log determinant derivatives
  double dlogdet1 = 0;
  MatrixXd ar_factor_deriv_chol(ar_factor_deriv.llt().matrixL());  
  for(int i = 0; i < grid.T; i++)
  {
    dlogdet1 += 2 * log(ar_factor_deriv_chol(i,i));
  } 
  dlogdet1 *= grid.N;

  double au, av, qf1;
  double ll1 = 0.0;
  int niter = u.cols();
  int idxlim;
  MatrixXd ar_mult = ar_factor_inverse * ar_factor_deriv * ar_factor_inverse;
  
#pragma omp parallel for reduction(+:ll1) private(au, av, qf1, idxlim) if(niter > 50)
  for(int k = 0; k < niter; k++)
  {
    VectorXd v(grid.N);    
    for(int t = 0; t < grid.T; t++)
    {
      v = ar_factor_inverse * u.col(k).segment(t*grid.N,grid.N);
      double qf = u(t*grid.N,k)*v(t*grid.N)/Dvec(0);
      for(int i = 1; i < grid.N; i++)
      {
        idxlim = i <= m ? i : m;
        VectorXd usec(idxlim);
        VectorXd vsec(idxlim);
        for(int j = 0; j < idxlim; j++) 
        {
          usec(j) = u(grid.NN(j,i)+t*grid.N,k);
          vsec(j) = v(grid.NN(j,i)+t*grid.N);
        }
        au = u(i+t*grid.N,k) - (A.col(i).segment(0,idxlim).transpose() * usec)(0);
        av = v(i+t*grid.N) - (A.col(i).segment(0,idxlim).transpose() * vsec)(0);
        qf += au*av/Dvec(i);
      }
      ll1 -= 0.5*qf; 
    }
  }

  grad(0) = -0.5 * dlogdet1 + ll1 / (double) niter;  
  return grad;
}

inline VectorMatrix rts::nngpCovariance::submatrix(int i)
{
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
  VectorMatrix result(idxlim);
  result.vec = Sv;
  result.mat = S;
  return result;
}

inline void rts::nngpCovariance::update_parameters(const dblvec& parameters)
{
  parameters_ = parameters;
  update_parameters_in_calculators();
  gen_AD();
}

inline void rts::nngpCovariance::update_parameters_extern(const dblvec& parameters)
{
  parameters_ = parameters;
  update_parameters_in_calculators();
  gen_AD();
}

inline void rts::nngpCovariance::update_parameters(const ArrayXd& parameters)
{
  if(parameters_.size()==0)
  {
    for(unsigned int i = 0; i < parameters.size(); i++)parameters_.push_back(parameters(i));
    update_parameters_in_calculators();
  } else if(parameters_.size() == parameters.size())
  {
    for(unsigned int i = 0; i < parameters.size(); i++)parameters_[i] = parameters(i);
    update_parameters_in_calculators();
  } 
  gen_AD();
};

inline void rts::nngpCovariance::update_parameters_d(const ArrayXd& parameters)
{
  if(parameters_.size()==0)
  {
    for(unsigned int i = 0; i < parameters.size(); i++) parameters_.push_back(parameters(i));
    update_parameters_in_calculators();
  } else if(parameters_.size() == parameters.size())
  {
    for(unsigned int i = 0; i < parameters.size(); i++) parameters_[i] = parameters(i);
    update_parameters_in_calculators();
  } else {
    throw std::runtime_error("NNGP received the wrong number of parameters");
  }
};


