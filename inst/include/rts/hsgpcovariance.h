#pragma once

#include <glmmr/covariance.hpp>
#include "griddata.h"
#include "rtsmaths.h"

namespace rts {

using namespace Eigen;
using namespace glmmr;

class hsgpCovariance : public Covariance {
public:
  double          rho;
  rts::griddata   grid;
  int             m;
  ArrayXd         L_boundary;

  hsgpCovariance(const str& formula,const ArrayXXd &data,const strvec& colnames,int T,int m_,const ArrayXd& L_boundary_);
  hsgpCovariance(const rts::hsgpCovariance& cov);

  double                spd_nD(int i);
  double                d_spd_nD(int i, int par, bool sqrt_lambda = true);
  ArrayXd               phi_nD(int i);
  MatrixXd              ZL() override;
  MatrixXd              ZL_deriv(int par, bool theta = true);
  MatrixXd              D(bool chol = true, bool upper = false) override;
  MatrixXd              LZWZL(const VectorXd& w) override;
  MatrixXd              ZLu(const MatrixXd& u) override;
  MatrixXd              Lu(const MatrixXd& u) override;
  sparse                ZL_sparse() override;
  int                   Q() const override;
  double                log_likelihood(const VectorXd &u) override;
  double                log_determinant() override;
  void                  update_rho(const double rho_);
  void                  update_grid(int T);
  void                  update_parameters(const dblvec& parameters) override;
  void                  update_parameters(const ArrayXd& parameters) override;
  void                  update_parameters_extern(const dblvec& parameters) override;
  void                  set_function(bool squared_exp);
  MatrixXd              PhiSPD(bool lambda = true, bool inverse = false);
  ArrayXd               LambdaSPD();
  MatrixXd              ar_matrix(bool chol = false);

protected:
  MatrixXd    L; // cholesky decomposition of Lambda + PhiTPhi m^2 * m^2
  ArrayXd     Lambda;
  MatrixXd    ar_factor;
  MatrixXd    ar_factor_chol;
  MatrixXd    ar_factor_inverse;
  MatrixXd    ar_factor_deriv;
  ArrayXXi    indices;
  MatrixXd    Phi;
  MatrixXd    PhiT;
  bool        sq_exp = false;
  void        gen_indices();
  void        gen_phi_prod();
  void        update_lambda();
};

}


inline rts::hsgpCovariance::hsgpCovariance(const str& formula,
                const ArrayXXd &data,
                const strvec& colnames,
                int T,
                int m_,
                const ArrayXd& L_boundary_) : Covariance(formula, data, colnames),
                grid(data, T), m(m_), L_boundary(L_boundary_),
                L(grid.N,m*m), Lambda(m*m),
                ar_factor(T,T), ar_factor_chol(T,T),ar_factor_inverse(T,T),ar_factor_deriv(T,T),
                indices(m*m,2), Phi(grid.N,m*m), PhiT(m*m,m*m) {
    gen_indices();
    gen_phi_prod();
    isSparse = false;
    update_rho(0.1);
  };

inline rts::hsgpCovariance::hsgpCovariance(const rts::hsgpCovariance& cov) : Covariance(cov.form_, cov.data_, cov.colnames_), grid(cov.grid),
  L_boundary(cov.L_boundary), L(cov.L), Lambda(cov.Lambda), ar_factor(cov.ar_factor), ar_factor_chol(cov.ar_factor_chol),
  ar_factor_inverse(cov.ar_factor_inverse), ar_factor_deriv(cov.ar_factor_deriv), indices(cov.indices), Phi(cov.Phi), PhiT(cov.PhiT) {
    isSparse = false;
    update_rho(cov.rho);
};

inline MatrixXd rts::hsgpCovariance::ar_matrix(bool chol)
{
  if(chol){
    return ar_factor_chol;
  } else {
    return ar_factor;
  }
}

inline double rts::hsgpCovariance::spd_nD(int i)
{
  Array2d w;
  w(0) = (indices(i,0)*M_PI)/(2*L_boundary(0));
  w(1) = (indices(i,1)*M_PI)/(2*L_boundary(1));
  w(0) = w(0)*w(0);
  w(1) = w(1)*w(1);
  double S;
  double phisq = parameters_[1] * parameters_[1];
  if(sq_exp){
    S = parameters_[0] * 2 * M_PI * phisq * exp(-0.5 * phisq * (w(0) + w(1)));
  } else {
    double S1 = parameters_[0] * 4 * M_PI * phisq;
    double S2 = 1 + phisq * (w(0) + w(1));
    S = S1 * pow(S2,-1.5);
  }
  return S;
}

inline double rts::hsgpCovariance::d_spd_nD(int i, int par, bool sqrt_lambda)
{
  Array2d w;
  w(0) = (indices(i,0)*M_PI)/(2*L_boundary(0));
  w(1) = (indices(i,1)*M_PI)/(2*L_boundary(1));
  w(0) = w(0)*w(0);
  w(1) = w(1)*w(1);
  double S;
  double phisq = parameters_[1] * parameters_[1];
  if(sq_exp){
    if(par == 0){
      S = 2 * M_PI * phisq * exp(-0.5 * phisq * (w(0) + w(1)));
    } else {
      S = -1.0 * (w(0) + w(1)) * parameters_[0] * 2 * M_PI * phisq * exp(-0.5 * phisq * (w(0) + w(1)));
    }
  } else {
    double S2 = 1 + phisq * (w(0) + w(1));
    if(par == 0){
      S = 4 * M_PI * phisq * pow(S2,-1.5);
    } else {
      S = parameters_[0] * 4 * M_PI * (2 * parameters_[1] * pow(S2,-1.5) - 3 * (w(0) + w(1)) * phisq *  parameters_[1] * pow(S2,-2.5));
    }
  }

  if(sqrt_lambda)
  {
    double SS = spd_nD(i);
    S *= 0.5 / sqrt(SS);
  }

  return S;
}

inline ArrayXd rts::hsgpCovariance::phi_nD(int i)
{
  ArrayXd fi1(grid.N);
  ArrayXd fi2(grid.N);
  fi1 = (1/sqrt(L_boundary(0))) * sin(indices(i,0)*M_PI*(grid.X.col(0)+L_boundary(0))/(2*L_boundary(0)));
  fi2 = (1/sqrt(L_boundary(1))) * sin(indices(i,1)*M_PI*(grid.X.col(1)+L_boundary(1))/(2*L_boundary(1)));
  fi1 *= fi2;
  return fi1;
}

inline MatrixXd rts::hsgpCovariance::D(bool chol, bool upper)
{
  MatrixXd As = PhiSPD();
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

inline void rts::hsgpCovariance::update_grid(int T)
{
  grid.setup(data_,T);
}

inline void rts::hsgpCovariance::update_parameters_extern(const dblvec& parameters)
{
  parameters_ = parameters;
  update_lambda();
};

inline void rts::hsgpCovariance::update_parameters(const dblvec& parameters)
{
  parameters_ = parameters;
  update_lambda();
};

inline void rts::hsgpCovariance::update_parameters(const ArrayXd& parameters)
{
  if(parameters_.size()==0){
    for(unsigned int i = 0; i < parameters.size(); i++){
      parameters_.push_back(parameters(i));
    }
  } else {
    for(unsigned int i = 0; i < parameters.size(); i++){
      parameters_[i] = parameters(i);
    }
  }
  update_lambda();
};

inline MatrixXd rts::hsgpCovariance::ZL()
{
  MatrixXd ZL = rts::kronecker(ar_factor_chol, PhiSPD(true, false));
  return ZL;
}

inline MatrixXd rts::hsgpCovariance::ZL_deriv(int par, bool theta)
{
  if(theta)
  {
    ArrayXd Lambda_deriv(Lambda.size());
  #pragma omp parallel for
    for(int i = 0; i < (m*m); i++)
    {
      Lambda_deriv(i) = d_spd_nD(i,par,true);
    }
    MatrixXd pnew = Phi;
    pnew *= Lambda_deriv.matrix().asDiagonal();
    MatrixXd ZL = rts::kronecker(ar_factor_chol, pnew);
    return ZL;
  } else {
    MatrixXd ZL = rts::kronecker(ar_factor_deriv, PhiSPD());
    return ZL;
  }
}

inline MatrixXd rts::hsgpCovariance::LZWZL(const VectorXd& w)
{
  MatrixXd ZL = rts::hsgpCovariance::ZL();
  MatrixXd LZWZL = ZL.transpose() * w.asDiagonal() * ZL;
  LZWZL += MatrixXd::Identity(LZWZL.rows(), LZWZL.cols());
  return LZWZL;
}

inline MatrixXd rts::hsgpCovariance::ZLu(const MatrixXd& u)
{
  MatrixXd ZL = rts::kronecker(ar_factor_chol, PhiSPD(true, false));
  return ZL * u;
}

inline MatrixXd rts::hsgpCovariance::Lu(const MatrixXd& u)
{
  MatrixXd ZL = rts::kronecker(ar_factor_chol, PhiSPD(true, false));
  return ZL * u;
}

inline sparse rts::hsgpCovariance::ZL_sparse()
{
  sparse dummy;
  return dummy;
}

inline int rts::hsgpCovariance::Q() const 
{
  return m * m * grid.T;
}

inline double rts::hsgpCovariance::log_likelihood(const VectorXd &u)
{
  double ll = 0;
  // need to collapse u to v
  MatrixXd umat(grid.N,grid.T);
  for(int t = 0; t< grid.T; t++){
    umat.col(t) = u.segment(t*grid.N,grid.N);
  }
  MatrixXd vmat = umat * ar_factor_inverse;
  double logdet = log_determinant();
  VectorXd uquad(grid.N);
  VectorXd vquad(grid.N);
  for(int t = 0; t< grid.T; t++){
    uquad = umat.col(t) * L;
    vquad = vmat.col(t) * L;
    ll += (-0.5*grid.N * log(2*M_PI) - 0.5*uquad.transpose()*vquad);
  }
  ll += 0.5*logdet;
  return -1.0*ll;
}

inline double rts::hsgpCovariance::log_determinant()
{
  double logdet = 0;
  for(int i = 0; i < (m*m); i++){
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

inline void rts::hsgpCovariance::update_rho(const double rho_)
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
    ar_factor_chol = MatrixXd(ar_factor.llt().matrixL());
    ar_factor_inverse = ar_factor.llt().solve(MatrixXd::Identity(grid.T,grid.T));
  } else {
    ar_factor_chol.setConstant(1.0);
    ar_factor_inverse.setConstant(1.0);
  }
}

inline void rts::hsgpCovariance::set_function(bool squared_exp)
{
  sq_exp = squared_exp;
}

inline void rts::hsgpCovariance::gen_indices()
{
  int counter = 0;
  for(int i = 1; i <= m; i++){
    for(int j = 1; j<= m; j++){
      indices(counter,0) = i;
      indices(counter,1) = j;
      counter++;
    }
  }
}

inline void rts::hsgpCovariance::gen_phi_prod()
{
  for(int i = 0; i < (m*m); i++){
    ArrayXd phi = phi_nD(i);
    Phi.col(i) = phi.matrix();
  }
  PhiT = Phi.transpose() * Phi;
}

inline MatrixXd rts::hsgpCovariance::PhiSPD(bool lambda, bool inverse)
{
  MatrixXd pnew = Phi;
  if(lambda){
    if(!inverse){
      pnew *= Lambda.sqrt().matrix().asDiagonal();
    } else {
      pnew *= Lambda.sqrt().inverse().matrix().asDiagonal();
    }
  }
  return pnew;
}

inline ArrayXd rts::hsgpCovariance::LambdaSPD()
{
  return Lambda;
}

inline void rts::hsgpCovariance::update_lambda()
{
  for(int i = 0; i < (m*m); i++) Lambda(i) = spd_nD(i);
  L = PhiSPD(true,true);
}
