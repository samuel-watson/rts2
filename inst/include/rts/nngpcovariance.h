#pragma once

#define _USE_MATH_DEFINES

#include <glmmr/covariance.hpp>
#include "griddata.h"
#include "rtsmaths.h"

namespace rts {

using namespace Eigen;
using namespace glmmr;

class nngpCovariance : public Covariance {
public:
  double rho;
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
                   A(m_,data.rows()), Dvec(data.rows()), m(m_), ar_factor(T,T), ar_factor_chol(T,T) {
    isSparse = false;
    grid.genNN(m);
    gen_AD();
    update_rho(0.1);
  }
  
  nngpCovariance(const str& formula,
                 const ArrayXXd &data,
                 const strvec& colnames,
                 int T, int m_,
                 const rts::griddata& grid_) : Covariance(formula, data, colnames),  
                 grid(grid_),
                 A(m_,data.rows()), Dvec(data.rows()), m(m_), ar_factor(T,T), ar_factor_chol(T,T) {
    isSparse = false;
    grid.genNN(m);
    update_rho(0.1);
  }
  
  nngpCovariance(const rts::nngpCovariance& cov) : Covariance(cov.form_, cov.data_, cov.colnames_, cov.parameters_),  
    grid(cov.grid), A(grid.m,grid.N), Dvec(grid.N), m(cov.m), ar_factor(grid.T,grid.T), ar_factor_chol(grid.T,grid.T) {
      isSparse = false;
      grid.genNN(m);
      gen_AD();
      update_rho(cov.rho);
    }
  
  MatrixXd inv_ldlt_AD(const MatrixXd &A, const VectorXd &D, const ArrayXXi &NN);
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
  
protected:
  MatrixXd ar_factor;    
  MatrixXd ar_factor_chol;
  // VectorXd tri_mult(const MatrixXd& M, const VectorXd& V);
};

}

// inline VectorXd rts::nngpCovariance::tri_mult(const MatrixXd& M, const VectorXd& V){
//   VectorXd MV = VectorXd::Zero(M.rows());
// #pragma omp parallel for
//   for(int i = 0; i < M.rows(); i++){
//     for(int j = 0; j < (i+1); j++){
//       MV(i) += M(i,j) * V(j);
//     }
//   }
//   return MV;
// }

inline MatrixXd rts::nngpCovariance::inv_ldlt_AD(const MatrixXd &A, 
                            const VectorXd &D,
                            const ArrayXXi &NN){
  int n = A.cols();
  int m = A.rows();
  MatrixXd y = MatrixXd::Zero(n,n);
  ArrayXd dsqrt = Dvec.array().sqrt();
#pragma omp parallel for  
  for(int k=0; k<n; k++){
    int idxlim;
    for (int i = 0; i < n; i++) {
      idxlim = i<=m ? i : m;
      double lsum = 0;
      for (int j = 0; j < idxlim; j++) {
        lsum += A(j,i) * y(NN(j,i),k);
      }
      y(i,k) = i==k ? (1+lsum)*dsqrt(k)  : lsum*dsqrt(k) ;
    }
  }
  
  return y;
}

inline MatrixXd rts::nngpCovariance::D(bool chol, bool upper){
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

inline MatrixXd rts::nngpCovariance::ZL(){
  MatrixXd L = inv_ldlt_AD(A,Dvec,grid.NN);//glmmr::sparse_to_dense(this->matL,false);
  MatrixXd ZL = rts::kronecker(ar_factor_chol, L);
  return ZL;
}

inline MatrixXd rts::nngpCovariance::LZWZL(const VectorXd& w){
  MatrixXd ZL = rts::nngpCovariance::ZL();
  MatrixXd LZWZL = ZL.transpose() * w.asDiagonal() * ZL;
  LZWZL += MatrixXd::Identity(LZWZL.rows(), LZWZL.cols());
  return LZWZL;
}

inline MatrixXd rts::nngpCovariance::ZLu(const MatrixXd& u){
  MatrixXd L = inv_ldlt_AD(A,Dvec,grid.NN);
  MatrixXd ZLU(grid.T*grid.N,u.cols());
  if(grid.T == 1){
    for(int i = 0; i < u.cols(); i++){
      ZLU = L * u;
      //ZLU.col(i) = tri_mult(L,u.col(i)); 
    }
  } else {
    MatrixXd umat(grid.N,grid.T);
    for(int i = 0; i < u.cols(); i++){
      for(int t = 0; t< grid.T; t++){
        umat.col(t) = L * u.segment(t*grid.N,grid.N);
        //umat.col(t) = tri_mult(L,u.col(i).segment(t*grid.N,grid.N));
      }
      umat = umat * ar_factor_chol.transpose();
      for(int t = 0; t< grid.T; t++){
        ZLU.block(t*grid.N,i,grid.N,1) = umat.col(t);
      }
    }
  }
  //MatrixXd ZLu = rts::nngpCovariance::ZL() * u;
  
  return ZLU;
}

inline MatrixXd rts::nngpCovariance::Lu(const MatrixXd& u){
  return rts::nngpCovariance::ZLu(u);
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
  int idxlim;
  double au,av;
  MatrixXd ar_factor_inverse = ar_factor.llt().solve(MatrixXd::Identity(grid.T,grid.T));
  // need to collapse u to v
  MatrixXd umat(grid.N,grid.T);
  for(int t = 0; t< grid.T; t++){
    umat.col(t) = u.segment(t*grid.N,grid.N);
  }
  MatrixXd vmat = umat * ar_factor_inverse;
  
  for(int t = 0; t<grid.T; t++){
    double qf = umat(0,t)*vmat(0,t)/Dvec(0);
#pragma omp parallel for reduction (+:qf)
    for(int i = 1; i < grid.N; i++){
      idxlim = i <= m ? i : m;
      VectorXd usec(idxlim);
      VectorXd vsec(idxlim);
      for(int j = 0; j < idxlim; j++) {
        usec(j) = umat(grid.NN(j,i),t);
        vsec(j) = vmat(grid.NN(j,i),t);
      }
      au = umat(i,t) - (A.col(i).segment(0,idxlim).transpose() * usec)(0);
      av = vmat(i,t) - (A.col(i).segment(0,idxlim).transpose() * vsec)(0);
      qf += au*av/Dvec(i);
    }
    ll1 -= 0.5*qf + 0.5*grid.N*log(2*M_PI); 
  }
  ll1 -= 0.5*logdet;
  return ll1;
}

inline double rts::nngpCovariance::log_determinant(){
  double logdet = Dvec.array().log().sum();
  logdet *= grid.T;
  double logdet_ar = 0;
  if(grid.T > 1){
    for(int t = 0; t < grid.T; t++) logdet_ar += 2*log(ar_factor_chol(t,t));
    logdet_ar *= grid.N;
  }
  
  return logdet + logdet_ar;
}

inline void rts::nngpCovariance::update_rho(const double rho_){
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

inline void rts::nngpCovariance::gen_AD(){
  A.setZero();
  Dvec.setZero();
  
  int idxlim;
  double val = Covariance::get_val(0,0,0);
  Dvec(0) = val;
  
#pragma omp parallel for
  for(int i = 1; i < grid.N; i++){
    idxlim = i <= m ? i : m;
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
    A.block(0,i,idxlim,1) = S.ldlt().solve(Sv);
    Dvec(i) = val - (A.col(i).segment(0,idxlim).transpose() * Sv)(0);
  }
  //this->matL = glmmr::dense_to_sparse(D(true,false),false);
}

inline vector_matrix rts::nngpCovariance::submatrix(int i){
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


