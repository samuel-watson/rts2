#ifndef NNGPDMATRIX_H
#define NNGPDMATRIX_H

#define _USE_MATH_DEFINES

#include <cmath> 
#include <glmmr.h>
#include <glmmrMCML.h>
#include <RcppEigen.h>
#include "eigenext.h"
#include "rtsmaths.h"
#ifdef _OPENMP
#include <omp.h>     
#else
// for machines with compilers void of openmp support
#define omp_get_num_threads()  1
#define omp_get_thread_num()   0
#define omp_get_max_threads()  1
#define omp_get_thread_limit() 1
#define omp_get_num_procs()    1
#define omp_set_nested(a)   // empty statement to remove the call
#define omp_get_wtime()        0
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

namespace rts {
class NNGPDmatrix {
public: 
  glmmr::DData data_;
  Eigen::ArrayXd gamma_;
  glmmr::DSubMatrix dblock_;
  Eigen::ArrayXXi NN_;
  int m_; // number of nearest neighbours
  int n_;
  int nT_;
  Eigen::MatrixXd A_;
  Eigen::VectorXd D_;
  
  NNGPDmatrix(
    const Eigen::ArrayXXi &cov,
    const Eigen::ArrayXd &data,
    const Eigen::ArrayXd &eff_range,
    const Eigen::ArrayXXi &NN,
    const Eigen::ArrayXd& gamma,
    int nT
  ) : data_(cov,data,eff_range),gamma_(gamma), dblock_(data_,gamma_),
  NN_(NN),  m_(NN.rows()), n_(NN.cols()),
  nT_(nT), A_(m_,n_), D_(n_), S(m_,m_), Sv(m_),
  ll1(0.0), ja(1){
    ja(0) = 0;
    genAD();
  }
  
  double loglik(const Eigen::MatrixXd &u){
    genAD();
    ll1 = 0.0;
    double qf = 0.0;
    double logdet = 0.0;
    
    for(int j=0; j<u.cols(); j++){
      for(int t = 0; t<nT_; t++){
        int idxlim;
        for(int i = 0; i < n_; i++){
          qf = 0.0;
          logdet = 0.0;
          idxlim = i <= m_ ? i - 1 : m_;
          Eigen::MatrixXd usec = glmmr::Eigen_ext::mat_indexing(u, NN_.col(i).segment(0,idxlim), ja);
          logdet += log(D_(i));//logdet
          if(i == 0){
            qf += u(t*n_,0)*u(t*n_,0)/D_(0); //qf
          } else {
            double au = u(t*n_ + i,0) - (A_.col(i).segment(0,idxlim).transpose() * usec)(0);
            qf += au*au/D_(i);//qf
          }
        }
        ll1 += -0.5*qf - 0.5*logdet - 0.5*n_*3.141593; 
      }
    }
    
    return ll1/u.cols();
  }
  
  void update_parameters(const Eigen::ArrayXd& gamma){
    gamma_ = gamma;
    dblock_.gamma_ = gamma.matrix();
  }
  
  
  void update_parameters(const Eigen::VectorXd& gamma){
    gamma_ = gamma.array();
    dblock_.gamma_ = gamma;
  }
  
  
  void update_parameters(const std::vector<double>& gamma){
    std::vector<double> par = gamma;
    Eigen::ArrayXd gammaa = Eigen::Map<Eigen::ArrayXd>(par.data(),par.size());
    gamma_ = gammaa;
    dblock_.gamma_ = gammaa.matrix();
  }
  
 
  
  void genAD(){
    A_.setZero();
    D_.setZero();
    S.setZero();
    Sv.setZero();
    int idxlim;
    double val = dblock_.get_val(0,0);
    D_(0) = val;
    
    for(int i = 1; i < n_; i++){
      idxlim = i <= m_ ? i : m_;
      S.setZero();
      Sv.setZero();
      for(int j = 0; j<idxlim; j++){
        S(j,j) = val;
      }
      if(idxlim > 1){
        for(int j = 0; j<(idxlim-1); j++){
          for(int k = j+1; k<idxlim; k++){
            S(j,k) = dblock_.get_val(NN_(j,i),NN_(k,i));
            S(k,j) = S(j,k);
          }
        }
      }
      for(int j = 0; j<idxlim; j++){
        Sv(j) = dblock_.get_val(i,NN_(j,i));
      }
      A_.block(0,i,idxlim,1) = S.block(0,0,idxlim,idxlim).ldlt().solve(Sv.segment(0,idxlim));
      D_(i) = val - (A_.col(i).segment(0,idxlim).transpose() * Sv.segment(0,idxlim))(0);
    }
  }
  
  Eigen::MatrixXd chol(){
    genAD();
    return rts::inv_ldlt_AD(A_,D_,NN_);
  }
private:
  Eigen::MatrixXd S;
  Eigen::VectorXd Sv;
  double ll1;
  Eigen::ArrayXi ja;
  
};

class rtsDMatrix {
public: 
  glmmr::DData data_;
  Eigen::VectorXd gamma_;
  glmmr::MCMLDmatrix dmat_;
  int nT_;
  int n_;
  
  rtsDMatrix(
    const Eigen::ArrayXXi &cov,
    const Eigen::ArrayXd &data,
    const Eigen::ArrayXd &eff_range,
    const Eigen::ArrayXd& gamma,
    int nT
  ) : data_(cov,data,eff_range), gamma_(gamma), 
  dmat_(cov,data,eff_range,gamma_), nT_(nT) {
    dmat_.data_.subdata(0);
  }
  
  void update_parameters(const Eigen::ArrayXd& gamma){
    gamma_ = gamma;
    dmat_.update_parameters(gamma);
  }
  
  void update_parameters(const Eigen::VectorXd& gamma){
    gamma_ = gamma.array();
    dmat_.update_parameters(gamma);
  }
  
  void update_parameters(const std::vector<double>& gamma){
    dmat_.update_parameters(gamma);
    std::vector<double> par = gamma;
    Eigen::ArrayXd gammaa = Eigen::Map<Eigen::ArrayXd>(par.data(),par.size());
    gamma_ = gammaa;
  }
  
  // log like + chol functions
  double loglik(const Eigen::MatrixXd &u){
    double ll = 0;
    int n = dmat_.data_.N();
    if(nT_==1){
      ll = dmat_.loglik(u);
    } else {
      Eigen::MatrixXd usub = u.block(0,0,n,u.cols());
      ll += dmat_.loglik(usub);
      for(int t=1; t<nT_; t++){
        usub = u.block(t*n,0,n,u.cols());
        ll += dmat_.loglik(usub);
      }
    }
    return ll;
  }
  
  Eigen::MatrixXd chol(){
    dmat_.update_parameters(gamma_);
    Eigen::MatrixXd L = dmat_.genD(true,false);
    return L;
  }
  
};

}

#endif