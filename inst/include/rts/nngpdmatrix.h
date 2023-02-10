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
#define omp_get_thread_num() 0
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

namespace rts {
class NNGPDmatrix {
public: 
  glmmr::DData* data_;
  Eigen::ArrayXd gamma_;
  Eigen::ArrayXXi NN_;
  int m_; // number of nearest neighbours
  int n_;
  int nT_;
  Eigen::MatrixXd A_;
  Eigen::VectorXd D_;
  
  NNGPDmatrix(
    glmmr::DData* data,
    const Eigen::ArrayXXi &NN,
    const Eigen::ArrayXd& gamma,
    int nT
  ) : data_(data), NN_(NN), gamma_(gamma), m_(NN.rows()), n_(NN.cols()),
  nT_(nT), A_(m_,n_), D_(n_) {}
  
  double loglik_col(const Eigen::MatrixXd &u){
    double logdet = 0.0;
    double qf = 0.0;
    // Eigen::ArrayXd logdet = Eigen::ArrayXd::Zero(u.cols()); 
    // Eigen::ArrayXd qf = Eigen::ArrayXd::Zero(u.cols());
    double ll = 0;
    genAD();
    
#pragma omp parallel for 
    for(int i = 0; i < n_; i++){
      for(int t = 0; t<nT_; t++){
        int idxlim = i <= m_ ? i - 1 : m_;
        double au;
        Eigen::ArrayXi ja(1);
        ja(0) = 0;
        Eigen::MatrixXd usec = glmmr::Eigen_ext::mat_indexing(u, NN_.col(i).segment(0,idxlim), ja);
        logdet += log(D_(i));
        if(i == 0){
          qf += u(t*n_,0)*u(t*n_,0)/D_(0);
        } else {
          au = u(t*n_ + i,0) - (A_.col(i).segment(0,idxlim).transpose() * usec)(0);
          qf += au*au/D_(i);
        }
      }
    }
    
    ll = -0.5*logdet -0.5*qf - 0.5*n_*3.141593;
    
    return ll;
  }
  
  double loglik(const Eigen::MatrixXd &u){
    double ll = 0;
    int niter = u.cols();
    for(int j=0; j<niter; j++){
      ll += loglik_col(u.col(j));
    }
    return ll/niter;
  }
  
  void update_parameters(const Eigen::ArrayXd& gamma){
    gamma_ = gamma;
  }
  
  
  void update_parameters(const Eigen::VectorXd& gamma){
    gamma_ = gamma.array();
  }
  
  
  void update_parameters(const std::vector<double>& gamma){
    std::vector<double> par = gamma;
    Eigen::ArrayXd gammaa = Eigen::Map<Eigen::ArrayXd>(par.data(),par.size());
    gamma_ = gammaa;
  }
  
 
  
  void genAD(){
    A_ = Eigen::MatrixXd::Zero(m_,n_);
    D_ = Eigen::VectorXd::Zero(n_);
    int idxlim;
    Eigen::MatrixXd S(m_,m_);
    Eigen::VectorXd Sv(m_);
    D_(0) = gamma_(0);
    
    glmmr::DSubMatrix* dblock;
    dblock = new glmmr::DSubMatrix(0, data_, gamma_);
    double val = dblock->get_val(0,0);
    
    for(int i = 1; i < n_; i++){
      idxlim = i<= m_ ? i-1 : m_;
      S = Eigen::MatrixXd::Zero(m_,m_);
      Sv = Eigen::VectorXd::Zero(m_);
      for(int j = 0; j<idxlim; j++){
        S(j,j) = gamma_(0);
      }
      if(idxlim > 1){
        for(int j = 0; j<(idxlim-1); j++){
          for(int k = j+1; k<idxlim; k++){
            S(j,k) = dblock->get_val(NN_(j,i),NN_(k,i));
            S(k,j) = S(j,k);
          }
        }
      }
      for(int j = 0; j<idxlim; j++){
        Sv(j) = dblock->get_val(i,NN_(j,i));
      }
      A_.block(0,i,idxlim,1) = S.block(0,0,idxlim,idxlim).ldlt().solve(Sv.segment(0,idxlim));
      D_(i) = gamma_(0) - (A_.col(i).segment(0,idxlim).transpose() * Sv.segment(0,idxlim))(0);
    }
    delete dblock;
  }
  
  Eigen::MatrixXd chol(){
    genAD();
    return rts::inv_ldlt_AD(A_,D_,NN_);
  }
  
};

class rtsDMatrix {
public: 
  glmmr::DData* data_;
  Eigen::VectorXd gamma_;
  glmmr::MCMLDmatrix dmat_;
  int nT_;
  int n_;
  
  rtsDMatrix(
    glmmr::DData* data,
    const Eigen::ArrayXd& gamma,
    int nT
  ) : data_(data), gamma_(gamma), dmat_(data,gamma), nT_(nT) {
    dmat_.data_->subdata(0);
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
    int n = dmat_.data_->N();
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
    Eigen::MatrixXd L = dmat_.genD(0,true,false);
    return L;
  }
  
};

}

#endif