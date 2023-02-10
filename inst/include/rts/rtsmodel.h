#ifndef RTSMODEL_H
#define RTSMODEL_H

#define _USE_MATH_DEFINES

#include <cmath> 
#include <glmmr.h>
#include <glmmrMCML.h>
#include <RcppEigen.h>
#include "eigenext.h"
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

namespace rts{

// a class for the region model LGCP
class rtsModel {
public:
  const Eigen::MatrixXd &X_;
  Eigen::VectorXd xb_;
  Eigen::MatrixXd zu_;
  Eigen::MatrixXd ZL_;
  Eigen::MatrixXd L_;
  Eigen::MatrixXd W_;
  Eigen::VectorXd y_;
  Eigen::MatrixXd u_;
  double var_par_;
  std::string family_; 
  std::string link_;
  Eigen::MatrixXd D_;
  Eigen::VectorXd offset_;
  int nT_;
  int nCell_;
  int n_;
  int Q_;
  int P_;
  int niter_;
  int flink;
  double rho_;
  bool useLflag_;
  
  rtsModel(
    const Eigen::MatrixXd &u,
    const Eigen::MatrixXd &X,
    Eigen::VectorXd y,
    Eigen::VectorXd beta,
    const Eigen::VectorXd &offset,
    int nT,
    double rho = 0.0
  ) : X_(X), 
  xb_(X.rows()),
  zu_(X.rows(),u.cols()),
  ZL_(Eigen::MatrixXd::Zero(1,1)),
  L_(Eigen::MatrixXd::Zero(1,1)),
  W_(Eigen::MatrixXd::Identity(X.rows(),X.rows())), y_(y),  
  u_(u), var_par_(1.0), 
  family_("poisson"), link_("log"),
  D_(u.rows(),u.rows()), offset_(offset),
  nT_(nT), rho_(rho), useLflag_(false) { 
    n_ = X.rows();
    nCell_ = n_/nT_;
    Q_ = u.rows();
    P_ = X.cols();
    update_beta(beta);
    niter_ = u.cols();
    flink = 1;
    update_u();
    update_W();
  }
  
  rtsModel(
    const Eigen::MatrixXd &L,
    int m,
    const Eigen::MatrixXd &X,
    Eigen::VectorXd y,
    Eigen::VectorXd beta,
    const Eigen::VectorXd &offset,
    int nT,
    double rho = 0.0
  ) : X_(X), 
  xb_(X.rows()),zu_(X.rows(),m), L_(L), 
  ZL_(Eigen::MatrixXd::Zero(X.rows(),X.rows())),
  W_(Eigen::MatrixXd::Identity(X.rows(),X.rows())), y_(y),  
  u_(Eigen::MatrixXd::Zero(X.rows(),m)), var_par_(1.0), 
  family_("poisson"), link_("log"), D_(L.rows(),L.rows()), 
  offset_(offset), nT_(nT), rho_(rho), useLflag_(true) { 
    D_ = L_ * L_.transpose();
    n_ = X.rows();
    nCell_ = n_/nT_;
    Q_ = n_;
    P_ = X.cols();
    niter_ = m;
    flink = 1;
    update_beta(beta);
    update_L();
    update_W();
  }
  
  
  void update_L(){
    D_ = L_ * L_.transpose();
    if(nT_==1){
      ZL_ = L_;
    } else {
      for(int t=0; t<nT_;t++){
        for(int s=t; s<nT_;s++){
          if(t==0){
            ZL_.block(s*nCell_,t*nCell_,nCell_,nCell_) = (pow(rho_,s)/(1-rho_*rho_))*L_;
          } else {
            ZL_.block(s*nCell_,t*nCell_,nCell_,nCell_) = pow(rho_,s-t)*L_;
          }
        }
        
      }
    }
    if(useLflag_){
      update_u();
    }
  }
  
  void update_L(const Eigen::MatrixXd &L){
    L_ = L;
    update_L();
  }
  
  void update_u(){
    if(useLflag_){
      if(nT_ == 1){
        zu_ = L_ * u_;
      } else {
        for(int t=0; t<nT_;t++){
          for(int s=t; s<nT_;s++){
            if(t==0){
              zu_.block(s*nCell_,0,nCell_,niter_) = (pow(rho_,s)/(1-rho_*rho_))*L_*u_.block(s*nCell_,0,nCell_,niter_);
            } else {
              zu_.block(s*nCell_,0,nCell_,niter_) += pow(rho_,s-t)*L_*u_.block(s*nCell_,0,nCell_,niter_);
            }
          }
        }
      }
    } else {
      if(nT_ == 1){
        zu_ = u_;
      } else {
        for(int t=0; t<nT_;t++){
          for(int s=t; s<nT_;s++){
            if(t==0){
              zu_.block(s*nCell_,0,nCell_,niter_) = (pow(rho_,s)/(1-rho_*rho_))*u_.block(s*nCell_,0,nCell_,niter_);
            } else {
              zu_.block(s*nCell_,0,nCell_,niter_) += pow(rho_,s-t)*u_.block(s*nCell_,0,nCell_,niter_);
            }
          }
        }
      }
    }
  }
  
  void update_u(const Eigen::MatrixXd &u){
    u_ = u;
    update_u();
  }
  
  void update_beta(const Eigen::VectorXd &beta){
    xb_ = X_* beta + offset_;
  }
  
  void update_rho(double rho){
    rho_ = rho;
    update_u();
  }
  
  void update_W(int i = 0){
    Eigen::VectorXd w = glmmr::maths::dhdmu(xb_ + zu_.col(i),family_,link_);
    for(int i = 0; i < n_; i++){
      W_(i,i) = 1/(w(i));
    }
  }
  
  // LOG GRADIENT
  Eigen::VectorXd log_grad(const Eigen::VectorXd &v){
    Eigen::VectorXd grad(v.size()); 
    Eigen::VectorXd mu = xb_;
    mu += ZL_*v;
    grad = -1.0*v;
    mu = (mu.array().exp()).matrix();
    grad.noalias() += ZL_.transpose()*(y_-mu);
    return grad;
  }
  
  Eigen::MatrixXd linpred(){
    Eigen::MatrixXd zd = zu_;
    for(int j=0; j<niter_ ; j++)zd.col(j) += xb_;
    return zd;
  }
  
  
  double log_likelihood() { 
    Eigen::ArrayXd ll = Eigen::ArrayXd::Zero(niter_);
#pragma omp parallel for
    for(int j=0; j<niter_ ; j++){
      for(int i = 0; i<n_; i++){
        ll(j) += glmmr::maths::log_likelihood(y_(i),xb_(i) + zu_.col(j)(i),var_par_,flink);
      }
    }
    return ll.mean();
  }
  
  
};

// a class for the region model LGCP
class rtsRegionModel {
public:
  const Eigen::MatrixXd &X_;
  Eigen::VectorXd xb_;
  Eigen::MatrixXd zu_;
  Eigen::MatrixXd ZL_;
  Eigen::MatrixXd L_;
  Eigen::MatrixXd W_;
  const Eigen::VectorXd &y_;
  Eigen::MatrixXd u_;
  Eigen::MatrixXd regionu_;
  double var_par_;
  std::string family_; 
  std::string link_;
  Eigen::MatrixXd D_;
  Eigen::VectorXd offset_;
  const Eigen::ArrayXi &n_cell_;
  const Eigen::ArrayXi &cell_id_;
  const Eigen::ArrayXd &q_weights_;
  int nT_;
  int nCell_;
  int nRegion_;
  int n_;
  int Q_;
  int P_;
  int niter_;
  int flink;
  double rho_;
  bool useLflag_;
  
  rtsRegionModel(
    const Eigen::MatrixXd &u,
    const Eigen::MatrixXd &X,
    const Eigen::VectorXd &y,
    Eigen::VectorXd beta,
    const Eigen::VectorXd &offset,
    int nT,
    const Eigen::ArrayXi &n_cell,
    const Eigen::ArrayXi &cell_id,
    const Eigen::ArrayXd &q_weights,
    double rho = 0.0
  ) : X_(X), 
  xb_(X.rows()),
  zu_(Eigen::MatrixXd::Zero(u.rows(),u.cols())),
  ZL_(Eigen::MatrixXd::Identity(1,1)),
  L_(Eigen::MatrixXd::Identity(1,1)),
  W_(Eigen::MatrixXd::Identity(X.rows(),X.rows())), 
  y_(y),  
  u_(u), 
  regionu_(Eigen::MatrixXd::Zero(y.size(),u.cols())), 
  var_par_(1.0), family_("poisson"), link_("log"), 
  D_(u.rows(),u.rows()), offset_(offset),
  n_cell_(n_cell), cell_id_(cell_id),
  q_weights_(q_weights), nT_(nT), rho_(rho), useLflag_(false) { 
    n_ = y.size();
    nRegion_ = n_cell.size()-1;
    nCell_ = u.rows()/nT_;
    Q_ = nCell_*nT_;
    P_ = X.cols();
    niter_ = u.cols();
    flink = 1;
    update_beta(beta);
    update_u();
    update_W();
  }
  
  rtsRegionModel(
    const Eigen::MatrixXd &L,
    int m,
    const Eigen::MatrixXd &X,
    const Eigen::VectorXd &y,
    Eigen::VectorXd beta,
    const Eigen::VectorXd &offset,
    int nT,
    const Eigen::ArrayXi &n_cell,
    const Eigen::ArrayXi &cell_id,
    const Eigen::ArrayXd &q_weights,
    double rho = 0.0
  ) : X_(X), 
  xb_(X.rows()),
  zu_(Eigen::MatrixXd::Zero(L.rows()*nT,m)), 
  L_(L), 
  ZL_(Eigen::MatrixXd::Zero(X.rows(),X.rows())),
  W_(Eigen::MatrixXd::Identity(X.rows(),X.rows())), y_(y),  
  u_(Eigen::MatrixXd::Zero(L.cols()*nT,m)), 
  regionu_(Eigen::MatrixXd::Zero(y.size(),m)), 
  D_(L.rows(),L.rows()), var_par_(1.0), 
  family_("poisson"), link_("log"), 
  offset_(offset), n_cell_(n_cell), cell_id_(cell_id),
  q_weights_(q_weights), nT_(nT), rho_(rho), useLflag_(true) { 
    n_ = y.size(); // total number of observations
    nRegion_ = n_cell.size()-1; //number of regions
    nCell_ = L.rows(); //number of cells
    Q_ = nCell_*nT_; //number of random effects
    P_ = X.cols();
    niter_ = m;
    flink = 1;
    var_par_=1.0;
    update_beta(beta);
    update_W();
    update_L();
  }
  
  
  void update_L(){
    D_ = L_*L_.transpose();
    if(nT_==1){
      ZL_ = L_;
    } else {
      for(int t=0; t<nT_;t++){
        for(int s=t; s<nT_;s++){
          if(t==0){
            ZL_.block(s*nCell_,t*nCell_,nCell_,nCell_) = (pow(rho_,s)/(1-rho_*rho_))*L_;
          } else {
            ZL_.block(s*nCell_,t*nCell_,nCell_,nCell_) = pow(rho_,s-t)*L_;
          }
        }
        
      }
    }
    if(useLflag_){
      update_u();
    }
  }
  
  void update_L(const Eigen::MatrixXd &L){
    L_ = L;
    update_L();
  }
  
  
  Eigen::MatrixXd get_ZL(){
    return ZL_;
  }
  
  void update_u(){
    if(useLflag_){
      if(nT_ == 1){
        zu_ = L_ * u_;
      } else {
        for(int t=0; t<nT_;t++){
          for(int s=t; s<nT_;s++){
            if(t==0){
              zu_.block(s*nCell_,0,nCell_,niter_) = (pow(rho_,s)/(1-rho_*rho_))*L_*u_.block(s*nCell_,0,nCell_,niter_);
            } else {
              zu_.block(s*nCell_,0,nCell_,niter_) += pow(rho_,s-t)*L_*u_.block(s*nCell_,0,nCell_,niter_);
            }
          }
        }
      }
    } else {
      if(nT_ == 1){
        zu_ = u_;
      } else {
        for(int t=0; t<nT_;t++){
          for(int s=t; s<nT_;s++){
            if(t==0){
              zu_.block(s*nCell_,0,nCell_,niter_) = (pow(rho_,s)/(1-rho_*rho_))*u_.block(s*nCell_,0,nCell_,niter_);
            } else {
              zu_.block(s*nCell_,0,nCell_,niter_) += pow(rho_,s-t)*u_.block(s*nCell_,0,nCell_,niter_);
            }
          }
        }
      }
    }
    
#pragma omp parallel for
    for(int r=0; r<nRegion_;r++){
      for(int t=0; t<nT_; t++){
        int nInter = n_cell_(r+1)-n_cell_(r);
        for(int j=0; j<niter_; j++){
          double accum = 0;
          for(int l=0; l<nInter; l++){
            accum += q_weights_(n_cell_(r)-1+l)*exp(zu_(cell_id_(n_cell_(r)-1+l) + t*nCell_,j));
          }
          regionu_(r + t*nRegion_,j) = accum;
        }
      }
    }
  }
  
  void update_u(const Eigen::MatrixXd &u){
    u_ = u;
    update_u();
  }
  
  void update_beta(const Eigen::VectorXd &beta){
    xb_ = X_* beta + offset_;
  }
  
  void update_rho(double rho){
    rho_ = rho;
    update_u();
  }
  
  
  Eigen::ArrayXXd region_intensity(bool uselog = true){
    Eigen::ArrayXXd intens = Eigen::ArrayXXd::Zero(n_,niter_);
    Eigen::ArrayXd expxb = xb_.array().exp();
    for(int j=0; j<niter_; j++){
      intens.col(j) = expxb * regionu_.col(j).array();
    }
    if(uselog){
      return intens.log();
    } else {
      return intens;
    }
  }
  
  // LOG GRADIENT
  Eigen::VectorXd log_grad(const Eigen::VectorXd &v){
    Eigen::VectorXd grad(v.size()); 
    Eigen::VectorXd mu = xb_;
    mu += ZL_*v;
    grad = -1.0*v;
    mu = (mu.array().exp()).matrix();
    grad.noalias() += ZL_.transpose()*(y_-mu);
    return grad;
  }
  
  void update_W(int i = 0){
    Eigen::ArrayXXd rintens(n_,niter_);
    rintens = region_intensity(true);
    Eigen::VectorXd w = glmmr::maths::dhdmu(rintens.col(i).matrix(),family_,link_);
    for(int i = 0; i < n_; i++){
      W_(i,i) = 1/(w(i));
    }
  }
  
  Eigen::MatrixXd linpred(){
    Eigen::ArrayXXd rintens(n_,niter_);
    rintens = region_intensity(true);
    return rintens.matrix();
  }
  
  
  double log_likelihood(bool useL = false) { 
    Eigen::ArrayXd ll = Eigen::ArrayXd::Zero(niter_);
    Eigen::MatrixXd rintens = linpred();
#pragma omp parallel for
    for(int j=0; j<niter_ ; j++){
      double accum = 0;
      for(int i = 0; i<n_; i++){
        accum += glmmr::maths::log_likelihood(y_(i),rintens(i,j),1.0,flink);
      }
      ll(j) = accum;
    }
    return ll.mean();
  }
  
  
};

}

#endif