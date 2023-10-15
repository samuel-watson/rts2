# pragma once

#include <glmmr/modeloptim.hpp>
#include "rtsmodelbits.h"
#include "griddata.h"
#include "regiondata.h"
#include "regionlinearpredictor.h"

namespace rts {

using namespace rminqa;
using namespace Eigen;
using namespace glmmr;

template<typename modeltype>
class rtsRegionModelOptim : public ModelOptim<modeltype> {
public:
  rts::RegionData& region;
  
  rtsRegionModelOptim(modeltype& model_, 
                glmmr::ModelMatrix<modeltype>& matrix_,
                glmmr::RandomEffects<modeltype>& re_,
                rts::RegionData& region_) : ModelOptim<modeltype>(model_,matrix_,re_), region(region_) {};
  
  rtsRegionModelOptim(const rts::rtsRegionModelOptim<modeltype>& optim) : ModelOptim<modeltype>(optim.model, optim.matrix, optim.re), region(optim.region) {};
  
  void update_theta(const dblvec &theta) override;
  void update_u(const MatrixXd& u) override;
  void update_rho(double rho);
  void ml_rho();
  void ml_theta() override;
  double log_likelihood() override;
  double full_log_likelihood() override;
  ArrayXXd region_intensity(bool uselog = true);
  ArrayXXd y_predicted(bool uselog = true);
  MatrixXd hessian(double tol = 1e-4);
  void set_bobyqa_control(int npt_, double rhobeg_, double rhoend_);
  
  private:
    int npt = 0;
    double rhobeg = 0.1;
    double rhoend = 1e-6;
    
    class rho_likelihood : public Functor<dblvec> {
      rtsRegionModelOptim<modeltype>& M_;
      double parrho;
      double ll;
    public:
      rho_likelihood(rtsRegionModelOptim<modeltype>& M) :  
      M_(M), parrho(0.0), ll(0.0) {};
      double operator()(const dblvec &par);
    };
    
    class D_likelihood_hsgp : public Functor<dblvec> {
      rtsRegionModelOptim<modeltype>& M;
      double logl;
    public:
      D_likelihood_hsgp(rtsRegionModelOptim<modeltype>& M_) :
      M(M_),
      logl(0.0) {};
      double operator()(const dblvec &par);
    };
  
};

}

template<typename modeltype>
inline void rts::rtsRegionModelOptim<modeltype>::set_bobyqa_control(int npt_, double rhobeg_, double rhoend_){
  npt = npt_;
  rhobeg = rhobeg_;
  rhoend = rhoend_;
}

template<typename modeltype>
inline void rts::rtsRegionModelOptim<modeltype>::update_theta(const dblvec &theta){
  this->model.covariance.update_parameters(theta);
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
}

template<typename modeltype>
inline void rts::rtsRegionModelOptim<modeltype>::update_rho(double rho){
  this->model.covariance.update_rho(rho);
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
}

template<typename modeltype>
inline void rts::rtsRegionModelOptim<modeltype>::ml_rho(){
  rho_likelihood ldl(*this);
  Rbobyqa<rho_likelihood,dblvec> opt;
  dblvec start;
  start.push_back(this->model.covariance.rho);
  dblvec lower;
  lower.push_back(-1.0);
  dblvec upper;
  upper.push_back(1.0);
  opt.set_lower(lower);
  opt.set_upper(upper);
  opt.control.iprint = this->trace;
  opt.minimize(ldl, start);
}

template<typename modeltype>
inline void rts::rtsRegionModelOptim<modeltype>::update_u(const MatrixXd& u_){
  if(u_.cols()!=this->re.u(false).cols()){
    this->re.u_.conservativeResize(this->model.covariance.Q(),u_.cols());
    this->re.zu_.resize(this->model.covariance.Q(),u_.cols());
  }
  this->re.u_ = u_;
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
}

template<typename modeltype>
inline double rts::rtsRegionModelOptim<modeltype>::log_likelihood(){
  double ll = 0;
  ArrayXXd xb = y_predicted(true);
  if(this->model.weighted){
#pragma omp parallel for reduction (+:ll) collapse(2)
    for(int j=0; j<xb.cols() ; j++){
      for(int i = 0; i<xb.rows(); i++){
        ll += this->model.data.weights(i)*glmmr::maths::log_likelihood(this->model.data.y(i),xb(i,j),this->model.data.variance(i),this->model.family.flink);
      }
    }
    ll *= this->model.data.weights.sum()/this->model.n();
  } else {
#pragma omp parallel for reduction (+:ll) collapse(2)
  for(int j=0; j<xb.cols() ; j++){
    for(int i = 0; i< xb.rows(); i++){
      ll += glmmr::maths::log_likelihood(this->model.data.y(i),xb(i,j),this->model.data.variance(i),this->model.family.flink);
    }
  }
}
  return ll/xb.cols();
}

template<typename modeltype>
inline ArrayXXd rts::rtsRegionModelOptim<modeltype>::y_predicted(bool uselog){
  ArrayXXd xb(this->model.n(), this->re.u_.cols());
  if constexpr (std::is_same_v<modeltype, BitsAR> || std::is_same_v<modeltype, BitsNNGP > || std::is_same_v<modeltype, BitsHSGP >){
    xb = region_intensity();
  } else if constexpr (std::is_same_v<modeltype, BitsARRegion > || std::is_same_v<modeltype, BitsNNGPRegion > || std::is_same_v<modeltype, BitsHSGPRegion >){
    xb = this->model.linear_predictor.xb_region(this->re.zu_);
  }
  xb.matrix().colwise() += this->model.data.offset;
  if(!uselog)xb = xb.exp();
  return xb;
}

template<typename modeltype>
inline ArrayXXd rts::rtsRegionModelOptim<modeltype>::region_intensity(bool uselog){
  MatrixXd regionu = region.grid_to_region(this->re.zu_);
  ArrayXXd intens = ArrayXXd::Zero(region.nRegion,this->re.u_.cols());
  ArrayXd expxb = this->model.linear_predictor.xb().array().exp();
  for(int j=0; j<intens.cols(); j++){
    intens.col(j) = expxb * regionu.col(j).array();
  }
  if(uselog){
    return intens.log();
  } else {
    return intens;
  }
}

template<typename modeltype>
inline double rts::rtsRegionModelOptim<modeltype>::full_log_likelihood(){
  double ll = rts::rtsRegionModelOptim<modeltype>::log_likelihood();
  double logl = 0;
  MatrixXd Lu = this->model.covariance.Lu(this->re.u_);
#pragma omp parallel for reduction (+:logl)
  for(int i = 0; i < Lu.cols(); i++){
    logl += this->model.covariance.log_likelihood(Lu.col(i));
  }
  logl *= 1/Lu.cols();
  return ll+logl;
}

template<typename modeltype>
inline double rts::rtsRegionModelOptim<modeltype>::rho_likelihood::operator()(const dblvec &par){
  parrho = par[0];
  M_.update_rho(parrho);
  ll = M_.log_likelihood();
  return -1*ll;
}

// template<typename modeltype>
// inline void rts::rtsRegionModelOptim<modeltype>::ml_theta(){
//   glmmr::ModelOptim<modeltype>::ml_theta();
// }

//BitsHSGP
template<typename modeltype>
inline void rts::rtsRegionModelOptim<modeltype>::ml_theta(){
  D_likelihood_hsgp ddl(*this);
  Rbobyqa<D_likelihood_hsgp,dblvec> opt;
  dblvec lower = this->get_lower_values(false,true,false);
  opt.set_lower(lower);
  opt.control.iprint = this->trace;
  opt.control.npt = npt;
  opt.control.rhobeg = rhobeg;
  opt.control.rhoend = rhoend;
  dblvec start_t = this->get_start_values(false,true,false);
  opt.minimize(ddl, start_t);
}

// template<>
// inline void rts::rtsRegionModelOptim<BitsHSGPRegion>::ml_theta(){
//   D_likelihood_hsgp ddl(*this);
//   Rbobyqa<D_likelihood_hsgp,dblvec> opt;
//   dblvec lower = this->get_lower_values(false,true,false);
//   opt.set_lower(lower);
//   opt.control.iprint = trace;
//   dblvec start_t = this->get_start_values(false,true,false);
//   opt.minimize(ddl, start_t);
// }

template<typename modeltype>
inline double rts::rtsRegionModelOptim<modeltype>::D_likelihood_hsgp::operator()(const dblvec &par) {
  M.update_theta(par);
  logl = M.log_likelihood();
  return -1*logl;
}

template<typename modeltype>
inline MatrixXd rts::rtsRegionModelOptim<modeltype>::hessian(double tol){
  typename ModelOptim<modeltype>::L_likelihood ldl(*this);
  dblvec hess(this->model.linear_predictor.P()*this->model.linear_predictor.P());
  std::fill(hess.begin(),hess.end(),0.0);
  dblvec ndeps(this->model.linear_predictor.P());
  std::fill(ndeps.begin(),ndeps.end(),tol);
  ldl.os.ndeps_ = ndeps;
  ldl.Hessian(this->model.linear_predictor.parameters,hess);
  MatrixXd H = Map<MatrixXd>(hess.data(),this->model.linear_predictor.P(),this->model.linear_predictor.P());
  return H;
}