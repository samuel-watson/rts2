# pragma once

#include <glmmr/modeloptim.hpp>
#include "rtsmodelbits.h"
#include "regionlinearpredictor.h"

namespace rts {

using namespace rminqa;
using namespace Eigen;
using namespace glmmr;

template<typename modeltype>
class rtsModelOptim : public ModelOptim<modeltype> {  public:
    
    rtsModelOptim(modeltype& model_, 
               glmmr::ModelMatrix<modeltype>& matrix_,
               glmmr::RandomEffects<modeltype>& re_) : ModelOptim<modeltype>(model_,matrix_,re_) {};
    
    rtsModelOptim(const rts::rtsModelOptim<modeltype>& optim) : ModelOptim<modeltype>(optim.model,optim.matrix,optim.re) {};
    
    void update_theta(const dblvec &theta) override;
    void update_u(const MatrixXd& u) override;
    void update_rho(const double rho_);
    void ml_rho();
    void ml_theta() override;
    
  private:
    class rho_likelihood : public Functor<dblvec> {
      rtsModelOptim<modeltype>& M_;
      double parrho;
      double ll;
    public:
      rho_likelihood(rtsModelOptim<modeltype>& M) :  
        M_(M), parrho(0.0), ll(0.0) {};
      double operator()(const dblvec &par);
    };
    
    class D_likelihood_hsgp : public Functor<dblvec> {
      rtsModelOptim<modeltype>& M;
      double logl;
    public:
      D_likelihood_hsgp(rtsModelOptim<modeltype>& M_) :
        M(M_),
        logl(0.0) {};
      double operator()(const dblvec &par);
    };
    
};

}

template<typename modeltype>
inline void rts::rtsModelOptim<modeltype>::update_theta(const dblvec &theta){
  this->model.covariance.update_parameters(theta);
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
}

template<typename modeltype>
inline void rts::rtsModelOptim<modeltype>::update_u(const MatrixXd& u_){
  if(u_.cols()!=this->re.u(false).cols()){
    this->re.u_.conservativeResize(this->model.covariance.Q(),u_.cols());
    this->re.zu_.resize(this->model.covariance.Q(),u_.cols());
  }
  this->re.u_ = u_;
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
}

template<typename modeltype>
inline void rts::rtsModelOptim<modeltype>::update_rho(const double rho_){
  this->model.covariance.update_rho(rho_);
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
}

template<typename modeltype>
inline void rts::rtsModelOptim<modeltype>::ml_rho(){
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
inline double rts::rtsModelOptim<modeltype>::rho_likelihood::operator()(const dblvec &par){
  parrho = par[0];
  M_.update_rho(parrho);
  ll = M_.log_likelihood();
  return -1*ll;
}

template<typename modeltype>
inline void rts::rtsModelOptim<modeltype>::ml_theta(){
  glmmr::ModelOptim<modeltype>::ml_theta();
}

template<>
inline void rts::rtsModelOptim<BitsHSGP>::ml_theta(){
  D_likelihood_hsgp ddl(*this);
  Rbobyqa<D_likelihood_hsgp,dblvec> opt;
  dblvec lower = this->get_lower_values(false,true,false);
  opt.set_lower(lower);
  opt.control.iprint = trace;
  dblvec start_t = this->get_start_values(false,true,false);
  opt.minimize(ddl, start_t);
}

template<typename modeltype>
inline double rts::rtsModelOptim<modeltype>::D_likelihood_hsgp::operator()(const dblvec &par) {
  M.update_theta(par);
  logl = M.log_likelihood();
  // ArrayXXd u = M.re.u_.array();
  // u = u.square();
  // logl -= u.matrix().sum();
  return -1*logl;
}
