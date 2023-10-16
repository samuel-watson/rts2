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
    MatrixXd hessian(double tol = 1e-4);
    void set_bobyqa_control(int npt_, double rhobeg_, double rhoend_);
    
  private:
    int npt = 0;
    double rhobeg = 0.1;
    double rhoend = 1e-6;
    
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
inline void rts::rtsModelOptim<modeltype>::set_bobyqa_control(int npt_, double rhobeg_, double rhoend_){
  npt = npt_;
  rhobeg = rhobeg_;
  rhoend = rhoend_;
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
  opt.control.iprint = this->trace;
  opt.control.npt = npt;
  opt.control.rhobeg = rhobeg;
  opt.control.rhoend = rhoend;
  dblvec start_t = this->get_start_values(false,true,false);
  opt.minimize(ddl, start_t);
}

template<>
inline void rts::rtsModelOptim<BitsNNGP>::ml_theta(){
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

template<typename modeltype>
inline double rts::rtsModelOptim<modeltype>::D_likelihood_hsgp::operator()(const dblvec &par) {
  M.update_theta(par);
  logl = M.log_likelihood();
  return -1*logl;
}

template<typename modeltype>
inline MatrixXd rts::rtsModelOptim<modeltype>::hessian(double tol){
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
