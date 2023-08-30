# pragma once

#include <glmmr/modeloptim.hpp>
#include "ar1covariance.h"
#include "nngpcovariance.h"
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
    
  private:
    class rho_likelihood : public Functor<dblvec> {
      rtsModelOptim<modeltype>& M_;
      double parrho;
      double ll;
    public:
      rho_likelihood(rtsModelOptim<modeltype>& M) :  
        M_(*this), parrho(0.0), ll(0.0) {};
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
}

template<typename modeltype>
inline double rts::rtsModelOptim<modeltype>::rho_likelihood::operator()(const dblvec &par){
  parrho = par[0];
  M_.update_rho(parrho);
  ll = M_.log_likelihood();
  return -1*ll;
}


