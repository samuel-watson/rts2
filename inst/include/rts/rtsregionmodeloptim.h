# pragma once

#include <glmmr/modeloptim.hpp>
#include "ar1covariance.h"
#include "nngpcovariance.h"
#include "griddata.h"
#include "regiondata.h"
#include "regionlinearpredictor.h"

namespace glmmr {

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
  
  void update_theta(const dblvec &theta) override;
  void update_u(const MatrixXd& u) override;
  double log_likelihood() override;
  double full_log_likelihood() override;
  ArrayXXd region_intensity(bool uselog = true);
  
  private:
    class rho_likelihood : public Functor<dblvec> {
      rtsRegionModelOptim<modeltype>& M_;
      double parrho;
      double ll;
    public:
      rho_likelihood(rtsRegionModelOptim<modeltype>& M) :  
      M_(*this), parrho(0.0), ll(0.0) {};
      double operator()(const dblvec &par);
    };    
  
};

}

template<typename modeltype>
inline void rts::rtsRegionModelOptim<modeltype>::update_theta(const dblvec &theta){
  this->model.covariance.update_parameters(theta);
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
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
  ArrayXXd xb(this->model.n(), this->re.u_.cols());
  if constexpr (std::is_same_v<modeltype, rts::rtsModelBits<rts::ar1Covariance, glmmr::LinearPredictor> > || std::is_same_v<modeltype, rts::rtsModelBits<rts::nngpCovariance, glmmr::LinearPredictor> >){
    xb = region_intensity();
  } else if constexpr (std::is_same_v<modeltype, rts::rtsModelBits<rts::ar1Covariance, rts::regionLinearPredictor> > || std::is_same_v<modeltype, rts::rtsModelBits<rts::nngpCovariance, rts::regionLinearPredictor> >){
    xb = this->model.linear_predictor.xb_region(this->re.u_);
  }
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
inline ArrayXXd rts::rtsRegionModelOptim<modeltype>::region_intensity(bool uselog){
  MatrixXd regionu = region.grid_to_region(this->re.u_);
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
inline double rts::rtsRegionModelOptim<modeltype>::rho_likelihood::operator()(const dblvec &par){
  parrho = par[0];
  M_.update_rho(parrho);
  ll = M_.log_likelihood();
  return -1*ll;
}