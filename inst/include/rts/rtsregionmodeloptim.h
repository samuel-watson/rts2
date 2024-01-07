# pragma once

#include <glmmr/modeloptim.hpp>
#include "rtsmodelbits.h"
#include "griddata.h"
#include "regiondata.h"
#include "regionlinearpredictor.h"

namespace rts {

using namespace Eigen;
using namespace glmmr;

template<typename modeltype>
class rtsRegionModelOptim : public ModelOptim<modeltype> {
public:
  rts::RegionData&  region;
  
  rtsRegionModelOptim(modeltype& model_, 
                glmmr::ModelMatrix<modeltype>& matrix_,
                glmmr::RandomEffects<modeltype>& re_,
                rts::RegionData& region_) : ModelOptim<modeltype>(model_,matrix_,re_), region(region_) {};
  
  rtsRegionModelOptim(const rts::rtsRegionModelOptim<modeltype>& optim) : ModelOptim<modeltype>(optim.model, optim.matrix, optim.re), region(optim.region) {};
  
  void        update_theta(const dblvec &theta) override;
  void        update_u(const MatrixXd& u) override;
  void        update_rho(double rho);
  template<class algo, typename = std::enable_if_t<std::is_base_of<optim_algo, algo>::value> >
  void        ml_beta();
  template<class algo, typename = std::enable_if_t<std::is_base_of<optim_algo, algo>::value> >
  void        ml_theta();
  template<class algo, typename = std::enable_if_t<std::is_base_of<optim_algo, algo>::value> >
  void        ml_rho();
  double      log_likelihood_rho(const dblvec &rho);
  double      log_likelihood_rho_with_gradient(const VectorXd &rho, VectorXd& g);
  double      log_likelihood_beta(const dblvec &beta);
  double      log_likelihood() override;
  double      full_log_likelihood() override;
  ArrayXXd    region_intensity(bool uselog = true);
  ArrayXXd    y_predicted(bool uselog = true);
};

}


template<typename modeltype>
template<class algo, typename>
inline void rts::rtsRegionModelOptim<modeltype>::ml_beta()
{  
  dblvec start = this->get_start_values(true,false,false);
  if constexpr (std::is_same_v<algo,LBFGS>)
  {
    throw std::runtime_error("L-BGFS not available with regional data model yet.");
  } else {
    optim<double(const std::vector<double>&),algo> op(start);
    if constexpr (std::is_same_v<algo,DIRECT>) {
      op.set_bounds(start,dblvec(start.size(),this->control.direct_range_beta),true);
      this->set_direct_control(op);
    } else if constexpr (std::is_same_v<algo,BOBYQA>) {
      this->set_bobyqa_control(op);
    } else if constexpr (std::is_same_v<algo,NEWUOA>) {
      this->set_newuoa_control(op);
    }
    if(this->beta_bounded) op.set_bounds(this->lower_bound,this->upper_bound);
    if constexpr (std::is_same_v<modeltype,BitsAR>)
    {
      op.template fn<&rts::rtsRegionModelOptim<BitsAR>::log_likelihood_beta, rts::rtsRegionModelOptim<BitsAR> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsNNGP>) {
      op.template fn<&rts::rtsRegionModelOptim<BitsNNGP>::log_likelihood_beta, rts::rtsRegionModelOptim<BitsNNGP> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsHSGP>) {
      op.template fn<&rts::rtsRegionModelOptim<BitsHSGP>::log_likelihood_beta, rts::rtsRegionModelOptim<BitsHSGP> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsARRegion>) {
      op.template fn<&rts::rtsRegionModelOptim<BitsARRegion>::log_likelihood_beta, rts::rtsRegionModelOptim<BitsARRegion> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsNNGPRegion>) {
      op.template fn<&rts::rtsRegionModelOptim<BitsNNGPRegion>::log_likelihood_beta, rts::rtsRegionModelOptim<BitsNNGPRegion> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsHSGPRegion>){
      op.template fn<&rts::rtsRegionModelOptim<BitsHSGPRegion>::log_likelihood_beta, rts::rtsRegionModelOptim<BitsHSGPRegion> >(this);
    }
    op.minimise();
  }
}

template<typename modeltype>
template<class algo, typename>
inline void rts::rtsRegionModelOptim<modeltype>::ml_theta(){  
  dblvec start = this->get_start_values(false,true,false);  
  dblvec lower = this->get_lower_values(false,true,false);
  dblvec upper = this->get_upper_values(false,true,false);
  if(this->re.scaled_u_.cols() != this->re.u_.cols())this->re.scaled_u_.conservativeResize(NoChange,this->re.u_.cols());
  this->re.scaled_u_ = this->model.covariance.Lu(this->re.u_);  
  if constexpr (std::is_same_v<algo,LBFGS>){
    VectorXd start_vec = Map<VectorXd>(start.data(),start.size());
    optim<double(const VectorXd&, VectorXd&),algo> op(start_vec); 
    op.set_bounds(lower,upper);
    this->set_lbfgs_control(op);
    if constexpr (std::is_same_v<modeltype,BitsAR>) 
    {
      op.template fn<&rts::rtsRegionModelOptim<BitsAR>::log_likelihood_theta_with_gradient, rts::rtsRegionModelOptim<BitsAR> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsARRegion>) {
      op.template fn<&rts::rtsRegionModelOptim<BitsARRegion>::log_likelihood_theta_with_gradient, rts::rtsRegionModelOptim<BitsARRegion> >(this);
    } else {
      throw std::runtime_error("L-BFGS not available for this model type");
    }
    op.minimise();
  } else {
    optim<double(const std::vector<double>&),algo> op(start);
    if constexpr (std::is_same_v<algo,DIRECT>) {      
      dblvec upper2(lower.size());
      std::fill(upper2.begin(),upper2.end(),1.0);
      op.set_bounds(lower,upper2,false);
      this->set_direct_control(op);
    } else if constexpr (std::is_same_v<algo,BOBYQA>) {
      this->set_bobyqa_control(op);
      op.set_bounds(lower,upper);
    } else if constexpr (std::is_same_v<algo,NEWUOA>) {
      this->set_newuoa_control(op);
      op.set_bounds(lower,upper);
    }
    if constexpr (std::is_same_v<modeltype,BitsAR>)
    {
      op.template fn<&rts::rtsRegionModelOptim<BitsAR>::log_likelihood_theta, rts::rtsRegionModelOptim<BitsAR> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsNNGP>) {
      op.template fn<&rts::rtsRegionModelOptim<BitsNNGP>::log_likelihood_theta, rts::rtsRegionModelOptim<BitsNNGP> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsHSGP>) {
      op.template fn<&rts::rtsRegionModelOptim<BitsHSGP>::log_likelihood_theta, rts::rtsRegionModelOptim<BitsHSGP> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsARRegion>) {
      op.template fn<&rts::rtsRegionModelOptim<BitsARRegion>::log_likelihood_theta, rts::rtsRegionModelOptim<BitsARRegion> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsNNGPRegion>) {
      op.template fn<&rts::rtsRegionModelOptim<BitsNNGPRegion>::log_likelihood_theta, rts::rtsRegionModelOptim<BitsNNGPRegion> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsHSGPRegion>){
      op.template fn<&rts::rtsRegionModelOptim<BitsHSGPRegion>::log_likelihood_theta, rts::rtsRegionModelOptim<BitsHSGPRegion> >(this);
    }
    op.minimise();
  }
}

template<typename modeltype>
template<class algo, typename>
inline void rts::rtsRegionModelOptim<modeltype>::ml_rho()
{  
  dblvec start;
  start.push_back(this->model.covariance.rho);
  dblvec lower;
  lower.push_back(-1.0);
  dblvec upper;
  upper.push_back(1.0);
  if(this->re.scaled_u_.cols() != this->re.u_.cols())this->re.scaled_u_.conservativeResize(NoChange,this->re.u_.cols());
  this->re.scaled_u_ = this->model.covariance.Lu(this->re.u_);  
  if constexpr (std::is_same_v<algo,LBFGS>){
    VectorXd start_vec = Map<VectorXd>(start.data(),start.size());
    optim<double(const VectorXd&, VectorXd&),algo> op(start_vec); 
    op.set_bounds(lower,upper);
    this->set_lbfgs_control(op);
    if constexpr (std::is_same_v<modeltype,BitsAR>) 
    {
      op.template fn<&rts::rtsRegionModelOptim<BitsAR>::log_likelihood_rho_with_gradient, rts::rtsRegionModelOptim<BitsAR> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsARRegion>) {
      op.template fn<&rts::rtsRegionModelOptim<BitsARRegion>::log_likelihood_rho_with_gradient, rts::rtsRegionModelOptim<BitsARRegion> >(this);
    }  else {
      throw std::runtime_error("L-BFGS not available for this model type");
    }
    op.minimise();
  } else {
    optim<double(const std::vector<double>&),algo> op(start);
    if constexpr (std::is_same_v<algo,DIRECT>) {      
      op.set_bounds(lower,upper,false);
      this->set_direct_control(op);
    } else if constexpr (std::is_same_v<algo,BOBYQA>) {
      this->set_bobyqa_control(op);
      op.set_bounds(lower,upper);
    } else if constexpr (std::is_same_v<algo,NEWUOA>) {
      this->set_newuoa_control(op);
      op.set_bounds(lower,upper);
    }
    if constexpr (std::is_same_v<modeltype,BitsAR>)
    {
      op.template fn<&rts::rtsRegionModelOptim<BitsAR>::log_likelihood_rho, rts::rtsRegionModelOptim<BitsAR> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsNNGP>) {
      op.template fn<&rts::rtsRegionModelOptim<BitsNNGP>::log_likelihood_rho, rts::rtsRegionModelOptim<BitsNNGP> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsHSGP>) {
      op.template fn<&rts::rtsRegionModelOptim<BitsHSGP>::log_likelihood_rho, rts::rtsRegionModelOptim<BitsHSGP> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsARRegion>) {
      op.template fn<&rts::rtsRegionModelOptim<BitsARRegion>::log_likelihood_rho, rts::rtsRegionModelOptim<BitsARRegion> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsNNGPRegion>) {
      op.template fn<&rts::rtsRegionModelOptim<BitsNNGPRegion>::log_likelihood_rho, rts::rtsRegionModelOptim<BitsNNGPRegion> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsHSGPRegion>){
      op.template fn<&rts::rtsRegionModelOptim<BitsHSGPRegion>::log_likelihood_rho, rts::rtsRegionModelOptim<BitsHSGPRegion> >(this);
    }
    op.minimise();
  }
}

template<typename modeltype>
inline double rts::rtsRegionModelOptim<modeltype>::log_likelihood_rho(const dblvec& rho){
  update_rho(rho[0]);
  double logl = 0;
  #pragma omp parallel for reduction (+:logl)
    for(int i = 0; i < this->re.scaled_u_.cols(); i++)
    {
      logl += this->model.covariance.log_likelihood(this->re.scaled_u_.col(i));
    }
  return -1*logl;
}

template<>
inline double rts::rtsRegionModelOptim<BitsHSGP>::log_likelihood_rho(const dblvec& rho){
  update_rho(rho[0]);
  double logl = log_likelihood();
  return -1*logl;
}

template<>
inline double rts::rtsRegionModelOptim<BitsHSGPRegion>::log_likelihood_rho(const dblvec& rho){
  update_rho(rho[0]);
  double logl = log_likelihood();
  return -1*logl;
}

template<typename modeltype>
inline double rts::rtsRegionModelOptim<modeltype>::log_likelihood_rho_with_gradient(const VectorXd& rho, VectorXd& g)
{
  update_rho(rho(0));
  double logl = 0;
  #pragma omp parallel for reduction (+:logl)
  for(int i = 0; i < this->re.scaled_u_.cols(); i++)
    {
      logl += this->model.covariance.log_likelihood(this->re.scaled_u_.col(i));
    }
  g = this->model.covariance.log_gradient_rho(this->re.scaled_u_);
  g.array() *= -1.0;
  return -1*logl;
}

template<typename modeltype>
inline double rts::rtsRegionModelOptim<modeltype>::log_likelihood_beta(const dblvec& beta)
{
  this->model.linear_predictor.update_parameters(beta);
  double ll = this->log_likelihood();
  return -1*ll;
}

template<typename modeltype>
inline void rts::rtsRegionModelOptim<modeltype>::update_theta(const dblvec &theta)
{
  this->model.covariance.update_parameters(theta);
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
}

template<typename modeltype>
inline void rts::rtsRegionModelOptim<modeltype>::update_rho(double rho)
{
  this->model.covariance.update_rho(rho);
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
}

template<typename modeltype>
inline void rts::rtsRegionModelOptim<modeltype>::update_u(const MatrixXd& u_)
{
  if(u_.cols()!=this->re.u(false).cols()){
    this->re.u_.conservativeResize(this->model.covariance.Q(),u_.cols());
    this->re.zu_.resize(this->model.covariance.Q(),u_.cols());
  }
  this->re.u_ = u_;
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
}

template<typename modeltype>
inline double rts::rtsRegionModelOptim<modeltype>::log_likelihood()
{
  double ll = 0;
  ArrayXXd xb = y_predicted(true);
  if(this->model.weighted){
#pragma omp parallel for reduction (+:ll) collapse(2)
    for(int j=0; j<xb.cols() ; j++){
      for(int i = 0; i<xb.rows(); i++){
        ll += this->model.data.weights(i)*glmmr::maths::log_likelihood(this->model.data.y(i),xb(i,j),this->model.data.variance(i),this->model.family.family,this->model.family.link);
      }
    }
    ll *= this->model.data.weights.sum()/this->model.n();
  } else {
#pragma omp parallel for reduction (+:ll) collapse(2)
  for(int j=0; j<xb.cols() ; j++){
    for(int i = 0; i< xb.rows(); i++){
      ll += glmmr::maths::log_likelihood(this->model.data.y(i),xb(i,j),this->model.data.variance(i),this->model.family.family,this->model.family.link);
    }
  }
}
  return ll/xb.cols();
}

template<typename modeltype>
inline ArrayXXd rts::rtsRegionModelOptim<modeltype>::y_predicted(bool uselog)
{
  ArrayXXd xb(this->model.n(), this->re.u_.cols());
  if constexpr (std::is_same_v<modeltype, BitsAR> || std::is_same_v<modeltype, BitsNNGP > || std::is_same_v<modeltype, BitsHSGP >){
    xb = region_intensity(true);
  } else if constexpr (std::is_same_v<modeltype, BitsARRegion > || std::is_same_v<modeltype, BitsNNGPRegion > || std::is_same_v<modeltype, BitsHSGPRegion >){
    xb = this->model.linear_predictor.xb_region(this->re.zu_);
  }
  xb.matrix().colwise() += this->model.data.offset;
  if(!uselog)xb = xb.exp();
  return xb;
}

template<typename modeltype>
inline ArrayXXd rts::rtsRegionModelOptim<modeltype>::region_intensity(bool uselog)
{
  MatrixXd regionu = region.grid_to_region(this->re.zu_);
  ArrayXXd intens = ArrayXXd::Zero(region.nRegion * region.gridT,this->re.u_.cols());
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
inline double rts::rtsRegionModelOptim<modeltype>::full_log_likelihood()
{
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

