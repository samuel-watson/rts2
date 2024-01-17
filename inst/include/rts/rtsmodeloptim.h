# pragma once

#include <glmmr/modeloptim.hpp>
#include "rtsmodelbits.h"
#include "regionlinearpredictor.h"

namespace rts {

using namespace Eigen;
using namespace glmmr;

template<typename modeltype>
class rtsModelOptim : public ModelOptim<modeltype> {  public:
    
    rtsModelOptim(modeltype& model_,  glmmr::ModelMatrix<modeltype>& matrix_, glmmr::RandomEffects<modeltype>& re_) : ModelOptim<modeltype>(model_,matrix_,re_) {};    
    rtsModelOptim(const rts::rtsModelOptim<modeltype>& optim) : ModelOptim<modeltype>(optim.model,optim.matrix,optim.re) {};    
    
    void            update_theta(const dblvec &theta) override;
    void            update_u(const MatrixXd& u) override;
    void            update_rho(const double rho_);
    template<class algo, typename = std::enable_if_t<std::is_base_of<optim_algo, algo>::value> >
    void            ml_beta();
    template<class algo, typename = std::enable_if_t<std::is_base_of<optim_algo, algo>::value> >
    void            ml_theta();
    template<class algo, typename = std::enable_if_t<std::is_base_of<optim_algo, algo>::value> >
    void            ml_rho();
    double          log_likelihood_rho(const dblvec &rho);
    double          log_likelihood_rho_with_gradient(const VectorXd &rho, VectorXd& g);
    double          log_likelihood_beta(const dblvec &beta);
    double          log_likelihood_theta(const dblvec &theta);
    double          log_likelihood_beta_with_gradient(const VectorXd &beta, VectorXd& g);
    double          log_likelihood_theta_with_gradient(const VectorXd& theta, VectorXd& g);
    
};

}


template<typename modeltype>
template<class algo, typename>
inline void rts::rtsModelOptim<modeltype>::ml_beta()
{  
  dblvec start = this->get_start_values(true,false,false);
  if constexpr (std::is_same_v<algo,LBFGS>){
    VectorXd start_vec = Map<VectorXd>(start.data(),start.size());
    optim<double(const VectorXd&, VectorXd&),algo> op(start_vec);
    this->set_lbfgs_control(op);
    if(this->beta_bounded) op.set_bounds(this->lower_bound,this->upper_bound);
      if constexpr (std::is_same_v<modeltype,BitsAR>) {
      op.template fn<&rts::rtsModelOptim<BitsAR>::log_likelihood_beta_with_gradient, rts::rtsModelOptim<BitsAR> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsNNGP>) {
      op.template fn<&rts::rtsModelOptim<BitsNNGP>::log_likelihood_beta_with_gradient, rts::rtsModelOptim<BitsNNGP> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsHSGP>){
      op.template fn<&rts::rtsModelOptim<BitsHSGP>::log_likelihood_beta_with_gradient, rts::rtsModelOptim<BitsHSGP> >(this);
    }
    op.minimise();
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
    if constexpr (std::is_same_v<modeltype,BitsAR>) {
      op.template fn<&rts::rtsModelOptim<BitsAR>::log_likelihood_beta, rts::rtsModelOptim<BitsAR> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsNNGP>) {
      op.template fn<&rts::rtsModelOptim<BitsNNGP>::log_likelihood_beta, rts::rtsModelOptim<BitsNNGP> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsHSGP>){
      op.template fn<&rts::rtsModelOptim<BitsHSGP>::log_likelihood_beta, rts::rtsModelOptim<BitsHSGP> >(this);
    }
    op.minimise();
  }
}

template<typename modeltype>
inline double rts::rtsModelOptim<modeltype>::log_likelihood_beta(const dblvec& beta)
{
  this->model.linear_predictor.update_parameters(beta);
  double ll = this->log_likelihood();
  return -1*ll;
}

template<typename modeltype>
inline double rts::rtsModelOptim<modeltype>::log_likelihood_beta_with_gradient(const VectorXd& beta, VectorXd& g)
{
  this->model.linear_predictor.update_parameters(beta.array());
  MatrixXd grad(g.size(),this->re.u_.cols());
  ArrayXd xb = this->model.xb();
#pragma omp parallel for
  for(int i = 0; i < this->re.u_.cols(); i++)
  { 
    ArrayXd mu = xb + (this->re.zu_.col(i)).array();
    mu = mu.exp();
    grad.col(i) = this->model.linear_predictor.X().transpose()*(this->model.data.y-mu.matrix());
  }  
  g = grad.rowwise().mean();
  g.array() *= -1.0;
  double ll = this->log_likelihood();
  return -1*ll;
}

template<typename modeltype>
template<class algo, typename>
inline void rts::rtsModelOptim<modeltype>::ml_theta()
{  
  dblvec start = this->get_start_values(false,true,false);  
  dblvec lower = this->get_lower_values(false,true,false);
  dblvec upper = this->get_upper_values(false,true,false);
  // if(this->model.covariance.grid.T > 1)
  // {
  //   start.push_back(this->model.covariance.rho);
  //   lower.push_back(0);
  //   upper.push_back(1);
  // }
  if(this->re.scaled_u_.cols() != this->re.u_.cols())this->re.scaled_u_.conservativeResize(NoChange,this->re.u_.cols());
  this->re.scaled_u_ = this->model.covariance.Lu(this->re.u_);  
  if constexpr (std::is_same_v<algo,LBFGS>){
    VectorXd start_vec = Map<VectorXd>(start.data(),start.size());
    optim<double(const VectorXd&, VectorXd&),algo> op(start_vec); 
    op.set_bounds(lower,upper);
    this->set_lbfgs_control(op);
    if constexpr (std::is_same_v<modeltype,BitsAR>) {
      op.template fn<&rts::rtsModelOptim<BitsAR>::log_likelihood_theta_with_gradient, rts::rtsModelOptim<BitsAR> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsNNGP>) {
      op.template fn<&rts::rtsModelOptim<BitsNNGP>::log_likelihood_theta_with_gradient, rts::rtsModelOptim<BitsNNGP> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsHSGP>){
      op.template fn<&rts::rtsModelOptim<BitsHSGP>::log_likelihood_theta_with_gradient, rts::rtsModelOptim<BitsHSGP> >(this);
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
      op.template fn<&rts::rtsModelOptim<BitsAR>::log_likelihood_theta, rts::rtsModelOptim<BitsAR> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsNNGP>) {
      op.template fn<&rts::rtsModelOptim<BitsNNGP>::log_likelihood_theta, rts::rtsModelOptim<BitsNNGP> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsHSGP>){
      op.template fn<&rts::rtsModelOptim<BitsHSGP>::log_likelihood_theta, rts::rtsModelOptim<BitsHSGP> >(this);
    }
    op.minimise();
  }
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
}

template<typename modeltype>
template<class algo, typename>
inline void rts::rtsModelOptim<modeltype>::ml_rho()
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
    if constexpr (std::is_same_v<modeltype,BitsAR>) {
      op.template fn<&rts::rtsModelOptim<BitsAR>::log_likelihood_rho_with_gradient, rts::rtsModelOptim<BitsAR> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsNNGP>) {
      op.template fn<&rts::rtsModelOptim<BitsNNGP>::log_likelihood_rho_with_gradient, rts::rtsModelOptim<BitsNNGP> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsHSGP>){
      op.template fn<&rts::rtsModelOptim<BitsHSGP>::log_likelihood_rho_with_gradient, rts::rtsModelOptim<BitsHSGP> >(this);
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
      op.template fn<&rts::rtsModelOptim<BitsAR>::log_likelihood_rho, rts::rtsModelOptim<BitsAR> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsNNGP>) {
      op.template fn<&rts::rtsModelOptim<BitsNNGP>::log_likelihood_rho, rts::rtsModelOptim<BitsNNGP> >(this);
    } else if constexpr (std::is_same_v<modeltype,BitsHSGP>){
      op.template fn<&rts::rtsModelOptim<BitsHSGP>::log_likelihood_rho, rts::rtsModelOptim<BitsHSGP> >(this);
    }
    op.minimise();
  }
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
}

template<typename modeltype>
inline void rts::rtsModelOptim<modeltype>::update_theta(const dblvec &theta)
{
  // if(this->model.covariance.grid.T == 1)
  //   {
  //     this->model.covariance.update_parameters(theta);
  //   } else {
  //     dblvec theta_p(2);
  //     theta_p[0] = theta[0];
  //     theta_p[1] = theta[1];
  //     this->model.covariance.update_parameters(theta_p);
  //     this->model.covariance.update_rho(theta[2]);
  //   }
    this->model.covariance.update_parameters(theta);
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
}

template<typename modeltype>
inline double rts::rtsModelOptim<modeltype>::log_likelihood_theta(const dblvec& theta)
{
    // if(this->model.covariance.grid.T == 1)
    // {
    //   this->model.covariance.update_parameters(theta);
    // } else {
    //   dblvec theta_p(2);
    //   theta_p[0] = theta[0];
    //   theta_p[1] = theta[1];
    //   this->model.covariance.update_parameters(theta_p);
    //   this->model.covariance.update_rho(theta[2]);
    // }
    this->model.covariance.update_parameters(theta);
    double logl = 0;
  #pragma omp parallel for reduction (+:logl)
    for(int i = 0; i < this->re.scaled_u_.cols(); i++)
    {
      logl += this->model.covariance.log_likelihood(this->re.scaled_u_.col(i));
    }
    return -1*logl/this->re.u_.cols();
}

template<>
inline double rts::rtsModelOptim<BitsHSGP>::log_likelihood_theta(const dblvec& theta){
  // if(this->model.covariance.grid.T == 1)
  //   {
  //     this->model.covariance.update_parameters(theta);
  //   } else {
  //     dblvec theta_p(2);
  //     theta_p[0] = theta[0];
  //     theta_p[1] = theta[1];
  //     this->model.covariance.update_parameters(theta_p);
  //     this->model.covariance.update_rho(theta[2]);
  //   }
  this->model.covariance.update_parameters(theta);
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
  double ll = this->log_likelihood();
  return -1*ll;
}

template<typename modeltype>
inline double rts::rtsModelOptim<modeltype>::log_likelihood_rho(const dblvec& rho)
{
  this->model.covariance.update_rho(rho[0]);
  double logl = 0;
  int niter = this->re.u_.cols();
  #pragma omp parallel for reduction (+:logl) if(niter > 30)
    for(int i = 0; i < niter; i++)
    {
      logl += this->model.covariance.log_likelihood(this->re.scaled_u_.col(i));
    }
  return -1.0*logl/niter;
}

template<>
inline double rts::rtsModelOptim<BitsHSGP>::log_likelihood_rho(const dblvec& rho)
{
  this->model.covariance.update_rho(rho[0]);
  //update_rho(rho[0]);
  double logl = this->log_likelihood();
  return -1*logl;
}

template<typename modeltype>
inline double rts::rtsModelOptim<modeltype>::log_likelihood_rho_with_gradient(const VectorXd& rho, VectorXd& g)
{
  this->model.covariance.update_rho(rho(0));
  //update_rho(rho(0));
  double logl = 0;
  int niter = this->re.u_.cols();
  #pragma omp parallel for reduction (+:logl) if(niter > 30)
    for(int i = 0; i < niter; i++)
    {
      logl += this->model.covariance.log_likelihood(this->re.scaled_u_.col(i));
    }
  g = this->model.covariance.log_gradient_rho(this->re.scaled_u_);
  g.array() *= -1.0;
  return -1.0*logl/niter;
}

template<>
inline double rts::rtsModelOptim<BitsHSGP>::log_likelihood_theta_with_gradient(const VectorXd& theta, VectorXd& g)
{
  // NOT SURE THIS WORKS, NOT INCLUDING IN FINAL VERSION
  throw std::runtime_error("L-BFGS-B not available for THETA with HSGP currently.");
//   this->model.covariance.update_parameters(theta.array());
//   double ll = this->log_likelihood();
//   ArrayXd xb = this->model.xb();
//   MatrixXd grad(2,this->re.u_.cols());
//   MatrixXd ZLd0 = this->model.covariance.ZL_deriv(0,true);
//   MatrixXd ZLd1 = this->model.covariance.ZL_deriv(1,true);
//   int niter = this->re.u_.cols();
// #pragma omp parallel for if(niter > 30)
//   for(int i = 0; i < niter; i++)
//   { 
//     ArrayXd mu = xb + this->re.zu_.col(i).array();
//     mu = mu.exp();    
//     grad(0,i) = ( (this->model.data.y-mu.matrix()) * (this->re.u_.col(i).transpose()) * ZLd0.transpose()).trace();
//     grad(1,i) = ( (this->model.data.y-mu.matrix()) * (this->re.u_.col(i).transpose()) * ZLd1.transpose()).trace();
//   }  
//   g = grad.rowwise().mean();
//   g.array() *= -1.0;
//   return -1.0*ll;
}

template<>
inline double rts::rtsModelOptim<BitsNNGP>::log_likelihood_theta_with_gradient(const VectorXd& theta, VectorXd& g)
{ 
  // if(this->model.covariance.grid.T == 1)
  // {
  //   this->model.covariance.update_parameters_d(theta.array());
  // } else {
  //   this->model.covariance.update_parameters_d(theta.head(2).array());
  //   this->model.covariance.update_rho(theta(2));
  // }
  this->model.covariance.update_parameters(theta);
  double logl = 0;
  g.head(2) = this->model.covariance.log_gradient(this->re.scaled_u_, logl);
  if(this->model.covariance.grid.T > 1)
  {
    g(2) = this->model.covariance.log_gradient_rho(this->re.scaled_u_)(0);
  }
  g.array() *= -1.0;
  return -1.0*logl;
}

template<>
inline double rts::rtsModelOptim<BitsAR>::log_likelihood_theta_with_gradient(const VectorXd& theta, VectorXd& g){
    // if(this->model.covariance.grid.T == 1)
    // {
    //   this->model.covariance.update_parameters(theta.array());
    // } else {
    //   this->model.covariance.update_parameters(theta.head(2).array());
    //   this->model.covariance.update_rho(theta(2));
    // }
    this->model.covariance.update_parameters(theta);
    double logl = 0;
    g.head(2) = this->model.covariance.log_gradient(this->re.scaled_u_, logl);
    if(this->model.covariance.grid.T > 1)
    {
      g(2) = this->model.covariance.log_gradient_rho(this->re.scaled_u_)(0);
    }
    g.array() *= -1.0;
    return -1*logl;
}

template<>
inline double rts::rtsModelOptim<BitsHSGP>::log_likelihood_rho_with_gradient(const VectorXd& rho, VectorXd& g)
{
  update_rho(rho(0));
  double ll = this->log_likelihood();
  double logl = 0;
  ArrayXd xb = this->model.xb();
  MatrixXd grad(1,this->re.u_.cols());
  MatrixXd ZLd0 = this->model.covariance.ZL_deriv(0,false);
  int niter = this->re.u_.cols();
#pragma omp parallel for if(niter > 30)
  for(int i = 0; i < niter; i++)
  { 
    ArrayXd mu = xb + this->re.zu_.col(i).array();
    mu = mu.exp();    
    grad(0,i) = ((this->model.data.y-mu.matrix()) * (this->re.u_.col(i).transpose()) * ZLd0.transpose()).trace();
  }  
  g = grad.rowwise().mean();
  g.array() *= -1.0;
  return -1.0*ll;
}

template<typename modeltype>
inline void rts::rtsModelOptim<modeltype>::update_u(const MatrixXd& u_)
{
  if(u_.cols()!=this->re.u(false).cols()){
    this->re.u_.conservativeResize(this->model.covariance.Q(),u_.cols());
    this->re.zu_.resize(this->model.covariance.Q(),u_.cols());
  }
  this->re.u_ = u_;
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
}

template<typename modeltype>
inline void rts::rtsModelOptim<modeltype>::update_rho(const double rho_)
{
  this->model.covariance.update_rho(rho_);
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
}



