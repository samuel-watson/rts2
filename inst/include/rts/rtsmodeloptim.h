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
    void            update_u(const MatrixXd& u, bool append) override;
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
  // store previous log likelihood values for convergence calculations
  this->previous_ll_values.first = this->current_ll_values.first;
  if(this->ll_previous.rows() != this->ll_current.rows()) this->ll_previous.resize(this->ll_current.rows(),NoChange);
  double old_ll = this->log_likelihood_beta(start);
  this->ll_previous.col(0) = this->ll_current.col(0);

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

  int eval_size = this->control.saem ? this->re.mcmc_block_size : this->ll_current.rows();
  this->current_ll_values.first = this->ll_current.col(0).tail(eval_size).mean();
}

template<typename modeltype>
inline double rts::rtsModelOptim<modeltype>::log_likelihood_beta(const dblvec& beta)
{
  this->model.linear_predictor.update_parameters(beta);
  double ll = this->log_likelihood();
  this->fn_counter.first += this->model.n() * this->re.scaled_u_.cols();
  if(this->control.saem)
  {
    int     iteration = std::max((int)this->re.zu_.cols() / this->re.mcmc_block_size, 1);
    double  gamma = pow(1.0/iteration,this->control.alpha);
    double  ll_t = 0;
    double  ll_pr = 0;
    for(int i = 0; i < iteration; i++){
      int lower_range = i * this->re.mcmc_block_size;
      int upper_range = (i + 1) * this->re.mcmc_block_size;
      if(i == (iteration - 1) && iteration > 1){
        double ll_t_c = ll_t;
        double ll_pr_c = ll_pr;
        ll_t = ll_t + gamma*(this->ll_current.col(0).segment(lower_range, this->re.mcmc_block_size).mean() - ll_t);
        if(this->control.pr_average) ll_pr += ll_t;
        for(int j = lower_range; j < upper_range; j++)
        {
          this->ll_current(j,0) = ll_t_c + gamma*(this->ll_current(j,0) - ll_t_c);
          if(this->control.pr_average) this->ll_current(j,0) = (this->ll_current(j,0) + ll_pr_c)/((double)iteration);
        }
      } else {
        ll_t = ll_t + gamma*(this->ll_current.col(0).segment(lower_range, this->re.mcmc_block_size).mean() - ll_t);
        if(this->control.pr_average) ll_pr += ll_t;
      }
    }
    if(this->control.pr_average){
      ll = ll_pr / (double)iteration;
    } else {
      ll = ll_t;
    }
  } else {
    ll = this->log_likelihood();
  }
  return -1*ll;
}

template<typename modeltype>
inline double rts::rtsModelOptim<modeltype>::log_likelihood_beta_with_gradient(const VectorXd& beta, VectorXd& g)
{
  if(this->control.saem){
    throw std::runtime_error("L-BFGS not available with SAEM");
  } else {
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
}

template<typename modeltype>
template<class algo, typename>
inline void rts::rtsModelOptim<modeltype>::ml_theta()
{  
  dblvec start = this->get_start_values(false,true,false);  
  dblvec lower = this->get_lower_values(false,true,false);
  dblvec upper = this->get_upper_values(false,true,false);
  
  // store previous log likelihood values for convergence calculations
  
  if(this->re.scaled_u_.cols() != this->re.u_.cols())this->re.scaled_u_.resize(NoChange,this->re.u_.cols());
  this->re.scaled_u_ = this->model.covariance.Lu(this->re.u_);  

  if(this->ll_previous.rows() != this->ll_current.rows()) this->ll_previous.resize(this->ll_current.rows(),NoChange);
  this->previous_ll_values.second = this->current_ll_values.second;
  double old_ll = this->log_likelihood_theta(start);
  this->ll_previous.col(1) = this->ll_current.col(1);

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

  // only store the updated log-likelihood if there is no rho step
  if(this->model.covariance.grid.T == 1)
  {
    int eval_size = this->control.saem ? this->re.mcmc_block_size : this->ll_current.rows();
    this->current_ll_values.second = this->ll_current.col(1).tail(eval_size).mean();
  }  
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

  int eval_size = this->control.saem ? this->re.mcmc_block_size : this->ll_current.rows();
  this->current_ll_values.second = this->ll_current.col(1).tail(eval_size).mean();
}

template<typename modeltype>
inline void rts::rtsModelOptim<modeltype>::update_theta(const dblvec &theta)
{
  this->model.covariance.update_parameters(theta);
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
}

template<typename modeltype>
inline double rts::rtsModelOptim<modeltype>::log_likelihood_theta(const dblvec& theta)
{
  this->model.covariance.update_parameters(theta);
  this->fn_counter.second += this->Q() * this->re.scaled_u_.cols();
#pragma omp parallel
  for(int i = 0; i < this->re.scaled_u_.cols(); i++)
  {
    this->ll_current(i,1) = this->model.covariance.log_likelihood(this->re.scaled_u_.col(i));
  }
  double ll = 0;
  if(this->control.saem)
  {
    int     iteration = std::max((int)this->re.zu_.cols() / this->re.mcmc_block_size, 1);
    double  gamma = pow(1.0/iteration,this->control.alpha);
    double  ll_t = 0;
    double  ll_pr = 0;
    for(int i = 0; i < iteration; i++){
      int lower_range = i * this->re.mcmc_block_size;
      int upper_range = (i + 1) * this->re.mcmc_block_size;
      if(i == (iteration - 1) && iteration > 1){
        double ll_t_c = ll_t;
        double ll_pr_c = ll_pr;
        ll_t = ll_t + gamma*(this->ll_current.col(1).segment(lower_range, this->re.mcmc_block_size).mean() - ll_t);
        if(this->control.pr_average) ll_pr += ll_t;
        for(int j = lower_range; j < upper_range; j++)
        {
          this->ll_current(j,1) = ll_t_c + gamma*(this->ll_current(j,1) - ll_t_c);
          if(this->control.pr_average) this->ll_current(j,1) = (this->ll_current(j,1) + ll_pr_c)/((double)iteration);
        }
      } else {
        ll_t = ll_t + gamma*(this->ll_current.col(1).segment(lower_range, this->re.mcmc_block_size).mean() - ll_t);
        if(this->control.pr_average) ll_pr += ll_t;
      }
    }
    if(this->control.pr_average){
      ll = ll_pr / (double)iteration;
    } else {
      ll = ll_t;
    }
  } else {
    ll = this->ll_current.col(1).mean();
  }
  return -1*ll;
}

template<>
inline double rts::rtsModelOptim<BitsHSGP>::log_likelihood_theta(const dblvec& theta){
  this->model.covariance.update_parameters(theta);
  this->fn_counter.second += this->Q() * this->re.scaled_u_.cols();
#pragma omp parallel
  for(int i = 0; i < this->re.scaled_u_.cols(); i++)
  {
    this->ll_current(i,1) = this->model.covariance.log_likelihood(this->re.scaled_u_.col(i));
  }
  double ll = 0;
  if(this->control.saem)
  {
    int     iteration = std::max((int)this->re.zu_.cols() / this->re.mcmc_block_size, 1);
    double  gamma = pow(1.0/iteration,this->control.alpha);
    double  ll_t = 0;
    double  ll_pr = 0;
    for(int i = 0; i < iteration; i++){
      int lower_range = i * this->re.mcmc_block_size;
      int upper_range = (i + 1) * this->re.mcmc_block_size;
      if(i == (iteration - 1) && iteration > 1){
        double ll_t_c = ll_t;
        double ll_pr_c = ll_pr;
        ll_t = ll_t + gamma*(this->ll_current.col(1).segment(lower_range, this->re.mcmc_block_size).mean() - ll_t);
        if(this->control.pr_average) ll_pr += ll_t;
        for(int j = lower_range; j < upper_range; j++)
        {
          this->ll_current(j,1) = ll_t_c + gamma*(this->ll_current(j,1) - ll_t_c);
          if(this->control.pr_average) this->ll_current(j,1) = (this->ll_current(j,1) + ll_pr_c)/((double)iteration);
        }
      } else {
        ll_t = ll_t + gamma*(this->ll_current.col(1).segment(lower_range, this->re.mcmc_block_size).mean() - ll_t);
        if(this->control.pr_average) ll_pr += ll_t;
      }
    }
    if(this->control.pr_average){
      ll = ll_pr / (double)iteration;
    } else {
      ll = ll_t;
    }
  } else {
    ll = this->ll_current.col(1).mean();
  }
  return -1*ll;
}

template<typename modeltype>
inline double rts::rtsModelOptim<modeltype>::log_likelihood_rho(const dblvec& rho)
{
  this->model.covariance.update_rho(rho[0]);
  this->fn_counter.second += this->Q() * this->re.scaled_u_.cols();
#pragma omp parallel
  for(int i = 0; i < this->re.scaled_u_.cols(); i++)
  {
    this->ll_current(i,1) = this->model.covariance.log_likelihood(this->re.scaled_u_.col(i));
  }
  double ll = 0;
  if(this->control.saem)
  {
    int     iteration = std::max((int)this->re.zu_.cols() / this->re.mcmc_block_size, 1);
    double  gamma = pow(1.0/iteration,this->control.alpha);
    double  ll_t = 0;
    double  ll_pr = 0;
    for(int i = 0; i < iteration; i++){
      int lower_range = i * this->re.mcmc_block_size;
      int upper_range = (i + 1) * this->re.mcmc_block_size;
      if(i == (iteration - 1) && iteration > 1){
        double ll_t_c = ll_t;
        double ll_pr_c = ll_pr;
        ll_t = ll_t + gamma*(this->ll_current.col(1).segment(lower_range, this->re.mcmc_block_size).mean() - ll_t);
        if(this->control.pr_average) ll_pr += ll_t;
        for(int j = lower_range; j < upper_range; j++)
        {
          this->ll_current(j,1) = ll_t_c + gamma*(this->ll_current(j,1) - ll_t_c);
          if(this->control.pr_average) this->ll_current(j,1) = (this->ll_current(j,1) + ll_pr_c)/((double)iteration);
        }
      } else {
        ll_t = ll_t + gamma*(this->ll_current.col(1).segment(lower_range, this->re.mcmc_block_size).mean() - ll_t);
        if(this->control.pr_average) ll_pr += ll_t;
      }
    }
    if(this->control.pr_average){
      ll = ll_pr / (double)iteration;
    } else {
      ll = ll_t;
    }
  } else {
    ll = this->ll_current.col(1).mean();
  }
  return -1*ll;
}

template<>
inline double rts::rtsModelOptim<BitsHSGP>::log_likelihood_rho(const dblvec& rho)
{
  this->model.covariance.update_rho(rho[0]);
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
  double ll = this->log_likelihood(false);
  this->fn_counter.first += this->model.n() * this->re.scaled_u_.cols();
  if(this->control.saem)
  {
    int     iteration = std::max((int)this->re.zu_.cols() / this->re.mcmc_block_size, 1);
    double  gamma = pow(1.0/iteration,this->control.alpha);
    double  ll_t = 0;
    double  ll_pr = 0;
    for(int i = 0; i < iteration; i++){
      int lower_range = i * this->re.mcmc_block_size;
      int upper_range = (i + 1) * this->re.mcmc_block_size;
      if(i == (iteration - 1) && iteration > 1){
        double ll_t_c = ll_t;
        double ll_pr_c = ll_pr;
        ll_t = ll_t + gamma*(this->ll_current.col(0).segment(lower_range, this->re.mcmc_block_size).mean() - ll_t);
        if(this->control.pr_average) ll_pr += ll_t;
        for(int j = lower_range; j < upper_range; j++)
        {
          this->ll_current(j,0) = ll_t_c + gamma*(this->ll_current(j,0) - ll_t_c);
          if(this->control.pr_average) this->ll_current(j,0) = (this->ll_current(j,0) + ll_pr_c)/((double)iteration);
        }
      } else {
        ll_t = ll_t + gamma*(this->ll_current.col(0).segment(lower_range, this->re.mcmc_block_size).mean() - ll_t);
        if(this->control.pr_average) ll_pr += ll_t;
      }
    }
    if(this->control.pr_average){
      ll = ll_pr / (double)iteration;
    } else {
      ll = ll_t;
    }
  } else {
    ll = this->log_likelihood(false);
  }
  return -1*ll;
}

template<typename modeltype>
inline double rts::rtsModelOptim<modeltype>::log_likelihood_rho_with_gradient(const VectorXd& rho, VectorXd& g)
{
  if(this->control.saem){
    throw std::runtime_error("L-BFGS-B not available with SAEM");
  } else {
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
  if(this->control.saem){
    throw std::runtime_error("L-BFGS-B not available with SAEM");
  } else {
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
}

template<>
inline double rts::rtsModelOptim<BitsAR>::log_likelihood_theta_with_gradient(const VectorXd& theta, VectorXd& g){
  if(this->control.saem){
    throw std::runtime_error("L-BFGS-B not available with SAEM");
  } else {
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
}

template<>
inline double rts::rtsModelOptim<BitsHSGP>::log_likelihood_rho_with_gradient(const VectorXd& rho, VectorXd& g)
{
  if(this->control.saem){
    throw std::runtime_error("L-BFGS-B not available with SAEM");
  } else {
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
}

template<typename modeltype>
inline void rts::rtsModelOptim<modeltype>::update_u(const MatrixXd& u_, bool append)
{
  int newcolsize = u_.cols();
  int currcolsize = this->re.u_.cols();

  if(append){
    this->re.u_.conservativeResize(NoChange,currcolsize + newcolsize);
    this->re.zu_.conservativeResize(NoChange,currcolsize + newcolsize);
    this->re.u_.rightCols(newcolsize) = u_;
    this->ll_current.resize(currcolsize + newcolsize,NoChange);
  } else {
    if(u_.cols()!=this->re.u_.cols()){
      this->re.u_.resize(NoChange,newcolsize);
      this->re.zu_.resize(NoChange,newcolsize);
      this->re.u_ = u_;
      if(newcolsize != this->ll_current.rows()) this->ll_current.resize(newcolsize,NoChange);
    }
  }
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
}

template<typename modeltype>
inline void rts::rtsModelOptim<modeltype>::update_rho(const double rho_)
{
  this->model.covariance.update_rho(rho_);
  this->re.zu_ = this->model.covariance.ZLu(this->re.u_);
}



