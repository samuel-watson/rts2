#pragma once

#include <chrono>
#include <glmmr/modelbits.hpp>


namespace rts {

using namespace Eigen;
using namespace glmmr;

template<typename cov>
class regionModel {
public:
  // region data and weights stored as a sparse matrix
  SparseMatrix<double>  weights;
  Formula               formula;
  cov                   covariance;
  MatrixXd              X;
  ArrayXd               y;
  int                   Q;
  int                   niter;
  Family                fam;
  VectorXd              beta;
  int                   trace = 1;
  
  template <typename C = cov, typename = std::enable_if_t<std::is_same_v<C, glmmr::Covariance> > >
  regionModel(const std::string& formula_,
              const ArrayXXd& data_,
              const strvec& colnames_,
              const MatrixXd& X_,
              const ArrayXd& y_,
              const int niter_) : formula(formula_), covariance(formula_,data_,colnames_),
              X(X_), y(y_), Q(covariance.Q()), niter(niter_), fam("poisson","log"), 
              beta(VectorXd::Zero(X_.cols())), 
              u_(MatrixXd::Zero(covariance.Q(), niter_)),
              scaled_u_(MatrixXd::Zero(covariance.Q(), niter_)), 
              u_mean_(VectorXd::Zero(covariance.Q())), 
              u_solve_(MatrixXd::Zero(covariance.Q(), niter_)), 
              u_weight_(ArrayXd::Zero(niter_)), 
              u_loglik_(VectorXd::Zero(niter_)),
              gradients(ArrayXd::Zero(X_.cols() + covariance.npar())), M(MatrixXd::Zero(X_.cols(),X_.cols())),
              offset(VectorXd::Zero(X_.rows())) {};
  
#ifdef GLMMR13
  template <typename C = cov, typename = std::enable_if_t<std::is_same_v<C, glmmr::ar1Covariance> > >
  regionModel(const std::string& formula_,
              const ArrayXXd& data_,
              const strvec& colnames_,
              const MatrixXd& X_,
              const ArrayXd& y_,
              const int niter_,
              const int T) : formula(formula_), covariance(formula_,data_,colnames_,T),
              X(X_), y(y_), Q(covariance.Q()), niter(niter_), fam("poisson","log"), 
              beta(VectorXd::Zero(X_.cols())), 
              u_(MatrixXd::Zero(covariance.Q(), niter_)),
              scaled_u_(MatrixXd::Zero(covariance.Q(), niter_)), 
              u_mean_(VectorXd::Zero(covariance.Q())), 
              u_solve_(MatrixXd::Zero(covariance.Q(), niter_)), 
              u_weight_(ArrayXd::Zero(niter_)), 
              u_loglik_(VectorXd::Zero(niter_)),
              gradients(ArrayXd::Zero(X_.cols() + covariance.npar())), M(MatrixXd::Zero(X_.cols(),X_.cols())),
              offset(VectorXd::Zero(X_.rows())) {};
#endif
  
  template <typename C = cov, typename = std::enable_if_t<std::is_same_v<C, glmmr::hsgpCovariance> > >
  regionModel(const std::string& formula_,
              const ArrayXXd& data_,
              const strvec& colnames_,
              const MatrixXd& X_,
              const ArrayXd& y_,
              const int niter_,
              const intvec& m_,
              const double L_boundary_) : formula(formula_), 
              covariance(formula_,data_,colnames_),
              X(X_), y(y_), Q(0), niter(niter_), fam("poisson","log"), 
              beta(VectorXd::Zero(X_.cols())), 
              u_(MatrixXd::Zero(1, niter_)),
              scaled_u_(MatrixXd::Zero(1, niter_)), 
              u_mean_(VectorXd::Zero(1)), 
              u_solve_(MatrixXd::Zero(1, niter_)), 
              u_weight_(ArrayXd::Zero(niter_)), 
              u_loglik_(VectorXd::Zero(niter_)),
              gradients(ArrayXd::Zero(X_.cols() + covariance.npar())), 
              M(MatrixXd::Zero(X_.cols(),X_.cols())),
              offset(VectorXd::Zero(X_.rows())) 
  {
    covariance.update_approx_parameters(m_, L_boundary_);
    Q = covariance.Q();
    // Resize to correct spectral dimension
    u_.resize(Q, niter_);         u_.setZero();
    scaled_u_.resize(X_.rows(), niter_); scaled_u_.setZero();
    u_mean_.resize(Q);            u_mean_.setZero();
    u_solve_.resize(Q, niter_);   u_solve_.setZero();
  };
  
  void              init_beta();
  void              usample(const int niter);
  MatrixXd          lambda_r();
  ArrayXd           log_likelihood();
  void              nr_beta();
  void              nr_theta();
  void              fit(const double tol, const int max_iter, const int hist, const int k0);
  bool             check_convergence(const double tol, const int hist, const int k, const int k0);
  MatrixXd         u() const;
  MatrixXd         information_matrix() const;
  ArrayXd          sampling_weights() const;
  double           total_log_likelihood() const;
  void             set_offset(const VectorXd& offset_);
  
private:
  
  MatrixXd    u_;
  MatrixXd    scaled_u_;
  VectorXd    u_mean_;
  MatrixXd    u_solve_;
  ArrayXd     u_weight_;
  VectorXd    u_loglik_;
  dblvec      logliks;
  double      ll_beta;
  double      ll_theta;
public:
  ArrayXd                           gradients;
private:
  std::deque<double>                gradient_history;
  dblvec                            converge_z;
  dblvec                            converge_bf;
  double      quantile;
  MatrixXd    M;
  VectorXd    offset;
  
};

}