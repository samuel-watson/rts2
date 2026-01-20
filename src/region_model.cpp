#include <glmmr/general.h>
#include <glmmr/griddata.hpp>
#include "region_model.h"


template<typename cov>
MatrixXd rts::regionModel<cov>::u() const
{
  return scaled_u_;
}

template<typename cov>
ArrayXd rts::regionModel<cov>::sampling_weights() const 
{
  return u_weight_;
}

template<typename cov>
MatrixXd rts::regionModel<cov>::information_matrix() const 
{
  return M;
}

template<typename cov>
double rts::regionModel<cov>::total_log_likelihood() const 
{
  return ll_beta + ll_theta;
}

template<typename cov>
void rts::regionModel<cov>::set_offset(const VectorXd& offset_) 
{
  if(offset_.size() != offset.size())Rcpp::stop("Offset wrong size");
  offset = offset_;
}

template<typename cov>
MatrixXd rts::regionModel<cov>::lambda_r()
{
  int n_regions = weights.rows();
  int n_cells = weights.cols();

#ifdef GLMMR13
  if constexpr (std::is_same_v<cov, glmmr::ar1Covariance>) {
    int T = covariance.Q() / covariance.Covariance::Q();
    int n_A = covariance.Covariance::Q();  // Spatial dimension
    MatrixXd lr(n_regions * T, u_.cols());
    
    for(int i = 0; i < u_.cols(); i++){
      ArrayXd eta_s = (X * beta + offset).array() + scaled_u_.col(i).array();
      ArrayXd lambda_s = eta_s.exp();
      for(int t = 0; t < T; t++){
        VectorXd lambda_s_t = lambda_s.segment(t * n_A, n_A).matrix();
        lr.col(i).segment(t * n_regions, n_regions) = weights * lambda_s_t;
      }
    }
    return lr;
  } else {
    // Standard spatial case
    MatrixXd lr(n_regions, u_.cols());
    for(int i = 0; i < u_.cols(); i++){
      ArrayXd eta_s = (X * beta + offset).array() + scaled_u_.col(i).array();
      ArrayXd lambda_s = eta_s.exp();
      lr.col(i) = weights * lambda_s.matrix();
    }
    return lr;
  }
#else
  // Standard spatial case
  MatrixXd lr(n_regions, u_.cols());
  for(int i = 0; i < u_.cols(); i++){
    ArrayXd eta_s = (X * beta + offset).array() + scaled_u_.col(i).array();
    ArrayXd lambda_s = eta_s.exp();
    lr.col(i) = weights * lambda_s.matrix();
  }
  return lr;
#endif
}

template<typename cov>
void rts::regionModel<cov>::init_beta() {
  // Initialize intercept to log of mean observed rate
  // Model: log(lambda) = X * beta + offset
  // So: lambda = exp(X * beta + offset)
  // We want: mean(y) ≈ mean_weight * exp(beta(0) + mean_offset)
  // Therefore: beta(0) = log(mean_y / mean_weight) - mean_offset
  
  double mean_y = y.mean();
  double mean_weight = weights.sum() / weights.rows();  // average aggregation factor
  double mean_offset = offset.mean();
  
  if(X.cols() > 0) {
    beta(0) = std::log(std::max(mean_y / mean_weight, 1e-6)) - mean_offset;
  }
}

template<typename cov>
void rts::regionModel<cov>::usample(const int niter_){
  int n_regions = weights.rows();
  int n_cells = weights.cols();
  int T = 1;
  int n_A = covariance.Q();

#ifdef GLMMR13
  if constexpr (std::is_same_v<cov, glmmr::ar1Covariance>) {
    T = covariance.Q() / covariance.Covariance::Q();
    n_A = covariance.Covariance::Q();
  } 
#endif
  
  ArrayXd xb = (X * beta + offset).array();
  ArrayXd Lu_umean = covariance.Lu(u_mean_).array();
  ArrayXd eta_s = xb + Lu_umean;
  
  const MatrixXd ZL = covariance.ZL();
  const MatrixXd ZLt = ZL.transpose();
  const int n_cols = ZL.cols();
  
  VectorXd Mb(n_cols);
  MatrixXd Vb(n_cols, n_cols);
  Vb.setIdentity();
  LLT<MatrixXd> llt_Pb;
  VectorXd yb(n_cols);
  MatrixXd LWL = MatrixXd::Identity(n_cols, n_cols);
  
  VectorXd b = u_mean_;
  VectorXd bnew(b);
  double diff = 1.0;
  int itero = 0;
  MatrixXd C(n_regions * T, n_cols);
  LLT<MatrixXd> llt_CCt;

#ifdef GLMMR13
  if constexpr (std::is_same_v<cov, glmmr::ar1Covariance>) {
    while(diff > 1e-6 && itero < 10) {
      VectorXd Zu = ZL * b;
      eta_s = xb + Zu.array();
      ArrayXd lambda_s = eta_s.exp();
      
      VectorXd lambda_r(n_regions * T);
      for(int t = 0; t < T; t++){
        VectorXd lambda_s_t = lambda_s.segment(t * n_A, n_A).matrix();
        VectorXd result = weights * lambda_s_t;
        lambda_r.segment(t * n_regions, n_regions) = result;
      }
      
      ArrayXd y_over_lambda_r = y / lambda_r.array();
      VectorXd resid_r = (y_over_lambda_r - 1.0).matrix();
      
      MatrixXd A = (ZL.array().colwise() * lambda_s).matrix();
      
      MatrixXd WA(n_regions * T, n_cols);
      for(int t = 0; t < T; t++){
        MatrixXd A_t = A.middleRows(t * n_A, n_A);
        WA.middleRows(t * n_regions, n_regions) = weights * A_t;
      }
      
      ArrayXd inv_sqrt_lambda_r = lambda_r.array().sqrt().inverse();
      C = (WA.array().colwise() * inv_sqrt_lambda_r).matrix();
      VectorXd ZLt_score = WA.transpose() * resid_r;
      
      // === WOODBURY: form CC' instead of C'C ===
      MatrixXd CCt = C * C.transpose();  // (n_regions*T) × (n_regions*T)
      CCt.diagonal().array() += 1.0;
      llt_CCt.compute(CCt);
      
      // Compute yb = C'Cb + ZLt_score
      VectorXd Cb = C * b;
      yb = C.transpose() * Cb + ZLt_score;
      
      // Solve (C'C + I)^{-1} yb using Woodbury:
      // x = yb - C'(CC' + I)^{-1}(C * yb)
      VectorXd Cyb = C * yb;
      VectorXd tmp = llt_CCt.solve(Cyb);
      bnew = yb - C.transpose() * tmp;
      
      diff = (b - bnew).array().abs().maxCoeff();
      itero++;
      b.swap(bnew);
    }
  } else {
    // Standard spatial case
    
    while(diff > 1e-6 && itero < 10) {
      
      VectorXd Zu = ZL * b;
      eta_s = xb + Zu.array();
      ArrayXd lambda_s = eta_s.exp();
      VectorXd lambda_r = weights * lambda_s.matrix();
      ArrayXd y_over_lambda_r = y / lambda_r.array();
      VectorXd resid_r = (y_over_lambda_r - 1.0).matrix();
      MatrixXd A = (ZL.array().colwise() * lambda_s).matrix();
      MatrixXd WA = weights * A;
      ArrayXd inv_sqrt_lambda_r = lambda_r.array().sqrt().inverse();
      C = (WA.array().colwise() * inv_sqrt_lambda_r).matrix();
      VectorXd ZLt_score = WA.transpose() * resid_r;
      
      // Woodbury
      MatrixXd CCt = C * C.transpose();  // n_regions × n_regions
      CCt.diagonal().array() += 1.0;
      llt_CCt.compute(CCt);
      
      VectorXd Cb = C * b;
      yb = C.transpose() * Cb + ZLt_score;
      
      VectorXd Cyb = C * yb;
      VectorXd tmp = llt_CCt.solve(Cyb);
      bnew = yb - C.transpose() * tmp;
      
      diff = (b - bnew).array().abs().maxCoeff();
      itero++;
      b.swap(bnew);
    }
  }
#else
  // Standard spatial case
  while(diff > 1e-6 && itero < 10) {
    
    VectorXd Zu = covariance.Lu(b);
    eta_s = xb + Zu.array();
    ArrayXd lambda_s = eta_s.exp();
    VectorXd lambda_r = weights * lambda_s.matrix();
    ArrayXd y_over_lambda_r = y / lambda_r.array();
    VectorXd resid_r = (y_over_lambda_r - 1.0).matrix();
    MatrixXd A = (ZL.array().colwise() * lambda_s).matrix();
    MatrixXd WA = weights * A;
    ArrayXd inv_sqrt_lambda_r = lambda_r.array().sqrt().inverse();
    C = (WA.array().colwise() * inv_sqrt_lambda_r).matrix();
    VectorXd ZLt_score = WA.transpose() * resid_r;
    
    // Woodbury
    MatrixXd CCt = C * C.transpose();  // n_regions × n_regions
    CCt.diagonal().array() += 1.0;
    llt_CCt.compute(CCt);
    
    VectorXd Cb = C * b;
    yb = C.transpose() * Cb + ZLt_score;
    VectorXd Cyb = C * yb;
    VectorXd tmp = llt_CCt.solve(Cyb);
    bnew = yb - C.transpose() * tmp;
    
    diff = (b - bnew).array().abs().maxCoeff();
    itero++;
    b.swap(bnew);
  }
#endif

  // Common code for both cases
  Mb = b;
  u_mean_ = Mb;
  //llt_Pb.solveInPlace(Vb);
  
  MatrixXd unew(u_.rows(), niter_);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<double> d(0.0, 1.0);
  double* data = unew.data();
  for(int i = 0; i < unew.size(); ++i) {
    data[i] = d(gen);
  }
  
  int m = C.rows();  // n_regions * T
  int n_samples = unew.cols();
  
  // Generate w ~ N(0, I)
  MatrixXd w(m, n_samples);
  // ... fill with N(0,1) random values ...
  data = w.data();
  for(int i = 0; i < w.size(); ++i) {
    data[i] = d(gen);
  }
  
  // RHS = z + C'w
  MatrixXd RHS = unew + C.transpose() * w;
  
  // Apply Woodbury: (C'C + I)^{-1} RHS = RHS - C'(CC'+I)^{-1}(C * RHS)
  MatrixXd CRHS = C * RHS;
  MatrixXd tmp = llt_CCt.solve(CRHS);
  unew = RHS - C.transpose() * tmp;
  
  // Add mean
  //unew.colwise() += bnew;
  
  // LLT<MatrixXd> llt(Vb);
  // MatrixXd LVb = llt.matrixL();
  // //LVb *= 1.5;
  // unew = LVb * unew;
  
  if(u_.cols() != niter_){
    u_.resize(NoChange, niter_);
    u_solve_.resize(NoChange, niter_);
    scaled_u_.resize(NoChange, niter_);
    u_weight_.resize(niter_);
    u_loglik_.resize(niter_);
  }
  
  u_.noalias() = unew;
  
  MatrixXd Cu = C * u_;
  ArrayXd Cu_norms = Cu.colwise().squaredNorm();
  ArrayXd u_norms = u_.colwise().squaredNorm();
  u_loglik_ = -0.5 * (Cu_norms + u_norms);
  
  
  // Compute proposal log-likelihood (no race condition)
// #pragma omp parallel for
//   for(int i = 0; i < u_.cols(); i++){
//     VectorXd v = llt.solve(u_.col(i));
//     u_loglik_(i) = -0.5 * v.dot(u_.col(i));
//   }
  
  u_.colwise() += u_mean_;
  
  // Compute importance weights
  scaled_u_ = covariance.Lu(u_);
#ifdef GLMMR13
  u_solve_ = covariance.solve(scaled_u_);
#else
  MatrixXd D = covariance.D();
  LLT<MatrixXd> llt_D;
  llt_D.compute(D);
  u_solve_ = llt_D.solve(scaled_u_);
#endif
  u_weight_.setZero();
  ArrayXd ll = log_likelihood();
#pragma omp parallel for
  for(int i = 0; i < scaled_u_.cols(); i++){
    double llmod = ll(i);
    double llprior = -0.5 * scaled_u_.col(i).dot(u_solve_.col(i));
    u_weight_(i) = llmod + llprior - u_loglik_(i);
  }
  u_weight_ -= u_weight_.maxCoeff();
  u_weight_ = u_weight_.exp();
  double weightsum = u_weight_.sum();
  u_weight_ *= 1.0 / weightsum;
  // Diagnostic checks
  if(u_mean_.hasNaN()) {
    Rcpp::Rcout << "=== DIAGNOSTIC: NaN in u_mean_ ===" << std::endl;
    Rcpp::Rcout << "y range: [" << y.minCoeff() << ", " << y.maxCoeff() << "]" << std::endl;
    Rcpp::Rcout << "beta: " << beta.transpose() << std::endl;
    Rcpp::Rcout << "offset range: [" << offset.minCoeff() << ", " << offset.maxCoeff() << "]" << std::endl;
    Rcpp::Rcout << "weights dims: " << weights.rows() << " x " << weights.cols() << std::endl;
    Rcpp::Rcout << "weights nonzeros: " << weights.nonZeros() << std::endl;
    MatrixXd lr = lambda_r();
    Rcpp::Rcout << "lambda_r range: [" << lr.minCoeff() << ", " << lr.maxCoeff() << "]" << std::endl;
    R_FlushConsole();
    Rcpp::stop("usample: NaN detected in u_mean_");
  }
  
  if(u_.hasNaN()) {
    Rcpp::Rcout << "=== DIAGNOSTIC: NaN in u_ ===" << std::endl;
    Rcpp::Rcout << "u_mean_ range: [" << u_mean_.minCoeff() << ", " << u_mean_.maxCoeff() << "]" << std::endl;
    Rcpp::Rcout << "Vb may be non-positive-definite" << std::endl;
    R_FlushConsole();
    Rcpp::stop("usample: NaN detected in u_");
  }
  
  if(u_weight_.hasNaN() || !u_weight_.isFinite().all()) {
    Rcpp::Rcout << "=== DIAGNOSTIC: NaN/Inf in u_weight_ ===" << std::endl;
    ArrayXd ll = log_likelihood();
    Rcpp::Rcout << "log_likelihood range: [" << ll.minCoeff() << ", " << ll.maxCoeff() << "]" << std::endl;
    Rcpp::Rcout << "u_loglik_ range: [" << u_loglik_.minCoeff() << ", " << u_loglik_.maxCoeff() << "]" << std::endl;
    Rcpp::Rcout << "scaled_u_ range: [" << scaled_u_.minCoeff() << ", " << scaled_u_.maxCoeff() << "]" << std::endl;
    Rcpp::Rcout << "u_solve_ range: [" << u_solve_.minCoeff() << ", " << u_solve_.maxCoeff() << "]" << std::endl;
    R_FlushConsole();
    Rcpp::stop("usample: NaN/Inf detected in u_weight_");
  }
  
  double ess = 1.0 / (u_weight_ * u_weight_).sum();
  if(ess < 2.0) {
    Rcpp::Rcout << "=== WARNING: Very low ESS ===" << std::endl;
    Rcpp::Rcout << "ESS: " << ess << " / " << u_weight_.size() << std::endl;
    R_FlushConsole();
  }
}

template<typename cov>
ArrayXd rts::regionModel<cov>::log_likelihood()
{
  ArrayXd ll(u_.cols());
  MatrixXd lr = lambda_r();
  
#pragma omp parallel for 
  for(int i = 0; i < scaled_u_.cols(); i++){
    
#ifdef GLMMR13
    ll(i) =  maths::log_likelihood(y.array(),lr.col(i).array().log(), ArrayXd::Constant(y.size(),1.0), fam);
#else
    ArrayXd mu = lr.col(i).array();
    ll(i) = (y.array() * mu - mu.exp()).sum();
#endif
  }
  return ll;
}

template<typename cov>
void rts::regionModel<cov>::nr_beta(){
  int niter = u_.cols();
  int n_regions = weights.rows();
  int n_cells = weights.cols();
  
  MatrixXd zd(X.rows(), niter);
  for(int i = 0; i < niter; i++){
    zd.col(i) = (X * beta + offset).array() + scaled_u_.col(i).array();
  }
  
  MatrixXd XtWXm = MatrixXd::Zero(beta.size(), beta.size());
  VectorXd score = VectorXd::Zero(beta.size());

#ifdef GLMMR13
  if constexpr (std::is_same_v<cov, glmmr::ar1Covariance>) {
    int T = covariance.Q() / covariance.Covariance::Q();
    int n_A = covariance.Covariance::Q();
#pragma omp parallel for   
    for(int i = 0; i < niter; ++i){
      ArrayXd lambda_s = zd.col(i).array().exp();
      
      // Aggregate to regions for each time period
      VectorXd lambda_r(n_regions * T);
      for(int t = 0; t < T; t++){
        VectorXd lambda_s_t = lambda_s.segment(t * n_A, n_A).matrix();
        lambda_r.segment(t * n_regions, n_regions) = weights * lambda_s_t;
      }
      
      VectorXd resid_r = (y / lambda_r.array() - 1.0).matrix();
      
      // B = W * diag(lambda_s) * X, but applied per time period
      // For each time t: B_t = W * diag(lambda_s_t) * X_t
      MatrixXd B(n_regions * T, beta.size());
      for(int t = 0; t < T; t++){
        ArrayXd lambda_s_t = lambda_s.segment(t * n_A, n_A);
        MatrixXd X_t = X.middleRows(t * n_A, n_A);
        MatrixXd lambda_s_X_t = (X_t.array().colwise() * lambda_s_t).matrix();
        B.middleRows(t * n_regions, n_regions) = weights * lambda_s_X_t;
      }
      
      score.noalias() += u_weight_(i) * B.transpose() * resid_r;
      
      ArrayXd inv_lambda_r = lambda_r.array().inverse();
      MatrixXd inv_lambda_r_B = (B.array().colwise() * inv_lambda_r).matrix();
      XtWXm.noalias() += u_weight_(i) * B.transpose() * inv_lambda_r_B;
    }
  } else {
    // Standard spatial case (original code)
#pragma omp parallel for
    for(int i = 0; i < niter; ++i){
      ArrayXd lambda_s = zd.col(i).array().exp();
      VectorXd lambda_r = weights * lambda_s.matrix();
      VectorXd resid_r = (y / lambda_r.array() - 1.0).matrix();
      
      MatrixXd lambda_s_X = (X.array().colwise() * lambda_s).matrix();
      MatrixXd B = weights * lambda_s_X;
      
      score.noalias() += u_weight_(i) * B.transpose() * resid_r;
      
      ArrayXd inv_lambda_r = lambda_r.array().inverse();
      MatrixXd inv_lambda_r_B = (B.array().colwise() * inv_lambda_r).matrix();
      XtWXm.noalias() += u_weight_(i) * B.transpose() * inv_lambda_r_B;
    }
  }
#else
  // Standard spatial case (original code)
#pragma omp parallel for
  for(int i = 0; i < niter; ++i){
    ArrayXd lambda_s = zd.col(i).array().exp();
    VectorXd lambda_r = weights * lambda_s.matrix();
    VectorXd resid_r = (y / lambda_r.array() - 1.0).matrix();
    
    MatrixXd lambda_s_X = (X.array().colwise() * lambda_s).matrix();
    MatrixXd B = weights * lambda_s_X;
    
    score.noalias() += u_weight_(i) * B.transpose() * resid_r;
    
    ArrayXd inv_lambda_r = lambda_r.array().inverse();
    MatrixXd inv_lambda_r_B = (B.array().colwise() * inv_lambda_r).matrix();
    XtWXm.noalias() += u_weight_(i) * B.transpose() * inv_lambda_r_B;
  }
#endif
  
  M = XtWXm;
  Eigen::LLT<MatrixXd> llt(XtWXm);
  gradients.head(X.cols()) = score;
  VectorXd bincr = llt.solve(score);
  beta += bincr;
  ll_beta = log_likelihood().mean();
  
  // Diagnostic checks
  if(beta.hasNaN()) {
    Rcpp::Rcout << "=== DIAGNOSTIC: NaN in beta ===" << std::endl;
    Rcpp::Rcout << "XtWXm diagonal: " << XtWXm.diagonal().transpose() << std::endl;
    Rcpp::Rcout << "XtWXm min eigenvalue: " << Eigen::SelfAdjointEigenSolver<MatrixXd>(XtWXm).eigenvalues().minCoeff() << std::endl;
    Rcpp::Rcout << "score: " << score.transpose() << std::endl;
    Rcpp::Rcout << "u_weight_ range: [" << u_weight_.minCoeff() << ", " << u_weight_.maxCoeff() << "]" << std::endl;
    Rcpp::Rcout << "u_weight_ sum: " << u_weight_.sum() << std::endl;
    Rcpp::Rcout << "scaled_u_ range: [" << scaled_u_.minCoeff() << ", " << scaled_u_.maxCoeff() << "]" << std::endl;
    MatrixXd lr = lambda_r();
    Rcpp::Rcout << "lambda_r range: [" << lr.minCoeff() << ", " << lr.maxCoeff() << "]" << std::endl;
    Rcpp::Rcout << "y range: [" << y.minCoeff() << ", " << y.maxCoeff() << "]" << std::endl;
    R_FlushConsole();
    Rcpp::stop("nr_beta: NaN detected in beta");
  }
  
  if(std::isnan(ll_beta) || std::isinf(ll_beta)) {
    Rcpp::Rcout << "=== DIAGNOSTIC: NaN/Inf in ll_beta ===" << std::endl;
    MatrixXd lr = lambda_r();
    Rcpp::Rcout << "lambda_r range: [" << lr.minCoeff() << ", " << lr.maxCoeff() << "]" << std::endl;
    Rcpp::Rcout << "lambda_r has zeros: " << ((lr.array() <= 0).any() ? "yes" : "no") << std::endl;
    Rcpp::Rcout << "beta: " << beta.transpose() << std::endl;
    Rcpp::Rcout << "y range: [" << y.minCoeff() << ", " << y.maxCoeff() << "]" << std::endl;
    R_FlushConsole();
    Rcpp::stop("nr_beta: NaN/Inf detected in ll_beta");
  }
}

template<typename cov>
bool rts::regionModel<cov>::check_convergence(const double tol, const int hist, const int k, const int k0){
  gradient_history.push_back(ll_theta + ll_beta);
  if(gradient_history.size() > hist) gradient_history.pop_front();
  double diffg = 1;
  double vardiffg = 0.1;
  double z = 10;
  int iter = 0;
  if(gradient_history.size()>1){
    diffg = (gradient_history.back() - gradient_history.front())/(gradient_history.size() - 1);
    if(gradient_history.size()>2){
      vardiffg = 0;
      double a = 0;
      for(const double x: gradient_history){
        if(iter > 0) {
          vardiffg += pow((x - a) - diffg, 2);
        } 
        a = x;
        iter++;
      }
      vardiffg *= 1.0/ (gradient_history.size() - 2);
      z = diffg * sqrt(gradient_history.size()) / sqrt(vardiffg);
      converge_z.push_back(z);
    }
  }
  if(trace > 0)Rcpp::Rcout << "\nMean: " << diffg << " sd: " << sqrt(vardiffg) << "\nZ (diff): " << z;
  double prior0 = 1.0 - exp(-(k*k/(k0*k0)));// squared weibull
  double p = maths::gaussian_cdf(z);
  double bf = (1-p)*prior0/(p*(1-prior0));
  if(gradient_history.size()>2)converge_bf.push_back(bf);
  if(trace > 0)Rcpp::Rcout << "\nBF: " << bf << " prior: " << prior0;
  double meang = 0;
  for(const double x: gradient_history) meang += x;
  meang *= 1.0/gradient_history.size();
  double diff = quantile == 0 ? 10 : meang - quantile;
  if(trace > 0)Rcpp::Rcout << "\nGradients: " << gradients.transpose() << " | ||G|| = " << gradients.matrix().norm();
  if(trace > 0)Rcpp::Rcout << "\nLog-likelihood running mean: " << meang << " diff: " << diff;
  quantile = meang;
  return bf > tol;
}

template<typename cov>
void rts::regionModel<cov>::nr_theta(){
  ArrayXd  tmp(u_.cols());
#ifdef GLMMR13
  covariance.nr_step(scaled_u_, u_solve_, tmp, gradients, u_weight_);
  ll_theta = tmp.mean();
#else
  Rcpp::Rcout << "If you're seeing this message you need to recompile against glmmrBase >=1.2.0\nThe covariance parameters are NOT updated" << std::endl;
  ll_theta = 0;
#endif
  
  
  // Diagnostic checks
  bool theta_has_nan = false;
  for(size_t i = 0; i < covariance.parameters_.size(); i++) {
    if(std::isnan(covariance.parameters_[i])) {
      theta_has_nan = true;
      break;
    }
  }
#ifdef GLMMR13
  if constexpr (std::is_same_v<cov, glmmr::ar1Covariance>) {
    if(std::isnan(covariance.rho)) {
      theta_has_nan = true;
    }
  }
#endif
  
  if(theta_has_nan) {
    Rcpp::Rcout << "=== DIAGNOSTIC: NaN in theta ===" << std::endl;
    Rcpp::Rcout << "theta: ";
    for(size_t i = 0; i < covariance.parameters_.size(); i++) {
      Rcpp::Rcout << covariance.parameters_[i] << " ";
    }
    Rcpp::Rcout << std::endl;
#ifdef GLMMR13
    if constexpr (std::is_same_v<cov, glmmr::ar1Covariance>) {
      Rcpp::Rcout << "rho: " << covariance.rho << std::endl;
    }
#endif
    Rcpp::Rcout << "gradients: " << gradients.transpose() << std::endl;
#ifdef GLMMR13
    Rcpp::Rcout << "infomat_theta diagonal: " << covariance.infomat_theta.diagonal().transpose() << std::endl;
#endif
    Rcpp::Rcout << "u_weight_ range: [" << u_weight_.minCoeff() << ", " << u_weight_.maxCoeff() << "]" << std::endl;
    Rcpp::Rcout << "scaled_u_ range: [" << scaled_u_.minCoeff() << ", " << scaled_u_.maxCoeff() << "]" << std::endl;
    Rcpp::Rcout << "u_solve_ range: [" << u_solve_.minCoeff() << ", " << u_solve_.maxCoeff() << "]" << std::endl;
    R_FlushConsole();
    Rcpp::stop("nr_theta: NaN detected in theta parameters");
  }
  
  if(std::isnan(ll_theta) || std::isinf(ll_theta)) {
    Rcpp::Rcout << "=== DIAGNOSTIC: NaN/Inf in ll_theta ===" << std::endl;
    Rcpp::Rcout << "tmp (log-likelihoods) may contain invalid values" << std::endl;
    R_FlushConsole();
    Rcpp::stop("nr_theta: NaN/Inf detected in ll_theta");
  }
}

template<typename cov>
void rts::regionModel<cov>::fit(const double tol, const int max_iter, const int hist, const int k0)
{
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;
  
  bool converged = false;    
  int iter = 1;
  dblpair ll, lldiff;
  double lltot, llvartot, prob;
  // start at ml values
  if(trace > 0)Rcpp::Rcout << "\nStarting values\nBeta: " << beta.transpose() << "\ntheta: ";
  if(trace > 0)for(const auto& i: covariance.parameters_) Rcpp::Rcout << " " << i;
  if(trace > 0)Rcpp::Rcout << " || " << std::endl;
  
  while (!converged && iter <= max_iter) {
    if(trace > 0)Rcpp::Rcout << "\n-------------- ITER: " << iter << " ------------" << std::endl;
    auto t1 = high_resolution_clock::now();
    usample(niter);
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    if(trace > 0)Rcpp::Rcout << "TIMING STEP 1 (posterior u sample): " << ms_double.count() << "ms" << std::endl;
    nr_beta();
    auto t3 = high_resolution_clock::now();
    ms_double = t3 - t2;
    if(trace > 0)Rcpp::Rcout << "TIMING STEP 2 (nr beta): " << ms_double.count() << "ms" << std::endl;
    nr_theta();
    auto t4 = high_resolution_clock::now();
    ms_double = t4 - t3;
    if(trace > 0)Rcpp::Rcout << "TIMING STEP 3 (nr theta): " << ms_double.count() << "ms" << std::endl;
    if(trace > 0)Rcpp::Rcout << "\nBeta: " << beta.transpose() << "\ntheta: ";
    if(trace > 0)for(const auto& i: covariance.parameters_) Rcpp::Rcout << " " << i;
#ifdef GLMMR13
    if constexpr (std::is_same_v<cov, glmmr::ar1Covariance>){
      if(trace>0) Rcpp::Rcout << "\nRho: " << covariance.rho;
    }
#endif
    if(trace > 0)Rcpp::Rcout << "\nu: " << u_.topRows(5).rowwise().mean().transpose();
    if(trace > 0)Rcpp::Rcout << "\nLog-likelihood: " << ll_beta << " | " << ll_theta << std::endl;
    if (iter > 2) {
      converged = check_convergence(tol,hist,iter,k0);
      if (converged && trace > 0) Rcpp::Rcout << "\nCONVERGED!";
    }
    iter++;
  }
}

// [[Rcpp::export]]
SEXP regionModel__new(const std::string& formula_,
                      const Eigen::ArrayXXd& data_,
                      const std::vector<std::string>& colnames_,
                      const Eigen::MatrixXd& X_,
                      const Eigen::ArrayXd& y_,
                      const int niter_,
                      int type = 0)
{
  using namespace rts;
  using namespace Rcpp;
  
  if(type == 0) {
    XPtr<regionModel<glmmr::Covariance>> ptr(
        new regionModel<glmmr::Covariance>(formula_, data_, colnames_, X_, y_, niter_), 
        true
    );
    return ptr;
  } else {
    // AR1 version needs T parameter - you'll need to add this to the function signature
    Rcpp::stop("AR1 version requires T parameter - use regionModel_ar__new instead");
  }
  
  return R_NilValue;
}

// Separate constructor for AR1 to avoid cluttering the signature
// [[Rcpp::export]]
SEXP regionModel_ar__new(const std::string& formula_,
                         const Eigen::ArrayXXd& data_,
                         const std::vector<std::string>& colnames_,
                         const Eigen::MatrixXd& X_,
                         const Eigen::ArrayXd& y_,
                         const int niter_,
                         const int T_)
{
  using namespace rts;
  using namespace Rcpp;

#ifdef GLMMR13
  XPtr<regionModel<glmmr::ar1Covariance>> ptr(
      new regionModel<glmmr::ar1Covariance>(formula_, data_, colnames_, X_, y_, niter_, T_), 
      true
  );
  return ptr;
#else
  throw std::runtime_error("Update to glmmrBase >= 1.2.0 for ar1");
  return 0;
#endif
  
}

// [[Rcpp::export]]
void regionModel__set_weights(SEXP xp, 
                              Rcpp::IntegerVector i, 
                              Rcpp::IntegerVector p, 
                              Rcpp::NumericVector x, 
                              int nrow, 
                              int ncol, 
                              int type = 0) {
  
  // Build triplet list from CSC components
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(x.size());
  
  for(int col = 0; col < ncol; col++) {
    for(int idx = p[col]; idx < p[col + 1]; idx++) {
      triplets.emplace_back(i[idx], col, x[idx]);
    }
  }
  
  Eigen::SparseMatrix<double> W(nrow, ncol);
  W.setFromTriplets(triplets.begin(), triplets.end());
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    ptr->weights = W;
    ptr->init_beta();
  } else {
#ifdef GLMMR13
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    ptr->weights = W;
    ptr->init_beta();
#endif
  }
  
}

// [[Rcpp::export]]
void regionModel__set_theta(SEXP xp, std::vector<double> theta, int type = 0) {
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    ptr->covariance.update_parameters_extern(theta);
  } else {
#ifdef GLMMR13
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    double rho = theta.back();
    theta.pop_back();
    ptr->covariance.update_parameters_extern(theta);
    ptr->covariance.update_rho(rho);
#else 
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    ptr->covariance.update_parameters(theta);
#endif
  }
}

// [[Rcpp::export]]
void regionModel__set_offset(SEXP xp, Eigen::VectorXd offset, int type = 0) {
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    ptr->set_offset(offset);
  } else {
#ifdef GLMMR13
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    ptr->set_offset(offset);
#endif
  }
}

// [[Rcpp::export]]
void regionModel__fit(SEXP xp, double tol, int max_iter, int hist, int k0, int type = 0){
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    ptr->fit(tol, max_iter, hist, k0);
  } else {
#ifdef GLMMR13
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    ptr->fit(tol, max_iter, hist, k0);
#endif
  }
}

// [[Rcpp::export]]
SEXP regionModel__information_matrix(SEXP xp, int type = 0){
#ifdef GLMMR13
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->information_matrix());
  } else {
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->information_matrix());
  }
#else
  Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
  return Rcpp::wrap(ptr->information_matrix());
#endif
}

// [[Rcpp::export]]
SEXP regionModel__information_matrix_theta(SEXP xp, int type = 0){
#ifdef GLMMR13
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->covariance.infomat_theta);
  } else {
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->covariance.infomat_theta);
  }
#else
  return Rcpp::wrap(Eigen::MatrixXd::Identity(2,2));
#endif
}

// [[Rcpp::export]]
SEXP regionModel__u(SEXP xp, int type = 0){
#ifdef GLMMR13
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->u());
  } 
  else {
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->u());
  }
#else
  Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
  return Rcpp::wrap(ptr->u());
#endif
}

// [[Rcpp::export]]
SEXP regionModel__get_beta(SEXP xp, int type = 0){
#ifdef GLMMR13
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->beta);
  } else {
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->beta);
  }
#else
  Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
  return Rcpp::wrap(ptr->beta);
#endif
}

// [[Rcpp::export]]
SEXP regionModel__get_weights(SEXP xp, int type = 0){
#ifdef GLMMR13
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->sampling_weights());
  } else {
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->sampling_weights());
  }
#else
  Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
  return Rcpp::wrap(ptr->sampling_weights());
#endif
}

// [[Rcpp::export]]
SEXP regionModel__log_likelihood(SEXP xp, int type = 0){
#ifdef GLMMR13
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->total_log_likelihood());
  } else {
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->total_log_likelihood());
  }
#else
  Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
  return Rcpp::wrap(ptr->total_log_likelihood());
#endif
}

// [[Rcpp::export]]
SEXP regionModel__get_theta(SEXP xp, int type = 0){
#ifdef GLMMR13
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->covariance.parameters_);
  } else {
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    std::vector<double> theta = ptr->covariance.parameters_;
    theta.push_back(ptr->covariance.rho);
    return Rcpp::wrap(theta);
  }
#else
  Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
  return Rcpp::wrap(ptr->covariance.parameters_);
#endif
}

// [[Rcpp::export]]
double max_dist(const Eigen::ArrayXXd &x){
  // this is a brute force algorithm for max distance
  // it can be improved by finding convex hull and then using rotating calipers method
  // but I haven't had the time to implement that!
  int n = x.rows();
  double maxdist = 0;
  double dist = 0;
  for(int i = 1; i < n; i++){
    for(int j = 0; j<(i-1); j++){
      dist = sqrt((x(i,0) - x(j,0))*(x(i,0) - x(j,0))+(x(i,1) - x(j,1))*(x(i,1) - x(j,1)));
      if(dist > maxdist) maxdist = dist;
    }
  }
  return maxdist;
}

// [[Rcpp::export]]
Eigen::ArrayXXd semivariogram(const Eigen::ArrayXXd &x,
                              const Eigen::ArrayXd &offs,
                              const Eigen::ArrayXd &y,
                              int nbins,
                              int nT){
  double maxd = max_dist(x);
  int n = y.size()/nT;
  Eigen::ArrayXd denoms = Eigen::ArrayXd::Zero(nbins);
  Eigen::ArrayXd sums = Eigen::ArrayXd::Zero(nbins);
  double binw = maxd/nbins;
  Eigen::ArrayXd z(y);
  for(int t = 0; t< nT; t++){
    z.segment(t*n,n) *= offs.inverse();
  }
  // Eigen::ArrayXd z = offs.inverse();
  // z *= y;
  double dist;
  // int n = x.rows();
  int binc;
  for(int t = 0; t < nT; t++){
    for(int i = 1; i < n; i++){
      for(int j = 0; j<(i-1); j++){
        dist = sqrt((x(i,0) - x(j,0))*(x(i,0) - x(j,0))+(x(i,1) - x(j,1))*(x(i,1) - x(j,1)));
        binc = static_cast<int>(std::floor(dist/binw));
        denoms(binc) += offs(i)*offs(j);
        sums(binc) += offs(i)*offs(j)*(z(i+n*t)-z(j+n*t))*(z(i+n*t)-z(j+n*t));
      }
    }
  }
  
  denoms *= 2;
  Eigen::ArrayXXd result(nbins,2);
  for(int i=0; i<nbins; i++)result(i,0) = i*binw + binw/2;
  result.col(1) = denoms.inverse()*sums;
  return result;
}

// [[Rcpp::export]]
SEXP GridData__NN(SEXP ptr_){
  Rcpp::XPtr<glmmr::griddata> ptr(ptr_);
  Eigen::ArrayXXi nn = ptr->NN;
  return Rcpp::wrap(nn);
}

// [[Rcpp::export]]
void GridData__gen_NN(SEXP ptr_, SEXP m_){
  Rcpp::XPtr<glmmr::griddata> ptr(ptr_);
  int m = Rcpp::as<int>(m_);
  ptr->genNN(m);
}

// [[Rcpp::export]]
SEXP GridData__new(SEXP x_, SEXP t_){
  Eigen::ArrayXXd x = Rcpp::as<Eigen::ArrayXXd>(x_);
  int t = Rcpp::as<int>(t_);
  Rcpp::XPtr<glmmr::griddata> ptr(new glmmr::griddata(x,t), true);
  return ptr;
}

