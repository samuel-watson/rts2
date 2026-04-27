#include <glmmr/general.h>
#include <glmmr/griddata.hpp>
#include "region_model.h"

template<typename cov>
MatrixXd rts::regionModel<cov>::u(bool scaled) const
{
  if(scaled ){
    return scaled_u_;
  } else { 
    return u_;
  }
  
}

template<typename cov>
ArrayXd rts::regionModel<cov>::sampling_weights() const 
{
  return u_weight_;
}

template<typename cov>
VectorXd rts::regionModel<cov>::zu_var() const 
{
  return zu_var_;
}

template<>
inline MatrixXd rts::regionModel<glmmr::hsgpCovariance>::information_matrix(bool monte_carlo)
{
  const int P_        = beta.size();
  const int Mdim      = covariance.Q();
  const int n_regions = weights.rows();
  const int n_cells   = weights.cols();
  const int total     = scaled_u_.rows();
  const int T         = total / n_cells;
  const int K         = scaled_u_.cols();
  
  MatrixXd Phi     = covariance.ZPhi();
  VectorXd xb_full = X * beta + offset;
  
  // ── Poisson marginalisation correction (Jensen) ──
  // E_{u|y}[exp(eta)] ≈ exp(eta_bar + 0.5 sigma2_i),
  // sigma2_i = (Phi diag(Lambda) Phi^T)_ii
  ArrayXd Lambda = covariance.LambdaSPD();
  MatrixXd PhiL = (Phi.array().rowwise() * Lambda.transpose().sqrt()).matrix();
  ArrayXd sigma2_i = PhiL.array().square().rowwise().sum();
  
  double wsum = 0.0;
  for(int k = 0; k < K; ++k) wsum += u_weight_(k);
  
  if (!monte_carlo) {
    // -------------------------------------------------------------
    // Plug-in (Laplace / GLS) information at posterior mean of u
    // -------------------------------------------------------------
    VectorXd u_bar = VectorXd::Zero(total);
    for(int k = 0; k < K; ++k) u_bar += u_weight_(k) * scaled_u_.col(k);
    u_bar /= wsum;
    u_bar.array() -= u_bar.mean();
    
    ArrayXd eta_s    = xb_full.array() + u_bar.array() - 0.5 * sigma2_i;;
    ArrayXd lambda_s = eta_s.exp();
    
    VectorXd lambda_r(n_regions * T);
    for(int t = 0; t < T; ++t){
      lambda_r.segment(t*n_regions, n_regions) =
        weights * lambda_s.segment(t*n_cells, n_cells).matrix();
    }
    
    MatrixXd tilde_X(n_regions * T, P_);
    MatrixXd tilde_A(n_regions * T, Mdim);
    for(int t = 0; t < T; ++t){
      ArrayXd  ls_t  = lambda_s.segment(t*n_cells, n_cells);
      MatrixXd X_t   = X.middleRows(t*n_cells, n_cells);
      MatrixXd Phi_t = Phi.middleRows(t*n_cells, n_cells);
      MatrixXd lsX   = (X_t.array().colwise()   * ls_t).matrix();
      MatrixXd lsPhi = (Phi_t.array().colwise() * ls_t).matrix();
      MatrixXd BX    = weights * lsX;
      MatrixXd BA    = weights * lsPhi;
      ArrayXd  inv_lr = lambda_r.segment(t*n_regions, n_regions).array().inverse();
      tilde_X.middleRows(t*n_regions, n_regions) = (BX.array().colwise() * inv_lr).matrix();
      tilde_A.middleRows(t*n_regions, n_regions) = (BA.array().colwise() * inv_lr).matrix();
    }
    
    VectorXd Wr = lambda_r;
    MatrixXd WrtildeX = (tilde_X.array().colwise() * Wr.array()).matrix();
    MatrixXd WrtildeA = (tilde_A.array().colwise() * Wr.array()).matrix();
    
    MatrixXd XtWX = tilde_X.transpose() * WrtildeX;
    MatrixXd AtWA = tilde_A.transpose() * WrtildeA;
    MatrixXd XtWA = tilde_X.transpose() * WrtildeA;
    
    ArrayXd inv_lambda = 1.0 / covariance.LambdaSPD();
    AtWA.diagonal() += inv_lambda.matrix();
    
    Eigen::LLT<MatrixXd> llt(AtWA);
    MatrixXd solve_XtWA = llt.solve(XtWA.transpose());
    return XtWX - XtWA * solve_XtWA;
  }
  
  // -------------------------------------------------------------
  // Louis-identity Monte Carlo information (β-block):
  //   J_ββ = Σ_k w_k  tilde_X_k' diag(λ_r_k) tilde_X_k
  //        - Σ_k w_k (g_k - ḡ)(g_k - ḡ)'
  //   g_k  = tilde_X_k' (y - λ_r_k)
  // Per-sample centering of Zu^(k) matches nr_beta.
  // -------------------------------------------------------------
  
  MatrixXd bread = MatrixXd::Zero(P_, P_);
  MatrixXd g_mat(P_, K);
  VectorXd g_bar = VectorXd::Zero(P_);
  
  for(int k = 0; k < K; ++k){
    VectorXd zu_k = scaled_u_.col(k);
    zu_k.array() -= zu_k.mean();                // per-sample centering, matches nr_beta
    
    ArrayXd eta_sk    = xb_full.array() + zu_k.array() - 0.5 * sigma2_i;;
    ArrayXd lambda_sk = eta_sk.exp(); 
    
    VectorXd lambda_rk(n_regions * T);
    for(int t = 0; t < T; ++t){
      lambda_rk.segment(t*n_regions, n_regions) =
        weights * lambda_sk.segment(t*n_cells, n_cells).matrix();
    }
    
    MatrixXd tilde_Xk(n_regions * T, P_);
    for(int t = 0; t < T; ++t){
      ArrayXd  ls_t = lambda_sk.segment(t*n_cells, n_cells);
      MatrixXd X_t  = X.middleRows(t*n_cells, n_cells);
      MatrixXd lsX  = (X_t.array().colwise() * ls_t).matrix();
      MatrixXd BX   = weights * lsX;
      ArrayXd  inv_lr = lambda_rk.segment(t*n_regions, n_regions).array().inverse();
      tilde_Xk.middleRows(t*n_regions, n_regions) = (BX.array().colwise() * inv_lr).matrix();
    }
    
    double wk = u_weight_(k) / wsum;
    
    // bread_k: tilde_Xk' diag(λ_rk) tilde_Xk
    MatrixXd WtXk = (tilde_Xk.array().colwise() * lambda_rk.array()).matrix();
    bread.noalias() += wk * (tilde_Xk.transpose() * WtXk);
    
    // score: g_k = tilde_Xk' (y - λ_rk)
    VectorXd gk = tilde_Xk.transpose() * (y.matrix() - lambda_rk);
    g_mat.col(k) = gk;
    g_bar.noalias() += wk * gk;
  }
  
  MatrixXd meat = MatrixXd::Zero(P_, P_);
  for(int k = 0; k < K; ++k){
    double wk  = u_weight_(k) / wsum;
    VectorXd d = g_mat.col(k) - g_bar;
    meat.noalias() += wk * (d * d.transpose());
  }
  
  return bread - meat;
}

template<>
inline MatrixXd rts::regionModel<glmmr::spdeCovariance>::information_matrix(bool monte_carlo)
{
  const int P_        = beta.size();
  const int nv        = covariance.Q();
  const int n_regions = weights.rows();
  const int n_cells   = weights.cols();
  const int K         = scaled_u_.cols();
  
  const SparseMatrix<double>& ZA = covariance.ZA_;   // n_cells × nv
  
  // Posterior mean of cell-level effects
  double wsum = 0.0;
  for(int k = 0; k < K; ++k) wsum += u_weight_(k);
  
  VectorXd u_bar_cell = VectorXd::Zero(n_cells);
  for(int k = 0; k < K; ++k) u_bar_cell += u_weight_(k) * scaled_u_.col(k);
  u_bar_cell /= wsum;
  u_bar_cell.array() -= u_bar_cell.mean();
  
  // ── Poisson marginalisation correction (Jensen) ──
  // E_{u|y}[exp(eta)] ≈ exp(eta_bar + 0.5 sigma^2_i) where
  // sigma^2_i = (ZA Q^{-1} ZA^T)_ii  approximated via Hutchinson on the prior
  ArrayXd  sigma2_i = covariance.diag_ZA_Qinv_ZAt_hutch();
  VectorXd xb_full  = X * beta + offset;
  ArrayXd  eta_s    = xb_full.array() + u_bar_cell.array() - 0.5 * sigma2_i;
  ArrayXd  lambda_s = eta_s.exp();
  
  if (!monte_carlo) {
    // -------------------------------------------------------------
    // Plug-in (Laplace / GLS) marginal information at corrected mean
    //   M_β = X̃' diag(λ_r) X̃  -  X̃' diag(λ_r) Ã (Ã' diag(λ_r) Ã + Q)^{-1} Ã' diag(λ_r) X̃
    // with X̃ = diag(1/λ_r) * weights * diag(λ_s) * X,
    //      Ã = diag(1/λ_r) * weights * diag(λ_s) * ZA
    // -------------------------------------------------------------
    VectorXd lambda_r = weights * lambda_s.matrix();
    ArrayXd  inv_lr   = lambda_r.array().inverse();
    
    MatrixXd lsX = (X.array().colwise() * lambda_s).matrix();
    MatrixXd BX  = weights * lsX;
    MatrixXd tilde_X = (BX.array().colwise() * inv_lr).matrix();          // n_regions × P
    
    SparseMatrix<double> Lam_s_ZA = lambda_s.matrix().asDiagonal() * ZA;
    SparseMatrix<double> BA_sp    = weights * Lam_s_ZA;                   // n_regions × nv (sparse)
    SparseMatrix<double> tilde_A  = inv_lr.matrix().asDiagonal() * BA_sp; // n_regions × nv (sparse)
    
    VectorXd Wr = lambda_r;
    MatrixXd WrtildeX = (tilde_X.array().colwise() * Wr.array()).matrix();
    MatrixXd XtWX     = tilde_X.transpose() * WrtildeX;                   // P × P
    
    // X̃' diag(Wr) Ã   (P × nv, kept dense since P is small)
    MatrixXd WtX_dense = WrtildeX;                                         // n_regions × P
    MatrixXd XtWA      = MatrixXd(tilde_A.transpose() * WtX_dense).transpose();  // P × nv
    
    // Posterior precision  P_post = Ã' diag(Wr) Ã + Q   (sparse)
    SparseMatrix<double> Wr_tildeA = Wr.asDiagonal() * tilde_A;
    SparseMatrix<double> AtWA      = SparseMatrix<double>(tilde_A.transpose()) * Wr_tildeA;
    SparseMatrix<double> P_post    = AtWA + covariance.Q_mat;
    P_post.makeCompressed();
    
    // Reuse cached chol_P (refactor at this W)
    covariance.refactor_P(Wr);
    MatrixXd solve_XtWA_t = covariance.chol_P.solve(MatrixXd(XtWA.transpose()));  // nv × P
    
    return XtWX - XtWA * solve_XtWA_t;
  }
  
  // -------------------------------------------------------------
  // Louis-identity Monte Carlo information (β-block)
  // J_ββ = Σ_k w_k  X̃_k' diag(λ_r_k) X̃_k  -  Σ_k w_k (g_k - ḡ)(g_k - ḡ)'
  // g_k  = X̃_k' (y - λ_r_k)
  // -------------------------------------------------------------
  MatrixXd bread = MatrixXd::Zero(P_, P_);
  MatrixXd g_mat(P_, K);
  VectorXd g_bar = VectorXd::Zero(P_);
  
  for(int k = 0; k < K; ++k){
    VectorXd zu_k = scaled_u_.col(k);
    zu_k.array() -= zu_k.mean();                                     // per-sample centering, matches nr_beta
    
    ArrayXd eta_sk    = xb_full.array() + zu_k.array() - 0.5 * sigma2_i;
    ArrayXd lambda_sk = eta_sk.exp();
    VectorXd lambda_rk = weights * lambda_sk.matrix();
    ArrayXd  inv_lr_k  = lambda_rk.array().inverse();
    
    MatrixXd lsX_k = (X.array().colwise() * lambda_sk).matrix();
    MatrixXd BX_k  = weights * lsX_k;
    MatrixXd tilde_Xk = (BX_k.array().colwise() * inv_lr_k).matrix();
    
    double wk = u_weight_(k) / wsum;
    
    MatrixXd WtXk = (tilde_Xk.array().colwise() * lambda_rk.array()).matrix();
    bread.noalias() += wk * (tilde_Xk.transpose() * WtXk);
    
    VectorXd gk = tilde_Xk.transpose() * (y.matrix() - lambda_rk);
    g_mat.col(k) = gk;
    g_bar.noalias() += wk * gk;
  }
  
  MatrixXd meat = MatrixXd::Zero(P_, P_);
  for(int k = 0; k < K; ++k){
    double wk  = u_weight_(k) / wsum;
    VectorXd d = g_mat.col(k) - g_bar;
    meat.noalias() += wk * (d * d.transpose());
  }
  
  return bread - meat;
}

template<typename cov>
MatrixXd rts::regionModel<cov>::information_matrix(bool monte_carlo)
{
  const int P_        = beta.size();
  const int n_regions = weights.rows();
  const int n_cells   = weights.cols();
  
  int T   = 1;
  int n_A = n_cells;
#ifdef GLMMR13
  if constexpr (std::is_same_v<cov, glmmr::ar1Covariance>) {
    T   = covariance.Q() / covariance.Covariance::Q();
    n_A = covariance.Covariance::Q();
  }
#endif
  
  const int K    = scaled_u_.cols();
  const int n_rT = n_regions * T;
  
  VectorXd xb_full = X * beta + offset;
  MatrixXd ZL      = covariance.ZL();
  const int Q_     = ZL.cols();
  
  double wsum = 0.0;
  for(int k = 0; k < K; ++k) wsum += u_weight_(k);
  
  if (!monte_carlo) {
    // ---- Laplace / GLS plug-in at posterior mean of u -----------------
    VectorXd u_bar = VectorXd::Zero(scaled_u_.rows());
    for(int k = 0; k < K; ++k) u_bar += u_weight_(k) * scaled_u_.col(k);
    u_bar /= wsum;
    u_bar.array() -= u_bar.mean();
    
    ArrayXd eta_s    = xb_full.array() + u_bar.array();
    ArrayXd lambda_s = eta_s.exp();
    
    VectorXd lambda_r(n_rT);
    for(int t = 0; t < T; ++t){
      lambda_r.segment(t*n_regions, n_regions)
      = weights * lambda_s.segment(t*n_A, n_A).matrix();
    }
    
    MatrixXd tilde_X(n_rT, P_);
    MatrixXd tilde_A(n_rT, Q_);
    for(int t = 0; t < T; ++t){
      ArrayXd  ls_t = lambda_s.segment(t*n_A, n_A);
      MatrixXd X_t  = X.middleRows(t*n_A, n_A);
      MatrixXd ZL_t = ZL.middleRows(t*n_A, n_A);
      MatrixXd lsX  = (X_t.array().colwise()  * ls_t).matrix();
      MatrixXd lsZL = (ZL_t.array().colwise() * ls_t).matrix();
      MatrixXd BX   = weights * lsX;
      MatrixXd BA   = weights * lsZL;
      ArrayXd  inv_lr = lambda_r.segment(t*n_regions, n_regions).array().inverse();
      tilde_X.middleRows(t*n_regions, n_regions) = (BX.array().colwise() * inv_lr).matrix();
      tilde_A.middleRows(t*n_regions, n_regions) = (BA.array().colwise() * inv_lr).matrix();
    }
    
    VectorXd Wr = lambda_r;
    MatrixXd WrtildeX = (tilde_X.array().colwise() * Wr.array()).matrix();
    MatrixXd WrtildeA = (tilde_A.array().colwise() * Wr.array()).matrix();
    MatrixXd XtWX = tilde_X.transpose() * WrtildeX;
    MatrixXd AtWA = tilde_A.transpose() * WrtildeA;
    MatrixXd XtWA = tilde_X.transpose() * WrtildeA;
    AtWA.diagonal().array() += 1.0;
    
    if (Q_ > n_rT) {
      ArrayXd sqrt_Wr = Wr.array().sqrt();
      MatrixXd C = (tilde_A.array().colwise() * sqrt_Wr).matrix();
      MatrixXd CCt = C * C.transpose();
      CCt.diagonal().array() += 1.0;
      Eigen::LLT<MatrixXd> llt_CCt(CCt);
      MatrixXd Ct_XtWAt = C * XtWA.transpose();
      MatrixXd sol      = llt_CCt.solve(Ct_XtWAt);
      return XtWX - XtWA * XtWA.transpose() + Ct_XtWAt.transpose() * sol;
    } else {
      Eigen::LLT<MatrixXd> llt(AtWA);
      MatrixXd solve_XtWA = llt.solve(XtWA.transpose());
      return XtWX - XtWA * solve_XtWA;
    }
  }
  
  // ---- Louis-identity Monte Carlo information (β-block) ---------------
  //   J_ββ = Σ_k w_k tilde_X_k' diag(λ_r_k) tilde_X_k
  //        - Σ_k w_k (g_k - ḡ)(g_k - ḡ)'
  //   g_k  = tilde_X_k' (y - λ_r_k)
  // Per-sample centering of ZL·v^(k) matches nr_beta.
  // ---------------------------------------------------------------------
  
  MatrixXd bread = MatrixXd::Zero(P_, P_);
  MatrixXd g_mat(P_, K);
  VectorXd g_bar = VectorXd::Zero(P_);
  
  for(int k = 0; k < K; ++k){
    VectorXd zu_k = scaled_u_.col(k);    // = ZL · v^(k), length n_A*T
    zu_k.array() -= zu_k.mean();         // per-sample centering
    
    ArrayXd eta_sk    = xb_full.array() + zu_k.array();
    ArrayXd lambda_sk = eta_sk.exp();
    
    VectorXd lambda_rk(n_rT);
    for(int t = 0; t < T; ++t){
      lambda_rk.segment(t*n_regions, n_regions)
      = weights * lambda_sk.segment(t*n_A, n_A).matrix();
    }
    
    MatrixXd tilde_Xk(n_rT, P_);
    for(int t = 0; t < T; ++t){
      ArrayXd  ls_t = lambda_sk.segment(t*n_A, n_A);
      MatrixXd X_t  = X.middleRows(t*n_A, n_A);
      MatrixXd lsX  = (X_t.array().colwise() * ls_t).matrix();
      MatrixXd BX   = weights * lsX;
      ArrayXd  inv_lr = lambda_rk.segment(t*n_regions, n_regions).array().inverse();
      tilde_Xk.middleRows(t*n_regions, n_regions) = (BX.array().colwise() * inv_lr).matrix();
    }
    
    double wk = u_weight_(k) / wsum;
    
    MatrixXd WtXk = (tilde_Xk.array().colwise() * lambda_rk.array()).matrix();
    bread.noalias() += wk * (tilde_Xk.transpose() * WtXk);
    
    VectorXd gk = tilde_Xk.transpose() * (y.matrix() - lambda_rk);
    g_mat.col(k) = gk;
    g_bar.noalias() += wk * gk;
  }
  
  MatrixXd meat = MatrixXd::Zero(P_, P_);
  for(int k = 0; k < K; ++k){
    double wk  = u_weight_(k) / wsum;
    VectorXd d = g_mat.col(k) - g_bar;
    meat.noalias() += wk * (d * d.transpose());
  }
  
  return bread - meat;
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

// ── lambda_r ──
template<>
inline MatrixXd rts::regionModel<glmmr::hsgpCovariance>::lambda_r()
{
  int n_regions = weights.rows();
  int n_cells   = weights.cols();
  int total     = scaled_u_.rows();     // n_cells * T
  int T         = total / n_cells;
  
  MatrixXd lr(n_regions * T, u_.cols());
  for(int i = 0; i < u_.cols(); i++){
    ArrayXd eta_s   = (X * beta + offset).array() + scaled_u_.col(i).array();
    ArrayXd lambda_s = eta_s.exp();
    for(int t = 0; t < T; t++){
      VectorXd lambda_s_t = lambda_s.segment(t * n_cells, n_cells).matrix();
      lr.col(i).segment(t * n_regions, n_regions) = weights * lambda_s_t;
    }
  }
  return lr;
}

template<>
inline MatrixXd rts::regionModel<glmmr::spdeCovariance>::lambda_r()
{
  MatrixXd lr(weights.rows(), u_.cols());
  for(int i = 0; i < u_.cols(); i++){
    ArrayXd eta_s    = (X * beta + offset).array() + scaled_u_.col(i).array();
    ArrayXd lambda_s = eta_s.exp();
    lr.col(i) = weights * lambda_s.matrix();
  }
  return lr;
}

// ── nr_beta ──
template<>
inline void rts::regionModel<glmmr::hsgpCovariance>::nr_beta()
{
  int niter_    = u_.cols();
  int n_regions = weights.rows();
  int n_cells   = weights.cols();
  int total     = scaled_u_.rows();
  int T         = total / n_cells;
  
  MatrixXd zd(X.rows(), niter_);
  for(int i = 0; i < niter_; i++){
    zd.col(i) = (X * beta + offset).array() + scaled_u_.col(i).array();
  }
  
  if(niter_ > 1){
    for(int i = 0; i < niter_; i++){
      double zu_mean = scaled_u_.col(i).mean();
      zd.col(i).array() -= zu_mean;
    }
  }
  
  
  MatrixXd XtWXm = MatrixXd::Zero(beta.size(), beta.size());
  VectorXd score = VectorXd::Zero(beta.size());
  
  for(int i = 0; i < niter_; ++i){
    ArrayXd lambda_s = zd.col(i).array().exp();
    
    VectorXd lambda_r(n_regions * T);
    for(int t = 0; t < T; t++){
      VectorXd ls_t = lambda_s.segment(t * n_cells, n_cells).matrix();
      lambda_r.segment(t * n_regions, n_regions) = weights * ls_t;
    }
    
    VectorXd resid_r = (y / lambda_r.array() - 1.0).matrix();
    
    // B = block-diagonal application of weights * diag(lambda_s_t) * X_t
    MatrixXd B(n_regions * T, beta.size());
    for(int t = 0; t < T; t++){
      ArrayXd ls_t = lambda_s.segment(t * n_cells, n_cells);
      MatrixXd X_t = X.middleRows(t * n_cells, n_cells);
      MatrixXd lsX = (X_t.array().colwise() * ls_t).matrix();
      B.middleRows(t * n_regions, n_regions) = weights * lsX;
    }
    
    score.noalias() += u_weight_(i) * B.transpose() * resid_r;
    
    ArrayXd inv_lr = lambda_r.array().inverse();
    MatrixXd inv_lr_B = (B.array().colwise() * inv_lr).matrix();
    XtWXm.noalias() += u_weight_(i) * B.transpose() * inv_lr_B;
  }
  
  M = XtWXm;
  Eigen::LLT<MatrixXd> llt(XtWXm);
  gradients.head(X.cols()) = score;
  VectorXd bincr = llt.solve(score);
  beta += bincr;
  ll_beta = log_likelihood().mean();
  
  if(beta.hasNaN()){
    Rcpp::Rcout << "DIAGNOSTIC: NaN in beta" << std::endl;
    Rcpp::stop("nr_beta (HSGP): NaN detected");
  }
  if(std::isnan(ll_beta) || std::isinf(ll_beta)){
    Rcpp::Rcout << "DIAGNOSTIC: NaN/Inf in ll_beta" << std::endl;
    Rcpp::Rcout << "lambda_r range: [" << lambda_r().minCoeff() 
                << ", " << lambda_r().maxCoeff() << "]" << std::endl;
    Rcpp::Rcout << "lambda_r has zeros: " 
                << ((lambda_r().array() <= 0).any() ? "yes" : "no") << std::endl;
    Rcpp::Rcout << "beta: " << beta.transpose() << std::endl;
    Rcpp::Rcout << "y range: [" << y.minCoeff() << ", " << y.maxCoeff() << "]" << std::endl;
    Rcpp::stop("nr_beta (HSGP): NaN/Inf in ll_beta");
  }
}

template<>
inline void rts::regionModel<glmmr::spdeCovariance>::nr_beta()
{
  int niter_    = u_.cols();
  int n_regions = weights.rows();
  int n_cells   = weights.cols();
  
  MatrixXd zd(X.rows(), niter_);
  for(int i = 0; i < niter_; i++){
    zd.col(i) = (X * beta + offset).array() + scaled_u_.col(i).array();
  }
  
  if(niter_ > 1){
    for(int i = 0; i < niter_; i++){
      double zu_mean = scaled_u_.col(i).mean();
      zd.col(i).array() -= zu_mean;
    }
  }
  
  MatrixXd XtWXm = MatrixXd::Zero(beta.size(), beta.size());
  VectorXd score = VectorXd::Zero(beta.size());
  
  for(int i = 0; i < niter_; ++i){
    ArrayXd  lambda_s = zd.col(i).array().exp();
    VectorXd lambda_r = weights * lambda_s.matrix();
    VectorXd resid_r  = (y / lambda_r.array() - 1.0).matrix();
    
    // B = weights * diag(lambda_s) * X
    MatrixXd lsX = (X.array().colwise() * lambda_s).matrix();
    MatrixXd B   = weights * lsX;
    
    score.noalias() += u_weight_(i) * B.transpose() * resid_r;
    
    ArrayXd  inv_lr   = lambda_r.array().inverse();
    MatrixXd inv_lr_B = (B.array().colwise() * inv_lr).matrix();
    XtWXm.noalias()  += u_weight_(i) * B.transpose() * inv_lr_B;
  }
  
  M = XtWXm;
  Eigen::LLT<MatrixXd> llt(XtWXm);
  gradients.head(X.cols()) = score;
  VectorXd bincr = llt.solve(score);
  beta += bincr;
  ll_beta = log_likelihood().mean();
  
  if(beta.hasNaN()){
    Rcpp::Rcout << "DIAGNOSTIC: NaN in beta" << std::endl;
    Rcpp::stop("nr_beta (SPDE): NaN detected");
  }
  if(std::isnan(ll_beta) || std::isinf(ll_beta)){
    Rcpp::Rcout << "DIAGNOSTIC: NaN/Inf in ll_beta" << std::endl;
    Rcpp::Rcout << "lambda_r range: [" << lambda_r().minCoeff()
                << ", " << lambda_r().maxCoeff() << "]" << std::endl;
    Rcpp::Rcout << "lambda_r has zeros: "
                << ((lambda_r().array() <= 0).any() ? "yes" : "no") << std::endl;
    Rcpp::Rcout << "beta: " << beta.transpose() << std::endl;
    Rcpp::Rcout << "y range: [" << y.minCoeff() << ", " << y.maxCoeff() << "]" << std::endl;
    Rcpp::stop("nr_beta (SPDE): NaN/Inf in ll_beta");
  }
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
  // Use row sums rather than total sum for a per-region average weight
  VectorXd row_sums = VectorXd::Zero(weights.rows());
  for(int k = 0; k < weights.outerSize(); ++k)
    for(SparseMatrix<double>::InnerIterator it(weights, k); it; ++it)
      row_sums(it.row()) += it.value();
  
  double mean_y = y.mean();
  double mean_row_sum = row_sums.mean();
  double mean_offset = offset.mean();
  
  if(X.cols() > 0){
    double ratio = mean_y / std::max(mean_row_sum, 1e-10);
    beta(0) = std::log(std::max(ratio, 1e-6)) - mean_offset;
    // Clamp to prevent extreme starting values
    beta(0) = std::max(std::min(beta(0), 10.0), -10.0);
  }
}

template<>
inline void rts::regionModel<glmmr::hsgpCovariance>::usample(const int niter_)
{
  const int n_regions = weights.rows();
  const int n_cells   = weights.cols();
  const int Mspec     = covariance.Q();
  const int total = scaled_u_.rows();  // n_cells * T
  const int T     = total / n_cells;
  
  MatrixXd ZPhi   = covariance.ZPhi();             // n_cells × M
  ArrayXd  Lambda = covariance.LambdaSPD();
  ArrayXd  sqrtLam  = Lambda.sqrt();
  ArrayXd  inv_lam  = 1.0 / Lambda;
  
  ArrayXd xb = (X * beta + offset).array();
  
  // ── IRLS for posterior mode in u-space ──
  // Precision:  P = C^T C + diag(1/Lambda)
  // Woodbury via C_tilde = C * diag(sqrt(Lambda)):
  //   P = S^{-1} (C_tilde^T C_tilde + I) S^{-1},  S = diag(sqrt(Lambda))
  //   P^{-1} x = S * (C_tilde^T C_tilde + I)^{-1} * S * x
  
  VectorXd b = u_mean_;
  if(b.size() != Mspec){ b.resize(Mspec); b.setZero(); }
  VectorXd bnew(Mspec);
  double diff = 1.0;
  int itero = 0;
  
  MatrixXd C_tilde(n_regions, Mspec);
  LLT<MatrixXd> llt_CCt;
  while(diff > 1e-6 && itero < 20){
    ArrayXd eta_s    = xb + (ZPhi * b).array();
    ArrayXd lambda_s = eta_s.exp();
    int T = total / n_cells;
    
    VectorXd lambda_r(n_regions * T);
    for(int t = 0; t < T; t++){
      VectorXd ls_t = lambda_s.segment(t * n_cells, n_cells).matrix();
      lambda_r.segment(t * n_regions, n_regions) = weights * ls_t;
    }
    VectorXd resid_r = (y / lambda_r.array() - 1.0).matrix();
    
    // B = weights * diag(lambda_s) * ZPhi, per time block
    MatrixXd A_scaled = (ZPhi.array().colwise() * lambda_s).matrix();
    MatrixXd B(n_regions * T, Mspec);
    for(int t = 0; t < T; t++){
      MatrixXd A_t = A_scaled.middleRows(t * n_cells, n_cells);
      B.middleRows(t * n_regions, n_regions) = weights * A_t;
    }
    
    ArrayXd inv_sqrt_lr = lambda_r.array().sqrt().inverse();
    MatrixXd C = (B.array().colwise() * inv_sqrt_lr).matrix();
    C_tilde = C * sqrtLam.matrix().asDiagonal();
    
    // Factorise C_tilde * C_tilde^T + I   (n_regions × n_regions)
    MatrixXd CCt = C_tilde * C_tilde.transpose();
    CCt.diagonal().array() += 1.0;
    llt_CCt.compute(CCt);
    
    // RHS:  yb = C^T C b + score,  where score = B^T * resid_r
    VectorXd Cb = C * b;
    VectorXd score = B.transpose() * resid_r;
    VectorXd yb = C.transpose() * Cb + score;
    
    // Solve P^{-1} yb = S * (C_tilde^T C_tilde + I)^{-1} * S * yb
    VectorXd yb_s = sqrtLam.matrix().asDiagonal() * yb;
    VectorXd Ct_yb_s = C_tilde * yb_s;
    VectorXd tmp = llt_CCt.solve(Ct_yb_s);
    VectorXd v = yb_s - C_tilde.transpose() * tmp;
    bnew = sqrtLam.matrix().asDiagonal() * v;
    
    diff = (b - bnew).array().abs().maxCoeff();
    itero++;
    b.swap(bnew);
  }
  
  u_mean_ = b;
  double log_det_CCt_I = 2.0 * llt_CCt.matrixL().toDenseMatrix()
                          .diagonal().array().log().sum();
  log_det_P_ = log_det_CCt_I - Lambda.log().sum();   // HSGP
  
  // ── Posterior variance diagonals ──
  // Var(u) = S * (C_tilde^T C_tilde + I)^{-1} * S,   S = diag(sqrtLam)
  // Woodbury: (C_tilde^T C_tilde + I)^{-1} = I - C_tilde^T (C_tilde C_tilde^T + I)^{-1} C_tilde
  //
  // diag(Var(u))_i = Lambda_i * (1 - c_i^T M^{-1} c_i),
  //   where c_i is the i-th column of C_tilde and M = C_tilde C_tilde^T + I
  
  MatrixXd MinvC = llt_CCt.solve(C_tilde);             // n_regions × Mspec
  ArrayXd  quad  = (C_tilde.array() * MinvC.array()).colwise().sum();   // length Mspec
  u_var_diag_    = (Lambda * (1.0 - quad)).matrix();
  
  // Var(Zu) diagonal on the cell grid:
  // Let A_s = ZPhi * S, so Var(Zu) = A_s (I - C_tilde^T M^{-1} C_tilde) A_s^T
  // diag = rowwise ||A_s||^2 - rowwise (K .* (M^{-1} K^T)^T).sum(),  K = A_s C_tilde^T
  MatrixXd A_s  = (ZPhi.array().rowwise() * sqrtLam.transpose()).matrix();   // n_cells × Mspec
  ArrayXd  diag_full = A_s.array().square().rowwise().sum();                 // n_cells
  MatrixXd K    = A_s * C_tilde.transpose();                                 // n_cells × n_regions
  MatrixXd MinvKt = llt_CCt.solve(K.transpose());                            // n_regions × n_cells
  ArrayXd  diag_corr = (K.array() * MinvKt.transpose().array()).rowwise().sum();
  zu_var_ = (diag_full - diag_corr).matrix();
  
  // ── Resize storage ──
  if(u_.cols() != niter_){
    u_.resize(Mspec, niter_);
    u_solve_.resize(Mspec, niter_);
    scaled_u_.resize(n_cells * T, niter_);
    u_weight_.resize(niter_);
    u_loglik_.resize(niter_);
  }
  
  if(niter_ == 1){
    // Laplace: posterior mode only
    // v = 0  =>  u = u_mean_,  proposal log-density const (set to 0)
    u_.col(0)        = u_mean_;
    u_loglik_(0)     = 0.0;
    scaled_u_        = ZPhi * u_;
    u_weight_(0)     = 1.0;
  } else {
    // ── Sample from N(u_mean, P^{-1}) ──
    // P^{-1} = S * (C_tilde^T C_tilde + I)^{-1} * S
    // Sample: u = u_mean + S * (C_tilde^T C_tilde + I)^{-1} (z + C_tilde^T w)
    //   where z, w ~ N(0, I)
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> d(0.0, 1.0);
    
    MatrixXd znorm(Mspec, niter_);
    for(int i = 0; i < znorm.size(); ++i) znorm.data()[i] = d(gen);
    
    MatrixXd w(n_regions * T, niter_);
    for(int i = 0; i < w.size(); ++i) w.data()[i] = d(gen);
    
    // v = (C_tilde^T C_tilde + I)^{-1} (z + C_tilde^T w)   via Woodbury
    MatrixXd RHS  = znorm + C_tilde.transpose() * w;
    MatrixXd CRHS = C_tilde * RHS;
    MatrixXd stmp = llt_CCt.solve(CRHS);
    MatrixXd v_samples = RHS - C_tilde.transpose() * stmp;
    
    // Proposal log-density (up to constant):
    //   -0.5 * v^T (C_tilde^T C_tilde + I) v = -0.5 * (||v||^2 + ||C_tilde v||^2)
    MatrixXd Ct_v = C_tilde * v_samples;
    u_loglik_ = -0.5 * (v_samples.colwise().squaredNorm().array()
                          + Ct_v.colwise().squaredNorm().array()).matrix();
    
    // u = u_mean + diag(sqrt(Lambda)) * v
    u_ = (sqrtLam.matrix().asDiagonal() * v_samples).colwise() + u_mean_;
    
    // Map to cell-space
    scaled_u_ = ZPhi * u_;
    
    // ── Importance weights ──
    // Prior: log f(u|theta) = -0.5 * sum_k u_k^2 / Lambda_k  (+ const)
    u_weight_.setZero();
    ArrayXd ll = log_likelihood();
    
    for(int i = 0; i < u_.cols(); i++){
      double llmod   = ll(i);
      double llprior = -0.5 * (u_.col(i).array().square() * inv_lam).sum();
      u_weight_(i)   = llmod + llprior - u_loglik_(i);
    }
    u_weight_ -= u_weight_.maxCoeff();
    u_weight_  = u_weight_.exp();
    u_weight_ *= 1.0 / u_weight_.sum();
    
    // ── Diagnostics ──
    if(u_mean_.hasNaN()){
      Rcpp::Rcout << "DIAGNOSTIC: NaN in u_mean_" << std::endl;
      Rcpp::stop("usample (HSGP): NaN detected in u_mean_");
    }
    double ess = 1.0 / (u_weight_ * u_weight_).sum();
    if(ess < 2.0){
      Rcpp::Rcout << "WARNING: Very low ESS = " << ess
                  << " / " << u_weight_.size() << std::endl;
    }
  }
}

template<>
inline void rts::regionModel<glmmr::spdeCovariance>::usample(const int niter_)
{
  const int n_regions = weights.rows();
  const int n_cells   = weights.cols();
  const int nv        = covariance.Q();
  
  // Cached prior structure
  if(!covariance.chol_Q_current) covariance.refactor_Q();
  auto& chol_Q = covariance.chol_Q;
  const SparseMatrix<double>& ZA = covariance.ZA_;          // n_cells × nv
  
  ArrayXd xb = (X * beta + offset).array();
  
  // ── IRLS for posterior mode in u-space ──
  // P = C^T C + Q,   C = diag(1/sqrt(lambda_r)) * weights * diag(lambda_s) * ZA
  // Woodbury:  P^{-1} y = Q^{-1} y - Y K^{-1} C Q^{-1} y,
  //            Y = Q^{-1} C^T,  K = I + C Y   (n_regions × n_regions dense)
  VectorXd b = u_mean_;
  if(b.size() != nv){ b.resize(nv); b.setZero(); }
  VectorXd bnew(nv);
  double diff = 1.0;
  int itero = 0;
  
  SparseMatrix<double> C_sp;
  MatrixXd Y_dense;          // nv × n_regions  = Q^{-1} C^T
  MatrixXd K_mat;
  LLT<MatrixXd> llt_K;
  
  while(diff > 1e-6 && itero < 20){
    ArrayXd  eta_s    = xb + (ZA * b).array();
    ArrayXd  lambda_s = eta_s.exp();
    VectorXd lambda_r = weights * lambda_s.matrix();
    VectorXd resid_r  = (y / lambda_r.array() - 1.0).matrix();
    ArrayXd  inv_sqrt_lr = lambda_r.array().sqrt().inverse();
    
    // B = weights * diag(lambda_s) * ZA   (n_regions × nv, sparse)
    SparseMatrix<double> Lam_s_ZA = lambda_s.matrix().asDiagonal() * ZA;
    SparseMatrix<double> B_sp     = weights * Lam_s_ZA;
    C_sp = inv_sqrt_lr.matrix().asDiagonal() * B_sp;
    
    VectorXd score_u = B_sp.transpose() * resid_r;
    VectorXd Cb      = C_sp * b;
    VectorXd yb      = SparseMatrix<double>(C_sp.transpose()) * Cb + score_u;
    
    // Woodbury setup (refresh each IRLS step since lambda_s changes)
    Y_dense = chol_Q.solve(MatrixXd(C_sp.transpose()));      // nv × n_regions
    K_mat   = MatrixXd::Identity(n_regions, n_regions) + C_sp * Y_dense;
    llt_K.compute(K_mat);
    
    VectorXd Qinv_yb = chol_Q.solve(yb);
    VectorXd Cy_yb   = C_sp * Qinv_yb;
    bnew             = Qinv_yb - Y_dense * llt_K.solve(Cy_yb);
    
    diff = (b - bnew).array().abs().maxCoeff();
    itero++;
    b.swap(bnew);
  }
  u_mean_ = b;
  
  // log det P = log det Q + log det K  (matrix determinant lemma)
  double log_det_K = 2.0 * llt_K.matrixL().toDenseMatrix()
                      .diagonal().array().log().sum();
  double log_det_Q = -covariance.log_determinant();          // log_determinant returns log det D = -log det Q
  log_det_P_       = log_det_Q + log_det_K;
  
  // ── Posterior variance diagonals ──
  // Var(u) diag = diag(Q^{-1}) - diag(Y K^{-1} Y^T)
  u_var_diag_ = covariance.diag_Qinv_hutch().matrix();
  MatrixXd Kinv_Yt = llt_K.solve(Y_dense.transpose());
  ArrayXd  corr_u  = (Y_dense.array() * Kinv_Yt.transpose().array()).rowwise().sum();
  u_var_diag_     -= corr_u.matrix();
  
  // Var(ZA u) diag at cells  ≈  diag_ZA_Qinv_ZAt - diag( (ZA Y) K^{-1} (ZA Y)^T )
  ArrayXd  diag_prior = covariance.diag_ZA_Qinv_ZAt_hutch();
  MatrixXd ZA_Y       = ZA * Y_dense;
  MatrixXd Kinv_ZAYt  = llt_K.solve(ZA_Y.transpose());
  ArrayXd  corr_zu    = (ZA_Y.array() * Kinv_ZAYt.transpose().array()).rowwise().sum();
  zu_var_             = (diag_prior - corr_zu).matrix();
  
  // ── Resize storage ──
  if(u_.cols() != niter_){
    u_.resize(nv, niter_);
    u_solve_.resize(nv, niter_);
    scaled_u_.resize(n_cells, niter_);
    u_weight_.resize(niter_);
    u_loglik_.resize(niter_);
  }
  
  if(niter_ == 1){
    // Laplace
    u_.col(0)    = u_mean_;
    u_loglik_(0) = 0.0;
    scaled_u_    = ZA * u_;
    u_weight_(0) = 1.0;
  } else {
    // ── Sample from N(u_mean, P^{-1}) via one-block construction ──
    // eps = L_Q z1 + C^T z2 ~ N(0, P);  u_offset = P^{-1} eps ~ N(0, P^{-1})
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> d(0.0, 1.0);
    
    MatrixXd Z1(nv,        niter_);
    MatrixXd Z2(n_regions, niter_);
    for(int i = 0; i < Z1.size(); ++i) Z1.data()[i] = d(gen);
    for(int i = 0; i < Z2.size(); ++i) Z2.data()[i] = d(gen);
    
    // v1 = P^T_perm * L_Q * z1   ~ N(0, Q),  using Q = P^T L L^T P
    MatrixXd LZ1 = chol_Q.matrixL() * Z1;
    MatrixXd v1  = chol_Q.permutationPinv() * LZ1;
    MatrixXd v2  = SparseMatrix<double>(C_sp.transpose()) * Z2;
    MatrixXd eps = v1 + v2;
    
    // u_offset = P^{-1} eps via Woodbury (reuse Y_dense, llt_K)
    MatrixXd Qinv_eps = chol_Q.solve(eps);
    MatrixXd C_Qe     = C_sp * Qinv_eps;
    MatrixXd u_offset = Qinv_eps - Y_dense * llt_K.solve(C_Qe);
    
    // Proposal log-density (up to const):  -0.5 (u-u_mean)^T P (u-u_mean) = -0.5 eps^T u_offset
    u_loglik_ = -0.5 * (eps.array() * u_offset.array()).colwise().sum().matrix();
    
    u_       = u_offset.colwise() + u_mean_;
    scaled_u_ = ZA * u_;
    
    // ── Importance weights ──
    // log f(u|θ) = 0.5 log det Q - 0.5 u^T Q u  (const drops in normalisation)
    u_weight_.setZero();
    ArrayXd ll = log_likelihood();
    const SparseMatrix<double>& Q_mat = covariance.Q_mat;
    for(int i = 0; i < u_.cols(); i++){
      double uQu     = u_.col(i).dot(Q_mat * u_.col(i));
      double llprior = -0.5 * uQu;
      u_weight_(i)   = ll(i) + llprior - u_loglik_(i);
    }
    u_weight_ -= u_weight_.maxCoeff();
    u_weight_  = u_weight_.exp();
    u_weight_ *= 1.0 / u_weight_.sum();
    
    if(u_mean_.hasNaN()){
      Rcpp::Rcout << "DIAGNOSTIC: NaN in u_mean_" << std::endl;
      Rcpp::stop("usample (SPDE): NaN detected in u_mean_");
    }
    double ess = 1.0 / (u_weight_ * u_weight_).sum();
    if(ess < 2.0){
      Rcpp::Rcout << "WARNING: Very low ESS = " << ess << " / " << u_weight_.size() << std::endl;
    }
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
  log_det_P_ = 2.0 * llt_CCt.matrixL().toDenseMatrix().diagonal().array().log().sum();
  
  // ── Posterior variance diagonals ──
  // Precision P = C^T C + I
  // Woodbury: P^{-1} = I - C^T (C C^T + I)^{-1} C
  // diag(P^{-1})_i = 1 - c_i^T M^{-1} c_i,   c_i = i-th column of C
  
  MatrixXd MinvC = llt_CCt.solve(C);                                    // (n_regions[*T]) × n_cols
  ArrayXd  quad  = (C.array() * MinvC.array()).colwise().sum();         // length n_cols
  u_var_diag_    = (1.0 - quad).matrix();
  
  // Var(Zu) = ZL * P^{-1} * ZL^T
  // diag = rowwise ||ZL||^2  -  rowwise (K .* (M^{-1} K^T)^T).sum(),   K = ZL * C^T
  ArrayXd  diag_full = ZL.array().square().rowwise().sum();             // n_A[*T]
  MatrixXd K         = ZL * C.transpose();                              // n_A[*T] × n_regions[*T]
  MatrixXd MinvKt    = llt_CCt.solve(K.transpose());                    // n_regions[*T] × n_A[*T]
  ArrayXd  diag_corr = (K.array() * MinvKt.transpose().array()).rowwise().sum();
  zu_var_            = (diag_full - diag_corr).matrix();
  
  // ── Resize storage ──
  if(u_.cols() != niter_){
    u_.resize(NoChange, niter_);
    u_solve_.resize(NoChange, niter_);
    scaled_u_.resize(NoChange, niter_);
    u_weight_.resize(niter_);
    u_loglik_.resize(niter_);
  }
  
  if(niter_ == 1){
    // Laplace: posterior mode only
    u_.col(0)    = u_mean_;
    u_loglik_(0) = 0.0;
    u_weight_(0) = 1.0;
    
    scaled_u_ = covariance.Lu(u_);
#ifdef GLMMR13
    u_solve_ = covariance.solve(scaled_u_);
#else
    MatrixXd D = covariance.D();
    LLT<MatrixXd> llt_D;
    llt_D.compute(D);
    u_solve_ = llt_D.solve(scaled_u_);
#endif
  } else {
    MatrixXd unew(u_.rows(), niter_);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> d(0.0, 1.0);
    double* data = unew.data();
    for(int i = 0; i < unew.size(); ++i) data[i] = d(gen);
    
    int m = C.rows();  // n_regions * T
    int n_samples = unew.cols();
    
    // Generate w ~ N(0, I)
    MatrixXd w(m, n_samples);
    data = w.data();
    for(int i = 0; i < w.size(); ++i) data[i] = d(gen);
    
    // RHS = z + C'w
    MatrixXd RHS = unew + C.transpose() * w;
    
    // Apply Woodbury: (C'C + I)^{-1} RHS = RHS - C'(CC'+I)^{-1}(C * RHS)
    MatrixXd CRHS = C * RHS;
    MatrixXd tmp  = llt_CCt.solve(CRHS);
    unew = RHS - C.transpose() * tmp;
    
    u_.noalias() = unew;
    
    MatrixXd Cu = C * u_;
    ArrayXd Cu_norms = Cu.colwise().squaredNorm();
    ArrayXd u_norms  = u_.colwise().squaredNorm();
    u_loglik_ = -0.5 * (Cu_norms + u_norms);
    
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
      double llmod   = ll(i);
      double llprior = -0.5 * scaled_u_.col(i).dot(u_solve_.col(i));
      u_weight_(i)   = llmod + llprior - u_loglik_(i);
    }
    u_weight_ -= u_weight_.maxCoeff();
    u_weight_  = u_weight_.exp();
    u_weight_ *= 1.0 / u_weight_.sum();
  }
  
  // ── Diagnostics (apply to both paths) ──
  if(u_mean_.hasNaN()) {
    Rcpp::Rcout << "=== DIAGNOSTIC: NaN in u_mean_ ===" << std::endl;
    // ... (rest of your diagnostic block unchanged) ...
    Rcpp::stop("usample: NaN detected in u_mean_");
  }
  
  if(u_.hasNaN()) {
    Rcpp::Rcout << "=== DIAGNOSTIC: NaN in u_ ===" << std::endl;
    Rcpp::stop("usample: NaN detected in u_");
  }
  
  if(u_weight_.hasNaN() || !u_weight_.isFinite().all()) {
    Rcpp::Rcout << "=== DIAGNOSTIC: NaN/Inf in u_weight_ ===" << std::endl;
    Rcpp::stop("usample: NaN/Inf detected in u_weight_");
  }
  
  double ess = 1.0 / (u_weight_ * u_weight_).sum();
  if(ess < 2.0 && niter_ > 1) {
    Rcpp::Rcout << "=== WARNING: Very low ESS === " << ess
                << " / " << u_weight_.size() << std::endl;
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
  
  for(int i = 0; i < niter; i++){
    double zu_mean = scaled_u_.col(i).mean();
    zd.col(i).array() -= zu_mean;
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

// ── nr_theta specialisation for hsgpCovariance ─────────────────────
template<>
inline void rts::regionModel<glmmr::hsgpCovariance>::nr_theta()
{
  const int n_regions = weights.rows();
  const int n_cells   = weights.cols();
  const int Mspec     = covariance.Q();
  const int npars     = covariance.npar();               // log(sigma^2), log(ell)
  const int n_iter    = u_.cols();
  const int total = scaled_u_.rows();  // n_cells * T
  const int T     = total / n_cells;
  
  MatrixXd Phi     = covariance.ZPhi();  // n_cells × M  (includes Z if present)
  ArrayXd  Lambda  = covariance.LambdaSPD();
  ArrayXd  inv_lam = 1.0 / Lambda;
  
  // SPD derivatives
  ArrayXXd dLambda(Mspec, npars);
  ArrayXXd dsqrtLam(Mspec, npars);
  ArrayXd  sqrtLam = Lambda.sqrt();
  for(int k = 0; k < Mspec; k++){
    dblvec deriv = covariance.d_spd_nD(k);
    for(int j = 0; j < npars; j++){
      dLambda(k, j)  = deriv[j];
      dsqrtLam(k, j) = deriv[j] / (2.0 * sqrtLam(k));
    }
  }
  
  // Precompute Phi * diag(dsqrt(Lambda)/dtheta_j) for obs-space Hessian
  std::vector<MatrixXd> Phi_dj(npars);
  for(int j = 0; j < npars; j++){
    Phi_dj[j] = Phi * dsqrtLam.col(j).matrix().asDiagonal();
  }
  
  // ── Prior gradient (deterministic + MC) ──
  ArrayXd inv_lam2 = inv_lam.square();
  VectorXd grad = VectorXd::Zero(npars);
  for(int j = 0; j < npars; j++){
    grad(j) = -0.5 * (dLambda.col(j) * inv_lam).sum();
  }
  
  if(n_iter == 1){
    ArrayXd second_moment = u_.col(0).array().square() + u_var_diag_.array();
    for(int j = 0; j < npars; j++){
      grad(j) += 0.5 * (second_moment * dLambda.col(j) * inv_lam2).sum();
    }
  } else {
    for(int i = 0; i < n_iter; i++){
      double w_i = u_weight_(i);
      ArrayXd u2 = u_.col(i).array().square();
      for(int j = 0; j < npars; j++){
        grad(j) += 0.5 * w_i * (u2 * dLambda.col(j) * inv_lam2).sum();
      }
    }
  }
  
  // ── Hessian: prior Fisher + observation-space Fisher ──
  MatrixXd Hess = MatrixXd::Zero(npars, npars);
  for(int j = 0; j < npars; j++){
    for(int l = j; l < npars; l++){
      double h = 0.5 * (dLambda.col(j) * dLambda.col(l) * inv_lam2).sum();
      Hess(j, l) = h;
      if(j != l) Hess(l, j) = h;
    }
  }
  
  // Precompute Phi * diag((dsqrt(Lambda)/dtheta_j) / sqrt(Lambda))
  // so that Phi_dj_u[j] * u gives Phi * (dsqrtLam/dtheta_j) * v, where v = u/sqrtLam
  std::vector<MatrixXd> Phi_dj_u(npars);
  for(int j = 0; j < npars; j++){
    ArrayXd scale = dsqrtLam.col(j) / sqrtLam;
    Phi_dj_u[j] = Phi * scale.matrix().asDiagonal();
  }
  
  // Observation-space Fisher
  ArrayXd xb = (X * beta + offset).array();
  for(int i = 0; i < n_iter; i++){
    ArrayXd eta_s    = xb + scaled_u_.col(i).array();
    ArrayXd lambda_s = eta_s.exp();
    VectorXd lambda_r = weights * lambda_s.matrix();
    ArrayXd  inv_lr   = lambda_r.array().inverse();
    double w_i = u_weight_(i);
    
    std::vector<VectorXd> dlr(npars);
    for(int j = 0; j < npars; j++){
      VectorXd deta_cell = Phi_dj_u[j] * u_.col(i);   // <-- was Phi_dj[j]
      ArrayXd scaled_arr = lambda_s * deta_cell.array();
      VectorXd dlr_j(n_regions * T);
      for(int t = 0; t < T; t++){
        VectorXd sc_t = scaled_arr.segment(t * n_cells, n_cells).matrix();
        dlr_j.segment(t * n_regions, n_regions) = weights * sc_t;
      }
      dlr[j] = dlr_j;
    }
    
    for(int j = 0; j < npars; j++){
      for(int l = j; l < npars; l++){
        double h_obs = (dlr[j].array() * dlr[l].array() * inv_lr).sum();
        Hess(j, l) += w_i * h_obs;
        if(j != l) Hess(l, j) += w_i * h_obs;
      }
    }
  }
  
  
  
  // ── Newton-Raphson step with damping ──
  covariance.infomat_theta = Hess;
  bool logpars = covariance.all_log_re();
  VectorXd logtheta(npars);
  if(logpars){
    logtheta = Map<VectorXd>(covariance.parameters_.data(), npars);
  } else {
    for(int j = 0; j < npars; j++) logtheta(j) = log(covariance.parameters_[j]);
  }
  
  double lambda_damp = 1e-4 * Hess.diagonal().array().abs().maxCoeff();
  MatrixXd Hess_reg  = Hess;
  Hess_reg.diagonal().array() += lambda_damp;
  
  VectorXd step;
  Eigen::LLT<MatrixXd> llt_H(Hess_reg);
  if(llt_H.info() == Eigen::Success){
    step = llt_H.solve(grad);
  } else {
    step = grad.array() / Hess.diagonal().array().abs().max(1e-6);
  }
  
  double max_step = 1.0;
  if(n_iter > 1 && step.array().abs().maxCoeff() > max_step)
    step *= max_step / step.array().abs().maxCoeff();
  
  logtheta += step;
  
  if(logpars){
    dblvec newpars(logtheta.data(), logtheta.data() + logtheta.size());
    covariance.update_parameters(newpars);
  } else {
    covariance.update_parameters(logtheta.array().exp());
  }
  
  gradients.tail(npars) = grad.array();
  
  ArrayXd new_lambda = covariance.LambdaSPD();
  ArrayXd new_inv_lam = 1.0 / new_lambda;
  double ll_det = -0.5 * new_lambda.log().sum();
  ll_theta = 0;
  for(int i = 0; i < n_iter; i++){
    double qf = (u_.col(i).array().square() * new_inv_lam).sum();
    ll_theta += u_weight_(i) * (ll_det - 0.5 * qf);
  }
}

template<>
inline void rts::regionModel<glmmr::spdeCovariance>::nr_theta()
{
  const int npars  = covariance.npar();   // {log sigma^2, log lambda}
  const int n_iter = u_.cols();
  const int nv     = covariance.Q();
  
  // Refresh prior structure at current parameters
  if(!covariance.chol_Q_current) covariance.refactor_Q();
  
  // Common random numbers across gradient and Hessian traces this iteration
  covariance.refresh_probes();
  
  // ── Prior gradient ──
  // ∂log f(u|θ)/∂θ_j = 0.5 tr(Q^{-1} ∂Q/∂θ_j) - 0.5 u^T (∂Q/∂θ_j) u
  //   ∂Q/∂log σ²  = -Q   ⇒ tr term = -nv/2,   u^T ∂Q u = -u^T Q u
  //   ∂Q/∂log λ : one Hutchinson trace + analytic quadratic form
  double tr_dQ_lambda = covariance.trace_Qinv_dQ_log_lambda_hutch(/*K=*/30);
  
  VectorXd grad = VectorXd::Zero(npars);
  // Deterministic trace contribution
  grad(0) =  0.5 * (-static_cast<double>(nv));
  grad(1) =  0.5 * tr_dQ_lambda;
  
  // MC quadratic form contribution
  if(n_iter == 1){
    // Laplace: u^T A u  →  u_mean^T A u_mean + tr(A * Var(u))
    // For SPDE we approximate using the posterior mode only; tr term contributes
    // to consistency but is dropped here matching the HSGP "second moment" treatment.
    // (If you want the trace correction, add it via Hutchinson on Var(u).)
    VectorXd u0 = u_.col(0);
    grad(0) += 0.5 * covariance.quad_form_Q(u0);            // -(-u^T Q u)/2 = +u^T Q u / 2
    grad(1) -= 0.5 * covariance.quad_form_dQ_log_lambda(u0);
  } else {
    for(int i = 0; i < n_iter; i++){
      VectorXd ui = u_.col(i);
      double w_i  = u_weight_(i);
      grad(0) += 0.5 * w_i * covariance.quad_form_Q(ui);
      grad(1) -= 0.5 * w_i * covariance.quad_form_dQ_log_lambda(ui);
    }
  }
  
  // ── Expected Hessian (prior Fisher information) ──
  // I_jl = 0.5 tr(Q^{-1} ∂Q/∂θ_j Q^{-1} ∂Q/∂θ_l)
  // Analytic: I(σ²,σ²) = nv/2,  I(σ²,λ) = -0.5 tr(Q^{-1} ∂Q/∂log λ)
  // Hutchinson: I(λ,λ) = 0.5 tr(Q^{-1} Qλ Q^{-1} Qλ)
  double I_lambda_lambda = 0.5 * covariance.trace_Qinv_dQ_Qinv_dQ_lambda(/*K=*/30);
  
  MatrixXd Hess(npars, npars);
  Hess(0, 0) = 0.5 * static_cast<double>(nv);
  Hess(0, 1) = -0.5 * tr_dQ_lambda;
  Hess(1, 0) = Hess(0, 1);
  Hess(1, 1) = I_lambda_lambda;
  
  // ── Newton-Raphson step with Levenberg damping and clamp ──
  covariance.infomat_theta = Hess;
  bool logpars = covariance.all_log_re();
  VectorXd logtheta(npars);
  if(logpars){
    logtheta = Map<VectorXd>(covariance.parameters_.data(), npars);
  } else {
    for(int j = 0; j < npars; j++) logtheta(j) = log(covariance.parameters_[j]);
  }
  
  double lambda_damp = 1e-4 * Hess.diagonal().array().abs().maxCoeff();
  MatrixXd Hess_reg  = Hess;
  Hess_reg.diagonal().array() += lambda_damp;
  
  VectorXd step;
  Eigen::LLT<MatrixXd> llt_H(Hess_reg);
  if(llt_H.info() == Eigen::Success){
    step = llt_H.solve(grad);
  } else {
    step = grad.array() / Hess.diagonal().array().abs().max(1e-6);
  }
  
  double max_step = 1.0;
  if(n_iter > 1 && step.array().abs().maxCoeff() > max_step)
    step *= max_step / step.array().abs().maxCoeff();
  
  logtheta += step;
  
  if(logpars){
    dblvec newpars(logtheta.data(), logtheta.data() + logtheta.size());
    covariance.update_parameters(newpars);
  } else {
    covariance.update_parameters(logtheta.array().exp());
  }
  
  gradients.tail(npars) = grad.array();
  
  // ── ll_theta at new parameters ──
  // log f(u|θ) = 0.5 log det Q - 0.5 u^T Q u + const
  // log_determinant() returns log det D = -log det Q (matching base class convention)
  double log_det_Q = -covariance.log_determinant();
  double ll_det    = 0.5 * log_det_Q;
  
  ll_theta = 0.0;
  for(int i = 0; i < n_iter; i++){
    double qf = covariance.quad_form_Q(u_.col(i));
    ll_theta += u_weight_(i) * (ll_det - 0.5 * qf);
  }
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
  double ll_old = -std::numeric_limits<double>::infinity();
  // start at ml values
  if(trace > 0)Rcpp::Rcout << "\nStarting values\nBeta: " << beta.transpose() << "\nTheta: ";
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
    if constexpr (std::is_same_v<cov, glmmr::ar1Covariance>){
      if(trace>0) Rcpp::Rcout << "\nRho: " << covariance.rho;
    }
    
    if(niter > 1){
      if(trace > 0)Rcpp::Rcout << "\nu: " << u_.topRows(5).rowwise().mean().transpose();
      if(trace > 0)Rcpp::Rcout << "\nLog-likelihood: " << ll_beta << " | " << ll_theta << std::endl;
      if (iter > 2) {
        converged = check_convergence(tol,hist,iter,k0);
        if (converged && trace > 0) Rcpp::Rcout << "\nCONVERGED!";
      }
    } else {
      double grad_norm = gradients.abs().maxCoeff();
      if(trace > 0) Rcpp::Rcout << "\nmax |grad|: " << grad_norm;
      converged = grad_norm < tol;
      double ll_new = ll_beta + ll_theta - 0.5 * log_det_P_;
      //converged = ll_new - ll_old < tol;
      if(trace > 0)Rcpp::Rcout << "\nLog likelihood (new | previous): " << ll_new << " | " << ll_old << " Diff: " << ll_new - ll_old << std::endl;
      ll_old = ll_new;
      if(trace > 0) if (converged) Rcpp::Rcout << "\nCONVERGED!";
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
SEXP regionModel_hsgp__new(const std::string& formula_,
                           const Eigen::ArrayXXd& data_,
                           const std::vector<std::string>& colnames_,
                           const Eigen::MatrixXd& X_,
                           const Eigen::ArrayXd& y_,
                           const int niter_,
                           const std::vector<int>& m_,
                           const double L_boundary_)
{
  using namespace rts;
  Rcpp::XPtr<regionModel<glmmr::hsgpCovariance>> ptr(
      new regionModel<glmmr::hsgpCovariance>(
          formula_, data_, colnames_, X_, y_, niter_, m_, L_boundary_), 
          true
  );
  if(data_.cols() == 3){
    ptr->covariance.set_anisotropic(true);
    ptr->gradients = ArrayXd::Zero(ptr->X.cols() + ptr->covariance.npar());
    Rcpp::Rcout << "\nUsing anisotropic HSGP model";
  }
  return ptr;
}

// [[Rcpp::export]]
SEXP regionModel_spde__new(const std::string& formula_,
                           const Eigen::ArrayXXd& data_,
                           const std::vector<std::string>& colnames_,
                           const Eigen::MatrixXd& X_,
                           const Eigen::ArrayXd& y_,
                           const int niter_,
                           const Eigen::SparseMatrix<double>& A_loc_,
                           const Eigen::VectorXd& C_diag_,
                           const Eigen::SparseMatrix<double>& G_,
                           const int alpha_ = 2)
{
  using namespace rts;
  Rcpp::XPtr<regionModel<glmmr::spdeCovariance>> ptr(
      new regionModel<glmmr::spdeCovariance>(
          formula_, data_, colnames_, X_, y_, niter_, A_loc_, C_diag_, G_),
          true
  );
  // After spde_data, Q is known; resize gradients accordingly
  ptr->gradients = Eigen::ArrayXd::Zero(ptr->X.cols() + ptr->covariance.npar());
  return ptr;
}

template<>
inline Eigen::VectorXd rts::regionModel<glmmr::ar1Covariance>::zu_variance_full()
{
  using Eigen::MatrixXd; using Eigen::VectorXd; using Eigen::ArrayXd;
  
  const int n_A       = weights.cols();         // cells per period
  const int n_regions = weights.rows();
  const int n_cells   = X.rows();               // = n_A * T
  const int T         = n_cells / n_A;
  const int Pdim      = X.cols();
  const int K         = u_.cols();
  
  MatrixXd ZL = covariance.ZL();                // n_cells × n_u (stacked by time)
  const int nu = ZL.cols();
  
  // ── Posterior-mean cell-level intensities (stacked) ──
  double wsum = 0.0; for(int k = 0; k < K; ++k) wsum += u_weight_(k);
  ArrayXd xb = (X * beta + offset).array();
  ArrayXd lambda_s_bar = ArrayXd::Zero(n_cells);
  for(int k = 0; k < K; ++k){
    ArrayXd lk = (xb + scaled_u_.col(k).array() - scaled_u_.col(k).mean()).exp();
    lambda_s_bar += (u_weight_(k) / wsum) * lk;
  }
  
  // λ_r per period (length n_regions * T)
  VectorXd lambda_r_bar(n_regions * T);
  for(int t = 0; t < T; ++t){
    VectorXd lst = lambda_s_bar.segment(t * n_A, n_A).matrix();
    lambda_r_bar.segment(t * n_regions, n_regions) = weights * lst;
  }
  ArrayXd inv_lr = lambda_r_bar.array().inverse();
  
  // ── tilde_X (n_regions*T × P) and tilde_ZL (n_regions*T × n_u) ──
  // Both are weights·diag(λ_s)·M applied per period, then row-scaled by 1/λ_r
  MatrixXd lsX  = (X.array().colwise()  * lambda_s_bar).matrix();   // n_cells × P
  MatrixXd lsZL = (ZL.array().colwise() * lambda_s_bar).matrix();   // n_cells × n_u
  
  MatrixXd tilde_X (n_regions * T, Pdim);
  MatrixXd tilde_ZL(n_regions * T, nu);
  for(int t = 0; t < T; ++t){
    tilde_X .middleRows(t * n_regions, n_regions) = weights * lsX .middleRows(t * n_A, n_A);
    tilde_ZL.middleRows(t * n_regions, n_regions) = weights * lsZL.middleRows(t * n_A, n_A);
  }
  tilde_X .array().colwise() *= inv_lr;
  tilde_ZL.array().colwise() *= inv_lr;
  
  // ── Schur ──
  VectorXd Wr       = lambda_r_bar;
  MatrixXd WrtildeX = (tilde_X.array().colwise() * Wr.array()).matrix();
  MatrixXd XtWX     = tilde_X.transpose() * WrtildeX;
  
  MatrixXd Wr_tildeZL = (tilde_ZL.array().colwise() * Wr.array()).matrix();
  MatrixXd A_mat = tilde_ZL.transpose() * Wr_tildeZL;
  A_mat.diagonal().array() += 1.0;
  Eigen::LLT<MatrixXd> llt_A(A_mat);
  
  MatrixXd B     = tilde_ZL.transpose() * WrtildeX;
  MatrixXd AinvB = llt_A.solve(B);
  MatrixXd Schur = XtWX - B.transpose() * AinvB;
  Eigen::LLT<MatrixXd> llt_Schur(Schur);
  
  // ── Var(η) at cells (stacked) ──
  // For AR1 the "evaluation projector" at cell i is row i of ZL, exactly as in dense.
  // ZL is n_cells × n_u; the stacking is implicit in ZL's rows.
  MatrixXd Ainv_ZLT = llt_A.solve(ZL.transpose());                  // n_u × n_cells
  VectorXd term1(n_cells);
  for(int i = 0; i < n_cells; ++i){
    term1(i) = ZL.row(i).dot(Ainv_ZLT.col(i));
  }
  
  MatrixXd Bt_Ainv_ZLT = B.transpose() * Ainv_ZLT;                  // P × n_cells
  MatrixXd diff_t      = X.transpose() - Bt_Ainv_ZLT;               // P × n_cells
  MatrixXd Sdiff       = llt_Schur.solve(diff_t);
  VectorXd term_beta   = (diff_t.array() * Sdiff.array()).colwise().sum();
  
  return term1 + term_beta;
}

template<>
inline VectorXd rts::regionModel<glmmr::spdeCovariance>::zu_variance_full(
    const Eigen::SparseMatrix<double>& A_pred,
    const Eigen::MatrixXd& X_pred)
{
  using Eigen::MatrixXd; using Eigen::VectorXd; using Eigen::ArrayXd;
  using Eigen::SparseMatrix;
  
  const int n_pred = A_pred.rows();
  const int K      = u_.cols();
  const SparseMatrix<double>& ZA = covariance.ZA_;   // mesh-side, used only for λ̄_s
  
  // ── Posterior-mean cell-level intensities (mesh side, Jensen-centered) ──
  double wsum = 0.0; for(int k = 0; k < K; ++k) wsum += u_weight_(k);
  ArrayXd xb = (X * beta + offset).array();
  ArrayXd lambda_s_bar = ArrayXd::Zero(X.rows());
  for(int k = 0; k < K; ++k){
    ArrayXd lk = (xb + scaled_u_.col(k).array() - scaled_u_.col(k).mean()).exp();
    lambda_s_bar += (u_weight_(k) / wsum) * lk;
  }
  VectorXd lambda_r_bar = weights * lambda_s_bar.matrix();
  ArrayXd  inv_lr       = lambda_r_bar.array().inverse();
  
  // ── Region-level working design + posterior precision (unchanged) ──
  MatrixXd lsX     = (X.array().colwise() * lambda_s_bar).matrix();
  MatrixXd tilde_X = ((weights * lsX).array().colwise() * inv_lr).matrix();
  
  SparseMatrix<double> Lam_s_ZA = lambda_s_bar.matrix().asDiagonal() * ZA;
  SparseMatrix<double> tilde_A  = inv_lr.matrix().asDiagonal() * (weights * Lam_s_ZA);
  
  VectorXd Wr = lambda_r_bar;
  MatrixXd WrtildeX = (tilde_X.array().colwise() * Wr.array()).matrix();
  MatrixXd XtWX     = tilde_X.transpose() * WrtildeX;
  
  // Form P_region = tilde_A^T diag(Wr) tilde_A + Q
  SparseMatrix<double> WrtA_sp = Wr.asDiagonal() * tilde_A;   // n_R × nv, sparse
  SparseMatrix<double> P_region = SparseMatrix<double>(tilde_A.transpose() * WrtA_sp)
    + covariance.Q_mat;
  
  // Factorise (pattern matches Q's pattern + tilde_A's column-wise outer product;
  // not identical to refactor_P's pattern, so use a fresh solver)
  Eigen::SimplicialLLT<SparseMatrix<double>, Eigen::Lower, Eigen::AMDOrdering<int>> chol_P_region;
  chol_P_region.compute(P_region);
  if (chol_P_region.info() != Eigen::Success) {
    Rcpp::stop("region-level P factorisation failed: Wr range [%.3e, %.3e], "
                 "Q diag range [%.3e, %.3e]",
                 Wr.minCoeff(), Wr.maxCoeff(),
                 covariance.Q_mat.diagonal().minCoeff(),
                 covariance.Q_mat.diagonal().maxCoeff());
  }
  
  // Replace your current uses of chol_P with chol_P_region
  MatrixXd B     = MatrixXd(tilde_A.transpose() * WrtildeX);     // nv × P
  MatrixXd PinvB = chol_P_region.solve(B);
  auto& chol_P = covariance.chol_P;
  
  
  MatrixXd Schur = XtWX - B.transpose() * PinvB;                 // P × P  ( = V_β⁻¹ )
  Eigen::LLT<MatrixXd> llt_Schur(Schur);
  
  // ── Evaluate at prediction points via A_pred, X_pred ──
  // term1_i = (A_pred P⁻¹ A_predᵀ)_ii
  MatrixXd Pinv_ApT = chol_P.solve(MatrixXd(A_pred.transpose())); // nv × n_pred
  VectorXd term1 = VectorXd::Zero(n_pred);
  for(int j = 0; j < A_pred.outerSize(); ++j){
    for(SparseMatrix<double>::InnerIterator it(A_pred, j); it; ++it){
      term1(it.row()) += it.value() * Pinv_ApT(j, it.row());
    }
  }
  
  // c_pred (n_pred × P) = A_pred * PinvB ;  diff = X_pred − c_pred
  MatrixXd c_pred = A_pred * PinvB;
  MatrixXd diff_t = (X_pred - c_pred).transpose();                 // P × n_pred
  MatrixXd Sdiff  = llt_Schur.solve(diff_t);                       // P × n_pred
  VectorXd term_beta = (diff_t.array() * Sdiff.array()).colwise().sum().transpose();
  
  return term1 + term_beta;   // length n_pred = Var(η) at grid points
}

template<>
inline Eigen::VectorXd rts::regionModel<glmmr::hsgpCovariance>::zu_variance_full()
{
  using Eigen::MatrixXd; using Eigen::VectorXd; using Eigen::ArrayXd;
  
  const int n_cells = X.rows();
  const int K       = u_.cols();
  
  MatrixXd Phi = covariance.ZPhi();              // n_cells × M
  ArrayXd  Lambda = covariance.LambdaSPD();      // length M
  const int Mdim = Phi.cols();
  
  // ── Posterior-mean cell-level intensities (Jensen-centered per sample) ──
  double wsum = 0.0; for(int k = 0; k < K; ++k) wsum += u_weight_(k);
  ArrayXd xb = (X * beta + offset).array();
  ArrayXd lambda_s_bar = ArrayXd::Zero(n_cells);
  for(int k = 0; k < K; ++k){
    ArrayXd lk = (xb + scaled_u_.col(k).array() - scaled_u_.col(k).mean()).exp();
    lambda_s_bar += (u_weight_(k) / wsum) * lk;
  }
  VectorXd lambda_r_bar = weights * lambda_s_bar.matrix();
  ArrayXd  inv_lr       = lambda_r_bar.array().inverse();
  
  // ── Region-level working design (same construction as SPDE version) ──
  MatrixXd lsX     = (X.array().colwise() * lambda_s_bar).matrix();
  MatrixXd tilde_X = ((weights * lsX).array().colwise() * inv_lr).matrix();          // n_regions × P
  
  MatrixXd Lam_s_Phi = (lambda_s_bar.matrix().asDiagonal() * Phi);                   // n_cells × M
  MatrixXd tilde_Phi = (inv_lr.matrix().asDiagonal() * (weights * Lam_s_Phi));       // n_regions × M
  
  VectorXd Wr       = lambda_r_bar;
  MatrixXd WrtildeX = (tilde_X.array().colwise() * Wr.array()).matrix();
  MatrixXd XtWX     = tilde_X.transpose() * WrtildeX;                                // P × P
  
  // A = tilde_Phi' diag(Wr) tilde_Phi + diag(1/Lambda)
  MatrixXd Wr_tildePhi = (tilde_Phi.array().colwise() * Wr.array()).matrix();
  MatrixXd A_mat = tilde_Phi.transpose() * Wr_tildePhi;
  A_mat.diagonal().array() += 1.0 / Lambda;
  Eigen::LLT<MatrixXd> llt_A(A_mat);
  
  MatrixXd B     = tilde_Phi.transpose() * WrtildeX;                                  // M × P
  MatrixXd AinvB = llt_A.solve(B);                                                    // M × P
  MatrixXd Schur = XtWX - B.transpose() * AinvB;                                      // P × P  (= V_β⁻¹)
  Eigen::LLT<MatrixXd> llt_Schur(Schur);
  
  // ── Evaluate Var(η) at cells: A_pred = Phi (dense), X_pred = X ──
  // term1_i = (Phi A⁻¹ Phiᵀ)_ii
  MatrixXd Ainv_PhiT = llt_A.solve(Phi.transpose());                                  // M × n_cells
  VectorXd term1(n_cells);
  for(int i = 0; i < n_cells; ++i){
    term1(i) = Phi.row(i).dot(Ainv_PhiT.col(i));
  }
  
  // c_i = (B' A⁻¹ Phiᵀ)_:,i  ;  diff = X - c
  MatrixXd Bt_Ainv_PhiT = B.transpose() * Ainv_PhiT;                                  // P × n_cells
  MatrixXd diff_t       = X.transpose() - Bt_Ainv_PhiT;                                // P × n_cells
  MatrixXd Sdiff        = llt_Schur.solve(diff_t);                                     // P × n_cells
  VectorXd term_beta    = (diff_t.array() * Sdiff.array()).colwise().sum();
  
  return term1 + term_beta;
}

template<>
inline Eigen::VectorXd rts::regionModel<glmmr::Covariance>::zu_variance_full()
{
  using Eigen::MatrixXd; using Eigen::VectorXd; using Eigen::ArrayXd;
  
  const int n_cells = X.rows();
  const int K       = u_.cols();
  
  MatrixXd ZL = covariance.ZL();                  // n_cells × n_u
  const int Mdim = ZL.cols();
  
  // ── Posterior-mean cell-level intensities ──
  double wsum = 0.0; for(int k = 0; k < K; ++k) wsum += u_weight_(k);
  ArrayXd xb = (X * beta + offset).array();
  ArrayXd lambda_s_bar = ArrayXd::Zero(n_cells);
  for(int k = 0; k < K; ++k){
    ArrayXd lk = (xb + scaled_u_.col(k).array() - scaled_u_.col(k).mean()).exp();
    lambda_s_bar += (u_weight_(k) / wsum) * lk;
  }
  VectorXd lambda_r_bar = weights * lambda_s_bar.matrix();
  ArrayXd  inv_lr       = lambda_r_bar.array().inverse();
  
  // ── Region-level working design ──
  MatrixXd lsX     = (X.array().colwise() * lambda_s_bar).matrix();
  MatrixXd tilde_X = ((weights * lsX).array().colwise() * inv_lr).matrix();
  
  MatrixXd Lam_s_ZL = (lambda_s_bar.matrix().asDiagonal() * ZL);
  MatrixXd tilde_ZL = (inv_lr.matrix().asDiagonal() * (weights * Lam_s_ZL));
  
  VectorXd Wr       = lambda_r_bar;
  MatrixXd WrtildeX = (tilde_X.array().colwise() * Wr.array()).matrix();
  MatrixXd XtWX     = tilde_X.transpose() * WrtildeX;
  
  // A = tilde_ZL' diag(Wr) tilde_ZL + I
  MatrixXd Wr_tildeZL = (tilde_ZL.array().colwise() * Wr.array()).matrix();
  MatrixXd A_mat = tilde_ZL.transpose() * Wr_tildeZL;
  A_mat.diagonal().array() += 1.0;
  Eigen::LLT<MatrixXd> llt_A(A_mat);
  
  MatrixXd B     = tilde_ZL.transpose() * WrtildeX;
  MatrixXd AinvB = llt_A.solve(B);
  MatrixXd Schur = XtWX - B.transpose() * AinvB;
  Eigen::LLT<MatrixXd> llt_Schur(Schur);
  
  // ── Var(η) at cells: A_pred = ZL, X_pred = X ──
  MatrixXd Ainv_ZLT = llt_A.solve(ZL.transpose());
  VectorXd term1(n_cells);
  for(int i = 0; i < n_cells; ++i){
    term1(i) = ZL.row(i).dot(Ainv_ZLT.col(i));
  }
  
  MatrixXd Bt_Ainv_ZLT = B.transpose() * Ainv_ZLT;
  MatrixXd diff_t      = X.transpose() - Bt_Ainv_ZLT;
  MatrixXd Sdiff       = llt_Schur.solve(diff_t);
  VectorXd term_beta   = (diff_t.array() * Sdiff.array()).colwise().sum();
  
  return term1 + term_beta;
}

using namespace Rcpp;

// [[Rcpp::export]]
SEXP regionModel__zu_variance_full(SEXP xp, int type)
{
  if(type == 2){          // HSGP
    XPtr<rts::regionModel<glmmr::hsgpCovariance>> ptr(xp);
    return wrap(ptr->zu_variance_full());
  } else {                // dense (type 0 or 1)
    XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    return wrap(ptr->zu_variance_full());
  }
}

// [[Rcpp::export]]
SEXP regionModel_spde__zu_variance_full(SEXP xp,
                                        SEXP A_pred_,
                                        SEXP X_pred_)
{
  XPtr<rts::regionModel<glmmr::spdeCovariance>> ptr(xp);
  Eigen::SparseMatrix<double> A_pred = as<Eigen::SparseMatrix<double>>(A_pred_);
  Eigen::MatrixXd            X_pred = as<Eigen::MatrixXd>(X_pred_);
  return wrap(ptr->zu_variance_full(A_pred, X_pred));
}

// [[Rcpp::export]]
void regionModel__set_weights(SEXP xp, 
                              Rcpp::IntegerVector i, 
                              Rcpp::IntegerVector p, 
                              Rcpp::NumericVector x, 
                              int nrow, 
                              int ncol, 
                              int type = 0) {
  
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
  } else if (type == 1) {
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    ptr->weights = W;
    ptr->init_beta();
  } else if (type == 2) {
    Rcpp::XPtr<rts::regionModel<glmmr::hsgpCovariance>> ptr(xp);
    ptr->weights = W;
    ptr->init_beta();
  } else {
    Rcpp::XPtr<rts::regionModel<glmmr::spdeCovariance>> ptr(xp);
    ptr->weights = W;
    ptr->init_beta();
  }
}

// [[Rcpp::export]]
void regionModel__set_theta(SEXP xp, std::vector<double> theta, int type = 0) {
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    ptr->covariance.update_parameters_extern(theta);
  } else if (type == 1){
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    double rho = theta.back();
    theta.pop_back();
    ptr->covariance.update_parameters_extern(theta);
    ptr->covariance.update_rho(rho);
  } else  if (type == 2){
    Rcpp::XPtr<rts::regionModel<glmmr::hsgpCovariance>> ptr(xp);
    ptr->covariance.update_parameters_extern(theta);
  } else if (type == 3){
    Rcpp::XPtr<rts::regionModel<glmmr::spdeCovariance>> ptr(xp);
    ptr->covariance.update_parameters_extern(theta);
  }
}

// [[Rcpp::export]]
void regionModel__set_offset(SEXP xp, Eigen::VectorXd offset, int type = 0) {
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    ptr->set_offset(offset);
  } else if (type == 1){
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    ptr->set_offset(offset);
  } else  if (type == 2){
    Rcpp::XPtr<rts::regionModel<glmmr::hsgpCovariance>> ptr(xp);
    ptr->set_offset(offset);
  } else if (type == 3){
    Rcpp::XPtr<rts::regionModel<glmmr::spdeCovariance>> ptr(xp);
    ptr->set_offset(offset);
  }
}

// [[Rcpp::export]]
void regionModel__fit(SEXP xp, double tol, int max_iter, int hist, int k0, int type = 0){
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    ptr->fit(tol, max_iter, hist, k0);
  } else if (type == 1) {
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    ptr->fit(tol, max_iter, hist, k0);
  } else if (type == 2){
    Rcpp::XPtr<rts::regionModel<glmmr::hsgpCovariance>> ptr(xp);
    ptr->fit(tol, max_iter, hist, k0);
  } else if (type == 3) {
    Rcpp::XPtr<rts::regionModel<glmmr::spdeCovariance>> ptr(xp);
    ptr->fit(tol, max_iter, hist, k0);
  }
}

// [[Rcpp::export]]
SEXP regionModel__information_matrix(SEXP xp, bool mc = true, int type = 0){
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->information_matrix(mc));
  } else if (type == 1){
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->information_matrix(mc));
  } else if (type == 2) {
    Rcpp::XPtr<rts::regionModel<glmmr::hsgpCovariance>> ptr(xp);
    return Rcpp::wrap(ptr->information_matrix(mc));
  } else if (type == 3) {
    Rcpp::XPtr<rts::regionModel<glmmr::spdeCovariance>> ptr(xp);
    return Rcpp::wrap(ptr->information_matrix(mc));
  }
}

// [[Rcpp::export]]
SEXP regionModel__information_matrix_theta(SEXP xp, int type = 0){
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->covariance.infomat_theta);
  } else if (type == 1){
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->covariance.infomat_theta);
  } else if (type == 2) {
    Rcpp::XPtr<rts::regionModel<glmmr::hsgpCovariance>> ptr(xp);
    return Rcpp::wrap(ptr->covariance.infomat_theta);
  } else if (type == 3) {
    Rcpp::XPtr<rts::regionModel<glmmr::spdeCovariance>> ptr(xp);
    return Rcpp::wrap(ptr->covariance.infomat_theta);
  }
}

// [[Rcpp::export]]
SEXP regionModel__u(SEXP xp, bool scaled = true, int type = 0){
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->u(scaled));
  } else if (type == 1){
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->u(scaled));
  } else if (type == 2) {
    Rcpp::XPtr<rts::regionModel<glmmr::hsgpCovariance>> ptr(xp);
    return Rcpp::wrap(ptr->u(scaled));
  } else if (type == 3) {
    Rcpp::XPtr<rts::regionModel<glmmr::spdeCovariance>> ptr(xp);
    return Rcpp::wrap(ptr->u(scaled));
  }
}

// [[Rcpp::export]]
SEXP regionModel__get_beta(SEXP xp, int type = 0){
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->beta);
  } else if (type ==1){
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->beta);
  } else if (type == 2) {
    Rcpp::XPtr<rts::regionModel<glmmr::hsgpCovariance>> ptr(xp);
    return Rcpp::wrap(ptr->beta);
  } else if (type == 3) {
    Rcpp::XPtr<rts::regionModel<glmmr::spdeCovariance>> ptr(xp);
    return Rcpp::wrap(ptr->beta);
  }
}

// [[Rcpp::export]]
SEXP regionModel__get_weights(SEXP xp, int type = 0){
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->sampling_weights());
  } else if (type == 1){
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->sampling_weights());
  } else if (type == 2) {
    Rcpp::XPtr<rts::regionModel<glmmr::hsgpCovariance>> ptr(xp);
    return Rcpp::wrap(ptr->sampling_weights());
  } else if (type == 3) {
    Rcpp::XPtr<rts::regionModel<glmmr::spdeCovariance>> ptr(xp);
    return Rcpp::wrap(ptr->sampling_weights());
  }
}

// [[Rcpp::export]]
SEXP regionModel__log_likelihood(SEXP xp, int type = 0){
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->total_log_likelihood());
  } else if (type == 1) {
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->total_log_likelihood());
  } else if (type == 2) {
    Rcpp::XPtr<rts::regionModel<glmmr::hsgpCovariance>> ptr(xp);
    return Rcpp::wrap(ptr->total_log_likelihood());
  } else if (type == 3) {
    Rcpp::XPtr<rts::regionModel<glmmr::spdeCovariance>> ptr(xp);
    return Rcpp::wrap(ptr->total_log_likelihood());
  }
}

// [[Rcpp::export]]
SEXP regionModel__get_theta(SEXP xp, int type = 0){
  if(type == 0) {
    Rcpp::XPtr<rts::regionModel<glmmr::Covariance>> ptr(xp);
    return Rcpp::wrap(ptr->covariance.parameters_);
  } else if (type == 1) {
    Rcpp::XPtr<rts::regionModel<glmmr::ar1Covariance>> ptr(xp);
    std::vector<double> theta = ptr->covariance.parameters_;
    theta.push_back(ptr->covariance.rho);
    return Rcpp::wrap(theta);
  } else if (type == 2) {
    Rcpp::XPtr<rts::regionModel<glmmr::hsgpCovariance>> ptr(xp);
    return Rcpp::wrap(ptr->covariance.get_parameters_extern());
  } else if (type == 3) {
    Rcpp::XPtr<rts::regionModel<glmmr::spdeCovariance>> ptr(xp);
    return Rcpp::wrap(ptr->covariance.get_parameters_extern());
  }
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

