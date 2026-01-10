disaggregate_positive <- function(W, y_r, pop, adj, lambda_smooth = 1, 
                                lambda_entropy = 0, weights = NULL) {
  # Pycnophylactic disaggregation with spatial smoothness on rates
  #
  
  # Solves: min ||W y_s - y_r||^2 
  #             + lambda_smooth * sum_{i~j} w_ij (r_i - r_j)^2
  #             + lambda_entropy * sum(y_s * log(y_s))
  #         s.t. y_s >= 0
  #
  # W:             m x n aggregation matrix
  
  # y_r:           length m vector (observed regional counts)
  # pop:           length n vector (population per cell)
  # adj:           n x n adjacency matrix (1 if neighbours, 0 otherwise)
  #                or a sparse Matrix from the Matrix package
  # lambda_smooth: weight on spatial smoothness
  # lambda_entropy: weight on entropy (0 = pure QP, >0 = nonlinear)
  # weights:       optional n x n matrix of edge weights (e.g., inverse distance)
  
  n <- ncol(W)
  eps <- 1e-10
  
  
  # Build the graph Laplacian for the adjacency structure
  if (is.null(weights)) {
    A <- as.matrix(adj)
  } else {
    A <- as.matrix(adj) * as.matrix(weights)
  }
  D <- diag(rowSums(A))
  L <- D - A  # Graph Laplacian
  
  
  # We want to smooth rates r = y/p, so penalty is r'Lr = y' P^{-1} L P^{-1} y
  # where P = diag(pop)
  P_inv <- diag(1 / pop)
  Q_smooth <- P_inv %*% L %*% P_inv
  
  # Precompute
  WtW <- crossprod(W)
  Wty <- as.vector(crossprod(W, y_r))
  
  if (lambda_entropy == 0) {
    # Pure QP solution
    D_qp <- 2 * (WtW + lambda_smooth * Q_smooth)
    d_qp <- 2 * Wty
    
    # Small ridge for numerical stability
    D_qp <- D_qp + 1e-8 * diag(n)
    
    sol <- solve.QP(Dmat = D_qp, dvec = d_qp, Amat = diag(n), bvec = rep(0, n))
    y_s <- sol$solution
    convergence <- 0
    
  } else {
    # Nonlinear optimisation with entropy term
    obj <- function(y) {
      resid <- as.vector(W %*% y) - y_r
      ls_term <- sum(resid^2)
      smooth_term <- as.numeric(t(y) %*% Q_smooth %*% y)
      entropy_term <- sum(y * log(y + eps))
      ls_term + lambda_smooth * smooth_term + lambda_entropy * entropy_term
    }
    
    grad <- function(y) {
      ls_grad <- 2 * (WtW %*% y - Wty)
      smooth_grad <- 2 * Q_smooth %*% y
      entropy_grad <- log(y + eps) + 1
      as.vector(ls_grad + lambda_smooth * smooth_grad + lambda_entropy * entropy_grad)
    }
    
    start <- rep(sum(y_r) / n, n)
    
    result <- optim(
      par = start, fn = obj, gr = grad,
      method = "L-BFGS-B",
      lower = rep(eps, n),
      control = list(maxit = 1000)
    )
    
    y_s <- result$par
    convergence <- result$convergence
  }
  
  names(y_s) <- colnames(W)
  rates <- y_s / pop
  
  # Compute Moran's I as a diagnostic for spatial autocorrelation
  rate_centered <- rates - mean(rates)
  moran_num <- as.numeric(t(rate_centered) %*% A %*% rate_centered)
  moran_denom <- sum(rate_centered^2)
  moran_I <- (n / sum(A)) * moran_num / moran_denom
  
  list(
    y_s = y_s,
    rates = rates,
    residual = sqrt(sum((W %*% y_s - y_r)^2)),
    rate_variance = var(rates),
    moran_I = moran_I,  # Spatial autocorrelation of fitted rates
    convergence = convergence
  )
}

disaggregate_covariate <- function(x_r, W, adj, pop = NULL, 
                                   lambda_smooth = 1, lambda_ridge = 1e-4) {
  # Disaggregate regional means to fine grid, preserving aggregate means
  #
  
  # x_r:           length m vector of regional means (or pop-weighted means)
  # W:             m x n aggregation matrix (1 if cell in region, 0 otherwise)
  # adj:           n x n adjacency matrix
  # pop:           optional length n population vector; if provided, preserves
  #                population-weighted means rather than simple means
  # lambda_smooth: weight on spatial smoothness
  # lambda_ridge:  small ridge for numerical stability (needed because L is singular)
  
  n <- ncol(W)
  m <- nrow(W)
  
  # Build the mean-computing matrix
  if (is.null(pop)) {
    # Simple mean: divide by count of cells per region
    n_per_region <- rowSums(W)
    W_mean <- W / n_per_region
  } else {
    # Population-weighted mean: weight by population share within region
    pop_per_region <- as.vector(W %*% pop)
    W_mean <- W * rep(pop, each = m) / pop_per_region
  }
  
  # Graph Laplacian for spatial smoothness
  A <- as.matrix(adj)
  D_graph <- diag(rowSums(A))
  L <- D_graph - A
  
  # Regularization: smooth + small ridge (ridge needed for invertibility)
  Q <- lambda_smooth * L + lambda_ridge * diag(n)
  
  # Analytical solution via Lagrange multipliers
  # We minimise x'Qx subject to W_mean x = x_r
  #
  # Solution: x = Q^{-1} W_mean' (W_mean Q^{-1} W_mean')^{-1} x_r
  
  Q_inv <- solve(Q)
  WQinvWt <- W_mean %*% Q_inv %*% t(W_mean)
  
  x_s <- as.vector(Q_inv %*% t(W_mean) %*% solve(WQinvWt, x_r))
  
  names(x_s) <- colnames(W)
  
  # Diagnostics
  fitted_means <- as.vector(W_mean %*% x_s)
  
  # Moran's I for spatial autocorrelation
  x_centered <- x_s - mean(x_s)
  moran_num <- as.numeric(t(x_centered) %*% A %*% x_centered)
  moran_denom <- sum(x_centered^2)
  moran_I <- (n / sum(A)) * moran_num / moran_denom
  
  list(
    x_s = x_s,
    fitted_means = fitted_means,
    target_means = x_r,
    mean_error = max(abs(fitted_means - x_r)),  # Should be ~0
    variance = var(x_s),
    moran_I = moran_I
  )
}

flat_disaggregate <- function(x_r, W) {
  W_colsums <- colSums(W)
  x_s <- as.vector(t(W) %*% x_r / pmax(W_colsums, 1e-10))
  x_s
}
