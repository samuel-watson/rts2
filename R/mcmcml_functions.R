#' MCMCML Algorithm for the LGCP
#' 
#' Processes data and runs Markov Chain Monte Carlo Maximum Likelihood algorithms to fit
#' Log Gaussian Cox Process (LGCP) model. Used internally. Runs NNGP and full LGCP and both
#' standard and "region" models. 
#' @param y Vector of outcome counts 
#' @param X Matrix of covariates for the linear predictor
#' @param coords A data frame providing the x and y coordinates of the computational grid
#' @param popdens Vector with the log population density
#' @param nT Number of time periods
#' @param region_data Optional. A list with the named elements `ncell`, `cell_id`, and `q_weights` giving, respectively
#' the number of cells overlapping each region, the indexes of the cells overlapping each region in order, and the weights
#' as the proportion of the area of the cell in each overlapping cell in order.
#' @param start Starting values of the model parameteters in the order c(beta, theta, rho). Rho is only for models with nT>1
#' @param mod Either "exp" for exponential or "sqexp" for squared exponential covariance function
#' @param use_cmdstanr Logical indicating whether to use cmdstanr
#' @param mcml_options List of options for the algorithm. See details.
#' @param mcmc_options List of options for the MCMC sampling. See details.
#' @param verbose Logical indicating whether to provide detailed feedback.
#' @details 
#' The argument `mcml_options` is a named list with the options of whether to use NNGP (`useNN`), whether to use 
#' Newton-Raphson (`mcnr`), whether to treat the covariance parameters as known (`known_theta`), whether to provide
#' more detailed output (`trace`), the number of nearest neighbours if using NNGP (`nNN`), the tolerance for when to 
#' terminate the algorithm (`tol`), and the maximum number of algorithm iterations (`maxiter`). The argument 
#' `mcmc_options` is a named list with the number of warmup and sampling iterations for the MCMC sampler.
#' @return A named list with the estimated parameters, number of iterations, whether the algorithm converged, and 
#' the random effect samples.
lgcp_mcmcml <- function(y,X,coords,
                        popdens,
                        nT,
                        region_data,
                        start,
                        mod = "exp",
                        use_cmdstanr = TRUE,
                        mcml_options = list(useNN=FALSE,mcnr=FALSE,known_theta=FALSE,trace=1,nNN=10,tol=1e-2, maxiter=10),
                        mcmc_options = list(warmup = 100, sampling = 100),
                        verbose = TRUE
                        ){
  nCell <- nrow(coords)
  useRegion <- !missing(region_data)
  
  beta <- start[1:ncol(X)]
  if(nT == 1){
    theta <- start[(ncol(X)+1):length(start)]
    rho <- 1
  } else {
    theta <- start[(ncol(X)+1):(length(start)-1)]
    rho <- start[length(start)]
  }
  newbeta <- beta
  newtheta <- theta
  newrho <- rho
  diff <- rep(1,length(start))
  
  f1 <- ifelse(mod=="exp","~(1|fexp(X,Y))","~(1|sqexp(X,Y))")
  
  cov <- glmmrBase::Covariance$new(
    formula =  formula(f1),
    parameters = theta,
    data=coords
  )
  
  ddata <- cov$get_D_data()
  
  if(mcml_options$useNN){
    NN <- genNN(as.matrix(coords),mcml_options$nNN)
    AD <- get_AD(cov = ddata$cov,data = ddata$data,NN = NN-1,theta = theta)
    L <- inv_ldlt(AD$A,AD$D,NN) # cholesky decomposition
    ZL <- get_ZL(L,nT,rho)
  } else {
    L <- cov$get_chol_D()
    ZL <- get_ZL(L,nT,rho)
  }
  
  ## MCMC sampling
  ## create data
  xb <- X%*%beta + popdens
  
  if(useRegion){
    dat <- list(
      N = nCell,
      nT = nT,
      nRegion = length(region_data$ncell)-1,
      n_Q = length(region_data$q_weights),
      Xb = drop(xb),
      ZL = ZL,
      y = y,
      n_cell = region_data$ncell,
      cell_id = region_data$cell_id,
      q_weights = region_data$q_weights
    )
    filen <- "mcml_poisson_region_cmd.stan"
    filer <- "mcml_poisson_region"
    filetmp <- "D:/Documents/R/rts2/inst/cmdstan/mcml_poisson_region_cmd.stan"
  } else {
    dat <- list(
      N = nrow(ZL),
      Xb = drop(xb),
      Z = ZL,
      y = y
    )
    filen <- "mcml_poisson_cmd.stan"
    filer <- "mcml_poisson"
    filetmp <- "D:/Documents/R/rts2/inst/cmdstan/mcml_poisson_cmd.stan"
  }
  
  ## function name
  str <- ifelse(mcml_options$useNN,"nngp_","lgcp_")
  if(useRegion)str <- paste0(str,"region_")
  str <- paste0(str,"optim")
  #if(mcml_options$known_theta)str <- paste0(str,"_known_cov")
  
  ## build arguments list
  args <- list(cov = matrix(ddata$cov,nrow=1),
               data = ddata$data,
               X = matrix(1,nrow=length(y),ncol=1),
               y = y)
  if(mcml_options$useNN){
    args <- append(args,list(NN= NN-1))
  }
  args <- append(args,list(u = matrix(0,nrow=nrow(L),ncol=mcmc_options$sampling)))
  start_par <- beta
  start_par <- c(start_par,theta)
  if(nT>1)start_par <- c(start_par,rho)
  args <- append(args,list(start= start_par))
  args <- append(args,list(offset = popdens))
  if(useRegion){
    args <- append(args,list(n_cell = region_data$ncell,
                             cell_id = region_data$cell_id,
                             q_weights = region_data$q_weights))
  }
  args <- append(args,list(nT = nT,
                           trace = mcml_options$trace,
                           mcnr = mcml_options$mcnr,
                           known_cov = mcml_options$known_theta))
  
  if(use_cmdstanr){
    if(!requireNamespace("cmdstanr"))stop("cmdstanr not available.")
    # model_file <- system.file("cmdstan",
    #                           filen,
    #                           package = "rts2",
    #                           mustWork = TRUE)
    model <- cmdstanr::cmdstan_model(filetmp)
  }
  
  iter <- 1
  converged <- FALSE
  maxdiff <- 1
  
  while(iter<=mcml_options$maxiter & !converged){
    if(verbose)cat(Reduce(paste0,rep("-",30)),"\nIter: ",iter,"\n")
    
    xb <- X%*%beta + popdens
    dat$Xb <- drop(xb)
    
    if(!mcml_options$known_theta){
      if(mcml_options$useNN){
        AD <- get_AD(cov = ddata$cov,data = ddata$data,NN = NN-1,theta = theta)
        L <- inv_ldlt(AD$A,AD$D,NN) # cholesky decomposition
        ZL <- get_ZL(L,nT,rho)
      } else {
        cov$update_parameters(theta)
        L <- cov$get_chol_D()
        ZL <- get_ZL(L,nT,rho)
      }
      dat$ZL <- ZL
    }
    
    start_par <- beta
    start_par <- c(start_par,theta)
    if(nT>1)start_par <- c(start_par,rho)
    args$start <- start_par
    
    if(use_cmdstanr){
      if(verbose){
        res <- model$sample(
          data=dat,
          chains=1,
          iter_warmup = mcmc_options$warmup,
          iter_sampling = mcmc_options$sampling,
          parallel_chains = 1,
        )
      } else {
        capture.output(suppressWarnings(res <- model$sample(
          data=dat,
          chains=1,
          iter_warmup = mcmc_options$warmup,
          iter_sampling = mcmc_options$sampling,
          parallel_chains = 1,
        )),file = tempfile())
      }
      
      
      u <- res$draws("gamma",format = "matrix")
      u <- t(u)
      dsamps <- get_Lu(L,u)
    } else {
      if(verbose){
        res <- rstan::sampling(stanmodels[[filer]],
                               data=dat,
                               chains=1,
                               iter = mcmc_options$warmup+mcmc_options$sampling,
                               warmup = mcmc_options$warmup,
                               cores = 1)
      } else {
        capture.output(
          suppressWarnings(
            res <- rstan::sampling(stanmodels[[filer]],
                                   data=dat,
                                   chains=1,
                                   iter = mcmc_options$warmup+mcmc_options$sampling,
                                   warmup = mcmc_options$warmup,
                                   cores = 1)
          ), file = tempfile()
        )
      }
      
      u <- rstan::extract(res,"gamma")
      u <- t(u$gamma)
      dsamps <- get_Lu(L,u)
    }
    
    args$u <- as.matrix(dsamps)
    if(mcml_options$trace > 0)cat("\nStarting parameter estimation...\n")
    out <- do.call(str,args)
    
    newbeta <- out$beta
    if(!mcml_options$known_theta)newtheta <- out$theta
    if(nT > 1)newrho <- out$rho
    diff <- c(abs(beta-newbeta),abs(theta-newtheta),abs(rho-newrho))
    maxdiff <- max(diff)
    converged <- maxdiff <= mcml_options$tol
    
    beta <- newbeta
    theta <- newtheta
    rho <- newrho
    iter <- iter + 1
    
    if(verbose){
      cat("\nbeta: ",beta,"\ntheta: ",theta)
      if(nT>1)cat("\nrho: ",rho)
      cat("\n",Reduce(paste0,rep("-",30)))
      if(converged)cat("\nCONVERGED!")
    }
  }
  
  
  # get standard errors
  
  if(!mcml_options$known_theta){
    if(mcml_options$useNN){
      AD <- get_AD(cov = ddata$cov,data = ddata$data,NN = NN-1,theta = theta)
      L <- inv_ldlt(AD$A,AD$D,NN) # cholesky decomposition
      ZL <- get_ZL(L,nT,rho)
    } else {
      cov$update_parameters(theta)
      L <- cov$get_chol_D()
      ZL <- get_ZL(L,nT,rho)
    }
  }
  S <- (ZL%*%t(ZL))
  XSX <- t(X)%*%solve(S)%*%X
  XSX <- solve(XSX)
  
  return(list(beta = beta, theta=theta, 
              rho=rho, iter = iter, 
              converged = converged, u = dsamps,
              se = XSX))
}

#' Laplace approximation Algorithm for the LGCP
#' 
#' Processes data and runs a Laplacian Approximation Maximum Likelihood algorithm to fit
#' Log Gaussian Cox Process (LGCP) model. Used internally. Runs NNGP and full LGCP for point data models.
#' @param y Vector of outcome counts 
#' @param X Matrix of covariates for the linear predictor
#' @param coords A data frame providing the x and y coordinates of the computational grid
#' @param popdens Vector with the log population density
#' @param nT Number of time periods
#' @param start Starting values of the model parameteters in the order c(beta, theta, rho). Rho is only for models with nT>1
#' @param mod Either "exp" for exponential or "sqexp" for squared exponential covariance function
#' @param la_options List of options for the algorithm. See details.
#' @param verbose Logical indicating whether to provide detailed feedback.
#' @details 
#' The argument `la_options` is a named list with the options of whether to use NNGP (`useNN`), whether to treat the covariance parameters as known (`known_theta`), whether to provide
#' more detailed output (`trace`), the number of nearest neighbours if using NNGP (`nNN`), the tolerance for when to 
#' terminate the algorithm (`tol`), and the maximum number of algorithm iterations (`maxiter`). 
#' @return A named list with the estimated parameters, number of iterations, whether the algorithm converged, and 
#' the random effect samples.
lgcp_la <- function(y,X,coords,
                        popdens,
                        nT,
                        start,
                      mod = "exp",
                        la_options = list(useNN=FALSE,known_theta=FALSE,trace=1,nNN=10,tol=1e-2, maxiter=10),
                        verbose = TRUE
){
  nCell <- nrow(coords)
  
  beta <- start[1:ncol(X)]
  if(nT == 1){
    theta <- start[(ncol(X)+1):length(start)]
    rho <- 1
  } else {
    theta <- start[(ncol(X)+1):(length(start)-1)]
    rho <- start[length(start)]
  }
  newbeta <- beta
  newtheta <- theta
  newrho <- rho
  diff <- rep(1,length(start))
  
  f1 <- ifelse(mod=="exp","~(1|fexp(X,Y))","~(1|sqexp(X,Y))")
  
  cov <- glmmrBase::Covariance$new(
    formula =  formula(f1),
    parameters = theta,
    data=coords
  )
  
  ddata <- cov$get_D_data()
  
  ## function name
  str <- ifelse(la_options$useNN,"nngp_optim_la","lgcp_optim_la")
  #if(la_options$known_theta)str <- paste0(str,"_known_cov")
  
  ## build arguments list
  args <- list(cov = matrix(ddata$cov,nrow=1),
               data = ddata$data,
               X = matrix(1,nrow=length(y),ncol=1),
               y = y)
  if(la_options$useNN){
    NN <- genNN(as.matrix(coords),la_options$nNN)
    args <- append(args,list(NN= NN-1))
  }
  start_par <- beta
  start_par <- c(start_par,theta)
  if(nT>1)start_par <- c(start_par,rho)
  args <- append(args,list(start= start_par))
  args <- append(args,list(offset = popdens,
                           nT = nT,
                           known_cov = la_options$known_theta,
                           trace = la_options$trace))
  
  out <- do.call(str,args)
  
  return(list(beta = out$beta, theta=out$theta, rho=out$rho, u = out$u))
}

#' Generates samples of the random effects
#' 
#' Using known model parameters generates samples of the random effects from either the standard or regional model
#' @param y Vector of outcome counts 
#' @param xb Vector of values of the linear predictor including any offsets
#' @param coords A data frame providing the x and y coordinates of the computational grid
#' @param nT Number of time periods
#' @param region_data Optional. A list with the named elements `ncell`, `cell_id`, and `q_weights` giving, respectively
#' the number of cells overlapping each region, the indexes of the cells overlapping each region in order, and the weights
#' as the proportion of the area of the cell in each overlapping cell in order.
#' @param start Values of the covariance parameters in the order c(theta, rho). Rho is only for models with nT>1
#' @param mod Either "exp" for exponential or "sqexp" for squared exponential covariance function
#' @param use_cmdstanr Logical indicating whether to use cmdstanr
#' @param mcmc_options List of options for the MCMC sampling, specifically the number of warmup and sampling iterations.
#' @return A matrix of samples
#' @export
sample_u <- function(y,xb,coords,
                        nT,
                        region_data,
                        start,
                        mod = "exp",
                        use_cmdstanr = TRUE,
                     verbose = TRUE,
                        mcmc_options = list(warmup = 100, sampling = 100)
){
  nCell <- nrow(coords)
  useRegion <- !missing(region_data)
  
  if(nT == 1){
    theta <- start
    rho <- 1
  } else {
    theta <- start[1:(length(start)-1)]
    rho <- start[length(start)]
  }
  
  f1 <- ifelse(mod=="exp","~(1|fexp(X,Y))","~(1|sqexp(X,Y))")
  
  cov <- glmmrBase::Covariance$new(
    formula =  formula(f1),
    parameters = theta,
    data=coords
  )
  
  ddata <- cov$get_D_data()
  
  L <- cov$get_chol_D()
  ZL <- get_ZL(L,nT,rho)
  
  if(useRegion){
    dat <- list(
      N = nCell,
      nT = nT,
      nRegion = length(region_data$ncell)-1,
      n_Q = length(region_data$q_weights),
      Xb = drop(xb),
      ZL = ZL,
      y = y,
      n_cell = region_data$ncell,
      cell_id = region_data$cell_id,
      q_weights = region_data$q_weights
    )
    filen <- "mcml_poisson_region_cmd.stan"
    filer <- "mcml_poisson_region"
    filetmp <- "D:/Documents/R/rts2/inst/cmdstan/mcml_poisson_region_cmd.stan"
  } else {
    dat <- list(
      N = nrow(ZL),
      Xb = drop(xb),
      Z = ZL,
      y = y
    )
    filen <- "mcml_poisson_cmd.stan"
    filer <- "mcml_poisson"
    filetmp <- "D:/Documents/R/rts2/inst/cmdstan/mcml_poisson_cmd.stan"
  }
  
  
  if(use_cmdstanr){
    # model_file <- system.file("cmdstan",
    #                           filen,
    #                           package = "rts2",
    #                           mustWork = TRUE)
    model <- cmdstanr::cmdstan_model(filetmp)
    if(verbose){
      res <- model$sample(
        data=dat,
        chains=1,
        iter_warmup = mcmc_options$warmup,
        iter_sampling = mcmc_options$sampling,
        parallel_chains = 1,
      )
    } else {
      capture.output(
        suppressWarnings(
          res <- model$sample(
            data=dat,
            chains=1,
            iter_warmup = mcmc_options$warmup,
            iter_sampling = mcmc_options$sampling,
            parallel_chains = 1,
          )
        ), file = tempfile()
      )
    }
    
    
    u <- res$draws("gamma",format = "matrix")
    u <- t(u)
  } else {
    if(verbose){
      res <- rstan::sampling(stanmodels[[filer]],
                             data=dat,
                             chains=1,
                             iter = mcmc_options$warmup+mcmc_options$sampling,
                             warmup = mcmc_options$warmup,
                             cores = 1)
    } else {
      capture.output(
        suppressWarnings(
          res <- rstan::sampling(stanmodels[[filer]],
                                 data=dat,
                                 chains=1,
                                 iter = mcmc_options$warmup+mcmc_options$sampling,
                                 warmup = mcmc_options$warmup,
                                 cores = 1)
        ), file = tempfile()
      )
    }
    
    u <- rstan::extract(res,"gamma")
    u <- u$gamma
  }
  
  return(u)
}