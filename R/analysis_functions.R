#' Fit an approximate log-Gaussian Cox Process model
#'
#' Fit an approximate log-Gaussian Cox Process model
#'
#' @details
#' Our statistical model is a Log Gaussian cox process,
#' whose realisation is observed on the Cartesian area of interest
#' A and time period T. The resulting data are relaisations of an inhomogeneous
#' Poisson process with stochastic intensity function \eqn{\{\lambda{s,t}:s\in A, t \in T\}}.
#' We specify a log-linear model for the intensity:
#'
#' \deqn{\lambda(s,t) = r(s,t)exp(X(s,t)'\gamma + Z(s,t))}
#'
#' where r(s,t) is a spatio-temporally varying Poisson offset.
#' X(s,t) is a length Q vector of covariates including an intercept and
#' Z(s,t) is a latent field. We use an auto-regressive specification for the
#' latent field, with spatial innovation in each field specified as a spatial
#' Gaussian process.
#'
#' We use the fast and accurate approximation for fully Bayesian Gaussian
#' Processes proposed by Solin and Särkkä (1), using basis function
#'  approximations based on approximation
#' via Laplace eigenfunctions for stationary covariance functions.
#' See references (1) and (2) for complete details. The approximation is a linear sum
#' of `m` eigenfunctions with the boundary conditions in each dimension `[-L,L]`.
#' Coordinates in each dimension are scaled to `[-1,1]`, so L represents the
#' proportionate extension of the analysis area.
#'
#' *Priors*
#'The priors should be provided as a list:
#'```
#'priors <- list(
#'   prior_lscale=c(0,0.5),
#'   prior_var=c(0,0.5),
#'   prior_linpred_mean=c(-5,rep(0,7)),
#'   prior_linpred_sd=c(3,rep(1,7))
#' )
#' ```
#' where these refer to the priors:
#' `prior_lscale`: the length scale parameter has a half-normal prior \eqn{N(a,b^2)I[0,\infty)}. The vector is `c(a,b)`.
#' `prior_var`: the standard deviation term has a half normal prior \eqn{\sigma ~ N(a,b^2)I[0,\infty)}. The vector is `c(a,b)`.
#' `prior_linpred_mean` and `prior_linpred_sd`: The parameters of the linear predictor.
#' If X is the nT x Q matrix of covariates, with the first column as ones for the intercept,
#' then the linear prediction contains the term \eqn{X'\gamma}. Each parameter in \eqn{\gamma} has prior
#' \eqn{\gamma_q ~ N(a_q,b_q^2)}.
#' `prior_linpred_mean` should be the vector `(a_1,a_2,...,a_Q)` and
#' `prior_linpred_sd` should be `(b_1,b_2,...,b_Q)`.
#' @references
#' (1) Solin A, Särkkä S. Hilbert space methods for reduced-rank Gaussian
#' process regression. Stat Comput. 2020;30:419–46.
#' doi:10.1007/s11222-019-09886-w.
#'
#' (2) Riutort-Mayol G, Bürkner P-C, Andersen MR, Solin A, Vehtari A.
#' Practical Hilbert space approximate Bayesian Gaussian processes for
#' probabilistic programming. 2020. http://arxiv.org/abs/2004.11408.
#' @param grid_data sf object. A regular grid covering the area of interest (see
#' \link[rts2]{create_grid}). Columns must include `t*`, giving the case
#' count in each time period, as well as any covariates to include in the model
#' (see \link[rts2]{add_covariates}) and the population density
#' @param popdens character string. Name of the population density column
#' @param covs vector of character string. Base names of the covariates to
#' include. For temporally-varying covariates only the stem is required and not
#' the individual column names for each time period (e.g. `dayMon` and not `dayMon1`,
#' `dayMon2`, etc.)
#' @param m integer. Number of basis functions. See Details.
#' @param L integer. Boundary condition as proportionate extension of area, e.g.
#' `L=2` is a doubling of the analysis area. See Details.
#' @param dir character string. Directory to save ouptut.
#' @param iter_warmup integer. Number of warmup iterations
#' @param iter_sampling integer. Number of sampling iterations
#' @param chains integer. Number of chains
#' @param parallel_chains integer. Number of parallel chains
#' @param priors list. See Details
#' @param verbose logical. Provide feedback on progress
#' @param use_cmdstanr logical. Defaults to false. If true then cmdstanr will be used
#' instead of rstan.
#' @param ... additional options to pass to `$sample()``, see \link[cmdstanr]{sample}
#' @return A \link[rstan]{stanfit} or a \link[cmdstanr]{CmdStanMCMC} object
#' @examples
#' \dontrun{
#' b1 <- sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
#' g1 <- create_grid(b1,0.5)
#' dp <- data.frame(y=runif(10,0,3),x=runif(10,0,3),date=paste0("2021-01-",11:20))
#' dp <- create_points(dp,pos_vars = c('y','x'),t_var='date')
#' cov1 <- create_grid(b1,0.8)
#' cov1$cov <- runif(nrow(cov1))
#' g1 <- add_covariates(g1,
#'                      cov1,
#'                      zcols="cov",
#'                     verbose = FALSE)
#' g1 <- points_to_grid(g1, dp, laglength=5)
#' priors <- list(
#'   prior_lscale=c(0,0.5),
#'   prior_var=c(0,0.5),
#'   prior_linpred_mean=c(0),
#'   prior_linpred_sd=c(5)
#'   )
#' res <- lgcp_fit(g1,
#'                 popdens="cov")
#' }
#' @importFrom utils stack capture.output
#' @export
lgcp_fit <- function(grid_data,
                     popdens,
                     covs=NULL,
                     m=10,
                     L=1.5,
                     dir=NULL,
                     iter_warmup=500,
                     iter_sampling=500,
                     chains=3,
                     parallel_chains=3,
                     priors=NULL,
                     verbose=TRUE,
                     use_cmdstanr = FALSE,
                     ...){

  if(verbose)if(is.null(dir))message("dir not set, files will be lost after session restart")
  if(m<1)stop("m must be positive")
  if(m >25)warning("m is very large, sampling may take a very long time.")

    #prepare data for model fit

  ind <- as.matrix(expand.grid(1:m,1:m))
  nT <- sum(grepl("\\bt[0-9]",colnames(grid_data)))
  if(nT==0){
    if("y"%in%colnames(grid_data)){
      nT <- 1
    } else {
      stop("case counts not defined in data")
    }
  }
  nCell <- nrow(grid_data)

  x_grid <- as.data.frame(suppressWarnings(sf::st_coordinates(
    sf::st_centroid(grid_data))))

  popd <- rep(as.data.frame(grid_data)[,popdens],nT)

  # scale to -1,1 in all dimensions
  xrange <- range(x_grid[,1])
  yrange <- range(x_grid[,2])
  std_val <- max(max(xrange - mean(xrange)),max(yrange - mean(yrange)))

  x_grid[,1] <- (x_grid[,1]- mean(xrange))/std_val
  x_grid[,2] <- (x_grid[,2]- mean(yrange))/std_val

  if(nT > 1){
    y <- stack(as.data.frame(grid_data)[,paste0("t",1:nT)])[,1]
  } else {
    y <- as.data.frame(grid_data)[,"y"]
  }


  #add covariates
  if(!is.null(covs)){
    nQ <- length(covs)
    X <- matrix(NA,nrow=length(y),ncol=nQ+1)
    X[,1] <- 1
    for(i in 1:nQ){
      nColV <- sum(grepl(covs[i],colnames(grid_data)))
      if(nColV==1){
        X[,i+1] <- rep(as.data.frame(grid_data)[,covs[i]],nT)
      } else if(nColV==0){
        stop(paste0(covs[i]," not found"))
      } else {
        if(nT>1){
          X[,i+1] <- stack(as.data.frame(grid_data)[,paste0(covs[i],1:nT)])[,1]
        } else {
          X[,i+1] <- as.data.frame(grid_data)[,covs[i]]
        }
      }
      Q <- nQ+1
    }
  } else {
    X <- matrix(1,nrow=length(y),ncol=1)
    Q <- 1
  }

  if(!is.null(priors)){
    if(length(priors$prior_linpred_mean)!=Q|length(priors$prior_linpred_sd)!=Q)
      stop("Prior mean or sd vector for linear predictior is not equal to number of covariates")
    if(length(priors$prior_lscale)!=2|length(priors$prior_var)!=2)
      stop("prior_lscale or prior_var not of length 2")
  } else {
    priors <- list(
      prior_lscale = c(0,0.5),
      prior_var = c(0,0.5),
      prior_linpred_mean = rep(0,Q),
      prior_linpred_sd = rep(5,Q)
    )
  }





  if(verbose)message(paste0(nCell," grid cells ",nT," time periods, and ",Q," covariates. Starting sampling..."))

  datlist <- list(
    D = 2,
    Q = Q,
    L = c(L,L),
    M = m,
    M_nD = m^2,
    nT= nT,
    Nsample = nCell,
    y = y,
    x_grid = x_grid[,1:2],
    indices = ind,
    popdens = popd,
    X= X,
    prior_lscale=priors$prior_lscale,
    prior_var=priors$prior_var,
    prior_linpred_mean = as.array(priors$prior_linpred_mean),
    prior_linpred_sd=as.array(priors$prior_linpred_sd)
  )

  if(use_cmdstanr){
    if(!requireNamespace("cmdstanr"))stop("cmdstanr not available.")
    model_file <- system.file("stan",
                              "approxlgcp.stan",
                              package = "rts2",
                              mustWork = TRUE)
    model <- cmdstanr::cmdstan_model(model_file)
    res <- model$sample(
      data=datlist,
      chains=chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      parallel_chains = parallel_chains,
      output_dir=dir,
      init=0.5,
      ...
    )
  } else {
    if(!verbose){
      capture.output(suppressWarnings( res <- rstan::sampling(stanmodels$approxlgcp,
                             data=datlist,
                             chains=chains,
                             iter = iter_warmup+iter_sampling,
                             warmup = iter_warmup,
                             cores = parallel_chains,
                             refresh = 0)), file=tempfile())
    } else {
      res <- rstan::sampling(stanmodels$approxlgcp,
                             data=datlist,
                             chains=chains,
                             iter = iter_warmup+iter_sampling,
                             warmup = iter_warmup,
                             cores = parallel_chains)
    }

  }


  #save output
  return(res)
}


#' Returns scale conversion factor
#'
#' Coordinates are scaled to `[-1,1]` for \link[rts2]{lgcp_fit}. This function
#' returns the scaling factor for this conversion.
#'
#' @param grid_data sf object. See \link[rts2]{create_grid}
#' @return numeric
#' @examples
#' b1 = sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
#' g1 <- create_grid(b1,0.5)
#' scale_conversion_factor(g1)
#' @export
scale_conversion_factor <- function(grid_data){
  x_grid <- as.data.frame(suppressWarnings(sf::st_coordinates(sf::st_centroid(grid_data))))
  xrange <- range(x_grid[,1])
  yrange <- range(x_grid[,2])
  std_val <- max(max(xrange - mean(xrange)),max(yrange - mean(yrange)))
  return(std_val)
}
