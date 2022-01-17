#' Extract predictions
#'
#' Extract incidence and relative risk predictions
#'
#' @param grid_data sf object describing a regular grid over the area of interest.
#' Should be the same object used in the model fit. See \link[rts2]{create_grid}
#' @param stan_fit A \link[rstan]{stanfit} or \link[cmdstanr]{CmdStanMCMC} object.
#' Output of \link[rts2]{lgcp_fit}
#' @param type Vector of character strings. Any combination of "pred", "rr", and "irr", which are,
#' posterior mean incidence (overall and population standardised), relative risk,
#' and incidence rate ratio, respectively.
#' @param irr.lag integer. If "irr" is requested as `type` then the number of time
#' periods lag previous the ratio is in comparison to
#' @param t.lag integer. Extract predictions for previous time periods.
#' @param popdens character string. Name of the column in `grid_data` with the
#' population density data
#' @return An `sf` object identical to `grid_data` but with extra columns. See Details
#' @details
#' Three outputs can be extracted from the model fit, which will be reported
#' in the columns:
#'
#' Predicted incidence: If type includes `pred` then `pred_mean_total` and
#' `pred_mean_total_sd` provide the
#' predicted mean total incidence and its standard deviation, respectively.
#' `pred_mean_pp` and `pred_mean_pp_sd` provide the predicted population
#' standardised incidence and its standard deviation.
#'
#' Relative risk: if type includes `rr` then the relative risk is reported in
#' the columns `rr` and `rr_sd`. The relative risk here is the exponential
#' of the latent field, which describes the relative difference between
#' expexted mean and predicted mean incidence.
#'
#' Incidence risk ratio: if type includes `irr` then the incidence rate ratio (IRR)
#' is reported in the columns `irr` and `irr_sd`. This is the ratio of the predicted
#' incidence in the last period (minus `t_lag`) to the predicted incidence in the
#' last period minus `irr_lag` (minus `t_lag`). For example, if the time period
#' is in days then setting `irr_lag` to 7 and leaving `t_lag=0` then the IRR
#' is the relative change in incidence in the present period compared to a week
#' prior.
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
#' o1 <- extract_preds(g1,
#'                     res,
#'                     type=c("pred","rr"),
#'                     popdens="cov")
#' }
#' @importFrom stats sd
#' @export
extract_preds <- function(grid_data,
                          stan_fit,
                          type=c("pred","rr","irr"),
                          irr.lag=NULL,
                          t.lag=0,
                          popdens=NULL){

  if("irr"%in%type&is.null(irr.lag))stop("For irr set irr.lag")
  if(!is(grid_data,"sf"))stop("grid_data not sf")
  if(!(is(stan_fit,"CmdStanMCMC")|is(stan_fit,"stanfit")))stop("stan fit required")
  if("pred"%in%type&is.null(popdens))stop("set popdens for pred")

  nCells <- nrow(grid_data)
  if(is(stan_fit,"stanfit")){
    ypred <- rstan::extract(stan_fit,"y_grid_predict")
    ypred <- ypred$y_grid_predict
    f <- rstan::extract(stan_fit,"f")
    f <- f$f
    nT <- dim(ypred)[2]/nCells
    cmdst <- FALSE
  } else if(is(stan_fit,"CmdStanMCMC")){
    if(requireNamespace("cmdstanr")){
      ypred <- stan_fit$draws("y_grid_predict")
      f <- stan_fit$draws("f")
      nT <- dim(ypred)[3]/nCells
      cmdst <- TRUE
    }
  }

  #print(nT)
  #print(nCells)

  if(nT>1){
    if("pred"%in%type){
      if(!cmdst){
        fmu <- ypred[,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells)]/as.data.frame(grid_data)[,popdens]
        grid_data$pred_mean_total <- apply(ypred[,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells)],2,mean)
        grid_data$pred_mean_total_sd <- apply(ypred[,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells)],2,sd)
        grid_data$pred_mean_pp <- apply(fmu,2,mean)
        grid_data$pred_mean_pp_sd <- apply(fmu,2,sd)
      } else {
        fmu <- ypred[,,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells)]/as.data.frame(grid_data)[,popdens]
        grid_data$pred_mean_total <- apply(ypred[,,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells)],3,mean)
        grid_data$pred_mean_total_sd <- apply(ypred[,,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells)],3,sd)
        grid_data$pred_mean_pp <- apply(fmu,3,mean)
        grid_data$pred_mean_pp_sd <- apply(fmu,3,sd)
      }

    }

    if("rr"%in%type){
      if(!cmdst){
        grid_data$rr <- exp(apply(f[,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells)],2,mean))
        grid_data$rr_sd <- exp(apply(f[,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells)],2,sd))
      } else {
        grid_data$rr <- exp(apply(f[,,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells)],3,mean))
        grid_data$rr_sd <- exp(apply(f[,,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells)],3,sd))
      }

    }

    if("irr"%in%type){
      if(!cmdst){
        grid_data$irr <- apply(ypred[,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells)]/
                                 ypred[,((nT-irr.lag-t.lag)*nCells+1):(((nT-t.lag)-irr.lag+1)*nCells)],2,mean)
        grid_data$irr_sd <- apply(ypred[,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells)]/
                                    ypred[,((nT-irr.lag-t.lag)*nCells+1):(((nT-t.lag)-irr.lag+1)*nCells)],2,sd)
      } else {
        grid_data$irr <- apply(ypred[,,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells)]/
                                 ypred[,,((nT-irr.lag-t.lag)*nCells+1):(((nT-t.lag)-irr.lag+1)*nCells)],3,mean)
        grid_data$irr_sd <- apply(ypred[,,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells)]/
                                    ypred[,,((nT-irr.lag-t.lag)*nCells+1):(((nT-t.lag)-irr.lag+1)*nCells)],3,sd)
      }

    }
  } else {

    if("irr"%in%type)stop("cannot estimate irr as only one time period")

    if("pred"%in%type){
      if(!cmdst){
        fmu <- ypred/as.data.frame(grid_data)[,popdens]
        grid_data$pred_mean_total <- apply(ypred,2,mean)
        grid_data$pred_mean_total_sd <- apply(ypred,2,sd)
        grid_data$pred_mean_pp <- apply(fmu,2,mean)
        grid_data$pred_mean_pp_sd <- apply(fmu,2,sd)
      } else {
        fmu <- ypred/as.data.frame(grid_data)[,popdens]
        grid_data$pred_mean_total <- apply(ypred,3,mean)
        grid_data$pred_mean_total_sd <- apply(ypred,3,sd)
        grid_data$pred_mean_pp <- apply(fmu,3,mean)
        grid_data$pred_mean_pp_sd <- apply(fmu,3,sd)
      }

    }

    if("rr"%in%type){
      if(!cmdst){
        grid_data$rr <- exp(apply(f,2,mean))
        grid_data$rr_sd <- exp(apply(f,2,sd))
      } else {
        grid_data$rr <- exp(apply(f,3,mean))
        grid_data$rr_sd <- exp(apply(f,3,sd))
      }

    }


  }


  return(grid_data)
}

#' Hotspots
#'
#' Generate hotspot probabilities
#'
#' Given a definition of a hotspot in terms of threshold(s) for incidence,
#' relative risk, and/or incidence rate ratio, returns the probabilities
#' each area is a "hotspot". See Details of \link[rts2]{extract_preds}
#'
#' @param grid_data sf object describing a regular grid over the area of interest.
#' Should be the same object used in the model fit. See \link[rts2]{create_grid}
#' @param stan_fit A \link[rstan]{stanfit} or \link[cmdstanr]{CmdStanMCMC} object.
#' Output of \link[rts2]{lgcp_fit}
#' @param incidence.threshold Numeric. Threshold of population standardised incidence
#' above which an area is a hotspot
#' @param irr.threshold Numeric. Threshold of incidence rate ratio
#' above which an area is a hotspot.
#' @param irr.lag integer. Lag of time period to calculate the incidence rate ratio.
#' Only required if `irr.threshold` is not `NULL`.
#' @param rr.threshold numeric. Threshold of local relative risk
#' above which an area is a hotspot
#' @param popdens character string. Name of variable in `grid_data`
#' specifying the population density. Needed if `incidence.threshold` is not
#' `NULL`
#' @param col_label character string. If not NULL then the name of the column
#' for the hotspot probabilities.
#' @return An `sf` object identical to `grid_data` with an extract column with the
#' probabilities each polygon is a hotspot
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
#' o1 <- hotspots(g1,
#'                res,
#'                incidence.threshold=1,
#'                popdens="cov")
#' }
#' @export
hotspots <- function(grid_data,
                     stan_fit,
                     incidence.threshold=NULL,
                     irr.threshold=NULL,
                     irr.lag=NULL,
                     rr.threshold=NULL,
                     popdens,
                     col_label=NULL){

  if(all(is.null(incidence.threshold),is.null(irr.threshold),is.null(rr.threshold)))stop("At least one criterion required.")
  if(!is.null(irr.threshold)&is.null(irr.lag))stop("irr.lag must be set")
  if(!is(grid_data,"sf"))stop("grid_data not sf")
  if(!(is(stan_fit,"CmdStanMCMC")|is(stan_fit,"stanfit")))stop("stan fit required")

  nCells <- nrow(grid_data)
  if(is(stan_fit,"stanfit")){
    ypred <- rstan::extract(stan_fit,"y_grid_predict")
    f <- rstan::extract(stan_fit,"f")
  } else if(is(stan_fit,"CmdStanMCMC")){
    if(requireNamespace("cmdstanr")){
      ypred <- stan_fit$draws("y_grid_predict")
      f <- stan_fit$draws("f")
    }
  }
  nT <- dim(ypred)[3]/nCells
  f <- f[,,((nT-1)*nCells+1):(nT*nCells)]

  nCr <- sum(c(!is.null(incidence.threshold),
               !is.null(irr.threshold),
               !is.null(rr.threshold)))


  inc1 <- array(0,dim(f))

  if(!is.null(incidence.threshold)){
    fmu <- ypred[,,((nT-1)*nCells+1):(nT*nCells)]/as.data.frame(grid_data)[,popdens]
    inc1 <- inc1 + I(fmu > incidence.threshold)*1

  }

  if(!is.null(irr.threshold)){
    if(nT==1)stop("cannot estimate irr as only one time period") else {
      inc1 <- inc1 + I(ypred[,,((nT-1)*nCells+1):(nT*nCells)]/
                         ypred[,,((nT-irr.lag)*nCells+1):((nT-irr.lag+1)*nCells)] > irr.threshold)*1
    }

  }

  if(!is.null(rr.threshold)){
    inc1 <- inc1 + I(exp(f) > rr.threshold)*1
  }

  inc1 <- I(inc1 == nCr)*1

  grid_data$hotspot_prob <- apply(inc1,3,mean)

  if(!is.null(col_label)){
    colnames(grid_data)[length(colnames(grid_data))] <- col_label
  }

  return(grid_data)
}

#' Aggregate output
#'
#' Aggregate `lgcp_fit` output to another geography
#'
#' @param grid_data sf object describing a regular grid over the area of interest.
#' Should be the same object used in the model fit. See \link[rts2]{create_grid} and
#' \link[rts2]{lgcp_fit}
#' @param new_geom sf object. A set of polygons covering the same area as `grid_data`
#' @param zcols vector of character strings. Names of the variables in `grid_data` to
#' map to the new geography
#' @param weight_type character string, either "area" or "pop" for area-weighted
#' or population weighted averaging, respectively
#' @param popdens character string. If `weight_type` is equal to "pop" then the
#' name of the column in `grid_data` with population density data
#' @param verbose logical. Whether to provide progress bar.
#' @return An `sf` object identical to `new_geom` with additional columns with the
#' variables specified in `zcols`
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
#' o1 <- extract_preds(g1,
#'                     res,
#'                     type=c("pred","rr"),
#'                     popdens="cov")
#' new1 <- aggregate_output(o1
#'                          cov1,
#'                          zcols="rr")
#' }
#' @importFrom stats weighted.mean
#' @export
aggregate_output <- function(grid_data,
                             new_geom,
                             zcols,
                             weight_type="area",
                             popdens=NULL,
                             verbose=TRUE){

  if(sf::st_crs(grid_data)!=sf::st_crs(new_geom)){
    sf::st_crs(grid_data) <- sf::st_crs(new_geom)
    warning("st_crs(grid_data)!=st_crs(new_geom) setting equal")
  }

  tmp <- sf::st_intersection(grid_data[,zcols],new_geom)
  tmp_len <- lengths(sf::st_intersects(new_geom,grid_data))
  tmp_len <- 1 - tmp_len[1] + cumsum(tmp_len)
  vals <- matrix(nrow=nrow(new_geom),
                 ncol=length(zcols))

  if(verbose)cat("Overlaying geographies\n")

  for(i in 1:nrow(new_geom)){
    if(i < nrow(new_geom)){
      idx_range <- tmp_len[i]:(tmp_len[i+1]-1)
    } else {
      idx_range <- tmp_len[i]:nrow(tmp)
    }

    tmp_range <- tmp[idx_range,]
    w <- as.numeric(sf::st_area(tmp_range))

    if(weight_type=="pop"){
      w <- w*as.data.frame(tmp_range)[,popdens]
    }

    for(j in 1:length(zcols)){
      vals[i,j] <- weighted.mean(as.data.frame(tmp_range)[,zcols[j]],
                                 w=w)
    }

    if(verbose)cat("\r",progress_bar(i,nrow(new_geom)))
  }

  for(j in 1:length(zcols)){
    new_geom$x <- vals[,j]
    colnames(new_geom)[length(colnames(new_geom))] <- zcols[j]
  }

  return(new_geom)

}

#' Generate priors from previous model fit
#'
#' Generate prior distributions from posterior distributions of previous model
#' fit by `lgcp_fit`
#'
#' @param stan_fit A \link[rstan]{stanfit} or \link[cmdstanr]{CmdStanMCMC} object.
#' Output of \link[rts2]{lgcp_fit}
#' @return list
#' @importFrom stats sd
#' @export
extract_priors <- function(stan_fit){
  if(!(is(stan_fit,"CmdStanMCMC")|is(stan_fit,"stanfit")))stop("stan fit required")


  if(is(stan_fit,"stanfit")){
    sig <- rstan::extract(stan_fit,"sigma")
    gam <- rstan::extract(stan_fit,"gamma")
    phi <- rstan::extract(stan_fit,"phi")
  } else if(is(stan_fit,"CmdStanMCMC")){
    if(requireNamespace("cmdstanr")){
      ypred <- stan_fit$draws("y_grid_predict")
      f <- stan_fit$draws("f")
    }
  }

  out <- list(
    prior_lscale=c(mean(phi$phi),sd(phi$phi)),
    prior_var=c(mean(sig$sig),sd(sig$sigma)),
    prior_linpred_mean=apply(gam$gamma,2,mean),
    prior_linpred_sd=apply(gam$gamma,2,sd)
  )

  return(out)
}
