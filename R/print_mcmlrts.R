#' Extracts the fixed effect estimates
#' 
#' Extracts the fixed effect estimates from an `rtsFit object` returned from call of `lgcp_ml()` or `lgcp_bayes()` in the \link[rts2]{grid} class.
#' @param object An `mcml` model fit.
#' @return A named, numeric vector of fixed-effects estimates.
#' @export
fixed.effects <- function(object){
  fixed <- object$coefficients$est[1:object$P]
  names(fixed) <- object$coefficients$par[1:object$P]
  return(fixed)
}

#' Extracts the random effect estimates
#' 
#' Extracts the random effect estimates or samples from an `rtsFit object` returned from call of `lgcp_ml()` or `lgcp_bayes()` in the \link[rts2]{grid} class.
#' @param object An `mcml` model fit.
#' @return A matrix of dimension (number of fixed effects ) x (number of MCMC samples). For Laplace approximation, the number of "samples" equals one.
#' @export
random.effects <- function(object){
  return(object$re.samps)
}

#' Extracts the estimates of the covariance parameters
#' 
#' Extracts the estimates of the covariance parameters an `rtsFit object` returned from call of `lgcp_ml()` or `lgcp_bayes()` in the \link[rts2]{grid} class.
#' @param object An `mcml` model fit.
#' @return A matrix of dimension (number of fixed effects ) x (number of MCMC samples). For Laplace approximation, the number of "samples" equals one.
#' @export
covariance.parameters <- function(object){
  cov.pars <- object$coefficients$est[(object$P+1):nrow(object$coefficients)]
  names(cov.pars) <- object$coefficients$par[(object$P+1):nrow(object$coefficients)]
  return(cov.pars)
}

#' Prints an rtsFit fit output
#' 
#' Print method for class "`rtsFit`"
#' 
#' @param x an object of class "`rtsFit`" 
#' @param ... Further arguments passed from other methods
#' @details 
#' `print.rtsFit` tries to replicate the output of other regression functions, such
#' as `lm` and `lmer` reporting parameters, standard errors, and z- and p- statistics for maximum 
#' likelihood esitmates, or posterior means, standard deviations and credible intervals for 
#' Bayesian models.
#' @return No return value, called for side effects.
#' @method print rtsFit
#' @export
print.rtsFit <- function(x, ...){
  cat("\nAn rts model fit\n")
  digits <- 4
  ml <- !x$method%in%c("mcmc","vb")
  cat(ifelse(!ml,
             "Bayesian Log Gaussian Cox Process Model\n",
             "Maximum Likelihood Log Gaussian Cox Process Model\n"))
  algo <- switch(as.character(x$method),
                 "1" = "MCMC Maximum Likelihood Expectation Maximisation",
                 "2" = "MCMC Maximum Likelihood Expectation Maximisation",
                 "3" = "MCMC Maximum Likelihood Expectation Maximisation",
                 "4" = "Stochastic Approximation Expectation Maximisation",
                 "5" = "Stochastic Approximation Expectation Maximisation with Ruppert-Polyak Averaging",
                 "6" = "MCMC Maximum Likelihood Expectation Maximisation with Adaptive Sample Sizes",
                 "7" = "MCMC Maximum Likelihood Expectation Maximisation with Adaptive Sample Sizes",
                 "8" = "MCMC Maximum Likelihood Expectation Maximisation with Adaptive Sample Sizes",
                 "mcmc" = "Markov Chain Monte Carlo",
                 "vb" = "Variational Bayes"
                 )
  cat("Using: ",algo)
  if(ml){
    if(x$method %in% c(1:3,6:8))
    cat("\nNumber of Monte Carlo simulations per iteration: ",x$m," with tolerance ",x$tol,"\n\n")
  } else {
    cat("\nMCMC sample size: ",x$m," with ",x$m/x$iter," chains")
  }
  approx <- switch(x$approx,
                   "none" = "None",
                   "lgcp" = "None",
                   "hsgp" = "Hilbert Space Gaussian Process",
                   "nngp" = "Nearest Neighbour Gaussian Process")
  cat("\nApproximation: ",approx)
  cat("\nFixed effects:",fixed.effects(x))
  cat("\nRandom effects:", covariance.parameters(x))
  cat("\n")
}

#' Summary method for class "rtsFit"
#' 
#' Summary method for class "`rtsFit`"
#' 
#' @param object an object of class "`rtsFit`" as a result of a call to `lgcp_ml()` or `lgcp_bayes()`
#' @param ... Further arguments passed from other methods
#' @details 
#' The summary methods aims to replicate the output of other regression model fitting functions and reports 
#' central point estimates, relevant test statistics, and uncertainty intervals. In addition, the returned 
#' summary object will also include time period specific relative risk and incidence predictions.
#' @return An rtsFitSummary object
#' @method summary rtsFit
#' @export
summary.rtsFit <- function(object,...){
  ml <- !object$method%in%c("mcmc","vb")
  algo <- switch(as.character(object$method),
                 "1" = "MCMC Maximum Likelihood Expectation Maximisation",
                 "2" = "MCMC Maximum Likelihood Expectation Maximisation",
                 "3" = "MCMC Maximum Likelihood Expectation Maximisation",
                 "4" = "Stochastic Approximation Expectation Maximisation",
                 "5" = "Stochastic Approximation Expectation Maximisation with Ruppert-Polyak Averaging",
                 "6" = "MCMC Maximum Likelihood Expectation Maximisation with Adaptive Sample Sizes",
                 "7" = "MCMC Maximum Likelihood Expectation Maximisation with Adaptive Sample Sizes",
                 "8" = "MCMC Maximum Likelihood Expectation Maximisation with Adaptive Sample Sizes",
                 "mcmc" = "Markov Chain Monte Carlo",
                 "vb" = "Variational Bayes"
  )
  approx <- switch(object$approx,
                   "none" = "None",
                   "hsgp" = "Hilbert Space Gaussian Process",
                   "nngp" = "Nearest Neighbour Gaussian Process")
  if(ml){
    colrange <- 2:7
  } else {
    colrange <- c(2,3,6,7)
  }
  pars <- object$coefficients[,colrange]
  if(ml){
    colnames(pars) <- c("Estimate","Std. Err.","z value","p value","2.5% CI","97.5% CI")
  } else {
    colnames(pars) <- c("Posterior mean","Posterior Std. Dev.","2.5% CrI","97.5% CrI")
  }
  rownames(pars) <- object$coefficients$par
  pars <- apply(pars,2,round,digits = 4)
  
  preds <- list()
  ncell <- nrow(object$re.samps)/object$nT
  ny <- nrow(object$y_predicted)/object$nT
  for(t in 1:object$nT){
    preds <- append(preds, 
                    list(
                      list(
                        rr = object$re.samps[((t-1)*ncell + 1):(t*ncell),],
                        y = object$y_predicted[((t-1)*ny + 1):(t*ny),]
                      )
                    ))
  }
  summout <- list(
    algo = algo,
    ml = ml,
    m = object$m,
    tol = object$tol,
    conv.criterion = object$conv_criterion,
    iter = object$iter,
    approx = approx,
    fixef = pars[1:object$P,],
    covpars = pars[(object$P+1):nrow(pars),1],
    P = object$P,
    Rsq = object$Rsq,
    aic = object$aic,
    logl = object$logl,
    preds = preds
  )
  class(summout) <- "rtsFitSummary"
  return(summout)
}

#' Prints an rtsFitSummary fit output
#' 
#' Print method for class "`rtsFitSummary`"
#' 
#' @param x an object of class "`rtsFitSummary`" 
#' @param ... Further arguments passed from other methods
#' @details 
#' `print.rtsFitSummary` prints the summary of an rtsFit, see \link[rts2]{summary.rtsFit}
#' @return No return value, called for side effects.
#' @method print rtsFitSummary
#' @export
print.rtsFitSummary <- function(x,...){
  cat("\nAn rts model fit summary\n")
  digits <- 4
  cat(ifelse(!x$ml,
             "Bayesian Log Gaussian Cox Process Model\n",
             "Maximum Likelihood Log Gaussian Cox Process Model\n"))
  cat("\nUsing: ",x$algo)
  if(x$ml){
    if(grepl("MCMC Maximum Likelihood Expectation Maximisation",x$algo)){
      cat("\nNumber of Monte Carlo simulations per iteration: ",x$m," with tolerance ",x$tol,"\n\n")
    } else {
      cat("\nFinal umber of Monte Carlo simulations: ",x$m)
    }
  } else {
    cat("\nMCMC sample size: ",x$m," with ",x$m/x$iter," chains")
  }
  cat("\nApproximation: ",x$approx)
  cat("\n\nRandom effects: \n")
  print(x$covpars)
  cat("\nFixed effects: \n")
  print(x$fixef)
  if(x$ml){
    cat("\ncAIC: ",round(x$aic,digits))
    cat("\nApproximate R-squared: Conditional: ",round(x$Rsq[1],digits)," Marginal: ",round(x$Rsq[2],digits))
    cat("\nLog-likelihood: ",round(x$logl,digits))
  }
  cat("\n\nModel predictions: \n")
  nT <- length(x$preds)
  for(t in 1:nT){
    cat("\n\U2BC8 Time period ",t,"\n     \U2BA1 Relative risk\n")
    print(summary(exp(rowMeans(x$preds[[t]]$rr))))
    cat("     \U2BA1 Predicted incidence\n")
    print(summary(rowMeans(x$preds[[t]]$y)))
  }
  cat("\n")
}


#' Extracts fixed effect coefficients from a rtsFit object
#' 
#' Extracts the fitted fixed effect coefficients from an `rtsFit` object returned from a call of `rtsFit` or `LA` in the \link[glmmrBase]{Model} class.
#' @param object An `rtsFit` model fit.
#' @param ... Further arguments passed from other methods
#' @return A named vector.
#' @method coef rtsFit
#' @export
coef.rtsFit <- function(object,...){
  pars <- object$coefficients$est[1:object$P]
  names(pars) <- object$coefficients$par[1:object$P]
  return(pars)
}

#' Extracts the log-likelihood from an rtsFit object
#' 
#' Extracts the final log-likelihood value from an rtsFit object. Only returns a value for maximum likelihood 
#' model fits, otherwise it produces an error.
#' @param object An `rtsFit` model fit.
#' @param ... Further arguments passed from other methods
#' @return An object of class `logLik` for maximum likelihood model fits, otherwise it returns an error.
#' @method logLik rtsFit
#' @export
logLik.rtsFit <- function(object, ...){
  if(object$method%in%c("mcmc","vb"))stop("Log likelihood only available for maximum likelihood models.")
  ll <- object$logl
  class(ll) <- "logLik"
  attr(ll,"df") <- object$P + 2 + ifelse(object$nT > 1,1,0)
  attr(ll,"nobs") <- length(object$y)
  attr(ll,"nall") <- length(object$y)
  return(ll)
}

#' Extracts the family from a `rtsFit` object. 
#' 
#' Extracts the \link[stats]{family} from a `rtsFit` object.
#' @param object A `rtsFit` object.
#' @param ... Further arguments passed from other methods
#' @return A \link[stats]{family} object.
#' @method family rtsFit
#' @export
family.rtsFit <- function(object,...){
  return(stats::poisson())
}

#' Extracts the formula from a `rtsFit` object. 
#' 
#' Extracts the \link[stats]{formula} from a `rtsFit` object. Only returns the top level formula. For region models 
#' this is the formula at the region level, otherwise the grid-level formula is returned. No random effects 
#' specifications are included in the returned formula.
#' @param x A `rtsFit` object.
#' @param ... Further arguments passed from other methods
#' @return A \link[stats]{formula} object.
#' @method formula rtsFit
#' @export
formula.rtsFit <- function(x,...){
  fo <- "~ "
  if(length(x$covs)>0){
    for(i in 1:length(covs))fo <- paste0(fo,ifelse(i==1,""," + "),covs[i])
  } else {
    fo <- paste0(fo,"1")
  }
  if(x$region)message("The returned formula does not include the grid level formula for the region model")
  message("The returned formula does not include random effects")
  return(as.formula(fo))
}

#' Extract the Variance-Covariance matrix for a `rtsFit` object
#' 
#' Returns the calculated variance-covariance matrix for a `rtsFit` object that was fit using maximum likelihood methods. 
#' Bayesian models will return an error.
#' @param object A `rtsFit` object.
#' @param ... Further arguments passed from other methods
#' @return A variance-covariance matrix.
#' @method vcov rtsFit
#' @export
vcov.rtsFit <- function(object,...){
  if(object$method%in%c("mcmc","vb"))stop("vcov only available currently for maximum likelihood model fits")
  M <- object$vcov
  rownames(M) <- colnames(M)  <- object$coefficients$par[1:object$P]
  return(object$vcov)
}

#' Predict from a `rtsFit` object
#' 
#' Predictions cannot be generated directly from an `rtsFit` object, rather new predictions should be
#' generated using the original `grid` object. A message is printed to the user.
#' @param object A `rtsFit` object.
#' @param ... Further arguments passed from other methods
#' @return Nothing. Called for effects.
#' @method predict rtsFit
#' @export
predict.rtsFit <- function(object,...){
  message("Predictions cannot be generated directly from a fitted rtsFit object. Predictions are generated by the grid object, see grid$predict().")
}

#' Fitted values from a `rtsFit` object
#' 
#' Fitted values should not be generated directly from an `rtsFit` object, rather fitted values should be
#' generated using the original `grid` object. A message is printed to the user. 
#' @param object A `rtsFit` object.
#' @param ... Further arguments passed from other methods
#' @return Nothing, called for effects.
#' @method fitted rtsFit
#' @export
fitted.rtsFit <- function(object,...){
  message("Fitted values cannot be generated directly from a fitted rtsFit object. See grid$predict() to generate relevant values.")
}

#' Fixed effect confidence intervals for a `rtsFit` object
#' 
#' Returns the computed confidence intervals for a `rtsFit` object.  
#' @param object A `rtsFit` object.
#' @param ... Further arguments passed from other methods
#' @return A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. 
#' @method confint rtsFit
#' @export
confint.rtsFit <- function(object, ...){
  out <- object$coefficients[1:object$P,c("lower","upper")]
  out <- as.matrix(out)
  colnames(out) <- c("2.5%","97.5%")
  return(out)
}

#' Residuals method for a `rtsFit` object
#' 
#' Conditional raw or standardised residuals for `rstFit` objects. The residuals are limited to conditional raw or standardised 
#' residuals currently to avoid copying the often large amount of model data stored in the associated grid object. 
#' @param object A `rtsFit` object.
#' @param type Either "standardized" or "raw"
#' @param ... Further arguments passed from other methods
#' @return A matrix with number of columns corresponding to the number of MCMC samples.
#' @method residuals rtsFit
#' @export
residuals.rtsFit <- function(object, type, ...){
  if(!type%in%c("raw","standardized"))stop("type must be either raw or standardized")
  resids <- matrix(NA,nrow = length(object$y), ncol = ncol(object$y_predicted))
  for(i in 1:ncol(object$y_predicted)) {
    resids[,i] <- object$y - object$y_predicted
    if(type == "standardized") resids[,i] <- resids[,i]/sd(resids[,i])
  }
  return(resids)
}

#' Summarizes a `grid` object
#' 
#' Summarizes `grid` object. 
#' @param object A `grid` object.
#' @param ... Further arguments passed from other methods
#' @return Nothing. Called for effects.
#' @method summary grid
#' @export
summary.grid <- function(object, ...){
  is_setup <- FALSE
  cat("\nAn rts2 grid object\n")
  if(is.null(object$region_data)){
    cat("\n \U2BC8 Boundary\n\n")
    print(head(object$boundary))
  }
  cat("\n \U2BC8 Computational grid\n     \U2BA1 ",nrow(object$grid_data)," cells")
  if(is.null(object$region_data)){
    nT <- sum(grepl("\\bt[0-9]",colnames(object$grid_data)))
    if(nT == 0) {
      nT <- 1
      if("y"%in%colnames(object$grid_data)){
        cat("\n     \U2BA1 Spatial-only case data")
        is_setup <- TRUE
      } else {
        cat("\n     \U2BA1 No case data currently mapped to grid: use points_to_grid() function or manually add a column y or t1, t2, t3, ... containing counts.") 
      }
    } else {
      cat("\n     \U2BA1 Case data for ",nT," time periods")
      is_setup <- TRUE
    }
    cat("\n\n")
    print(head(object$grid_data))
  } else {
    cat("\n \U2BC8 Spatially aggregated count data:\n     \U2BA1 ",nrow(object$region_data)," regions")
    nT <- sum(grepl("\\bt[0-9]",colnames(object$region_data)))
    if(nT == 0){
      nT <- 1
      if("y"%in%colnames(object$region_data)){
        cat("\n     \U2BA1 Spatial-only case data")
        is_setup <- TRUE
      } else {
        cat("\n     \U2BA1 No case data currently in region data: manually add a column y or t1, t2, t3, ... containing counts.") 
      }
    } else {
      cat("\n     \U2BA1 Case data for ",nT," time periods")
      is_setup <- TRUE
    }
    cat("\n     \U2BA1 ",nrow(private$intersection_data)," region-grid intersection areas\n")
    print(head(object$region_data))
  }
  cat("\n \U2BC8 Last model fit \n")
  if(!is.null(private$last_model_fit)){
    print(private$last_model_fit)
  } else {
    cat("     \U2BA1 No model has been fit to these data: see lgcp_ml() and lgcp_bayes()")
  }
  cat("\n")
}

#' Extract predictions from a `grid` object
#'
#' Extract incidence and relative risk predictions. The predictions will be extracted from the last model fit in the `grid` object. 
#' If no previous model fit then use either `grid$lgcp_ml()` or `grid$lgcp_bayes()`, or see `grid$model_fit()` to update the stored model fit.
#'
#' @param object A `grid` object.
#' @param type Vector of character strings. Any combination of "pred", "rr", and "irr", which are,
#' posterior mean incidence (overall and population standardised), relative risk,
#' and incidence rate ratio, respectively.
#' @param irr.lag integer. If "irr" is requested as `type` then the number of time
#' periods lag previous the ratio is in comparison to
#' @param t.lag integer. Extract predictions for previous time periods.
#' @param popdens character string. Name of the column in `grid_data` with the
#' population density data
#' @param verbose Logical indicating whether to print messages to the console
#' @param ... Further arguments passed from other methods
#' @return An `sf` object in which the predictions are stored.
#' @details 
#' Three outputs can be extracted from the model fit:
#'
#' Predicted incidence: If type includes `pred` then `pred_mean_total` and
#' `pred_mean_total_sd` provide the
#' predicted mean total incidence and its standard deviation, respectively.
#' `pred_mean_pp` and `pred_mean_pp_sd` provide the predicted population
#' standardised incidence and its standard deviation. These are added to the grid data or to the 
#' regional data for spatially-aggregated data.
#'
#' Relative risk: if type includes `rr` then the relative risk is reported in
#' the columns `rr` and `rr_sd`. The relative risk here is the exponential
#' of the latent field, which describes the relative difference between
#' expected mean and predicted mean incidence. These are added to the grid data.
#'
#' Incidence risk ratio: if type includes `irr` then the incidence rate ratio (IRR)
#' is reported in the columns `irr` and `irr_sd`. This is the ratio of the predicted
#' incidence in the last period (minus `t_lag`) to the predicted incidence in the
#' last period minus `irr_lag` (minus `t_lag`). For example, if the time period
#' is in days then setting `irr_lag` to 7 and leaving `t_lag=0` then the IRR
#' is the relative change in incidence in the present period compared to a week
#' prior. These are added to the grid data or to the 
#' regional data for spatially-aggregated data.
#' @examples
#' # See examples for grid$lgcp_bayes() and grid$lgcp_ml()
#' @importFrom stats sd
#' @method predict grid
predict.grid <- function(object, 
                         type=c("pred","rr","irr"),
                         irr.lag=NULL,
                         t.lag=0,
                         popdens=NULL,
                         verbose = TRUE,
                         ...){
  object$extract_preds(type,irr.lag,t.lag,popdens,verbose)
  if(type %in% c("pred") & !is.null(object$region_data)){
    return(invisible(object$region_data))
  } else {
    return(invisible(object$grid_data))
  }
}

#' Residuals method for a `grid` object
#' 
#' Conditional raw or standardised residuals are returned for a stored `rtsFit` objects. If no prior model fit is stored,
#' then an error is returned. 
#' @param object A `grid` object.
#' @param type Either "standardized" or "raw"
#' @param ... Further arguments passed from other methods
#' @return A matrix with number of columns corresponding to the number of MCMC samples.
#' @method residuals grid
#' @export
residuals.grid <- function(object, type, ...){
  return(residuals(object$model_fit(),type))
}

#' Calculate Variance-Covariance matrix for a maximum likelihood object stored in `grid`
#' 
#' Returns the variance-covariance matrix for a LGCP object fit using maximum likelihood. If no relevant 
#' model is stored then the function returns an error
#' @param object A `grid` object.
#' @param ... Further arguments passed from other methods
#' @return A variance-covariance matrix.
#' @method vcov grid
#' @export
vcov.grid <- function(object,...){
  return(vcov(object$model_fit()))
}

#' Extracts the family from a `grid` object. 
#' 
#' Extracts the \link[stats]{family} from a `grid` object.
#' @param object A `grid` object.
#' @param ... Further arguments passed from other methods
#' @return A \link[stats]{family} object.
#' @method family rtsFit
#' @export
family.grid <- function(object,...){
  return(stats::poisson())
}

#' Extracts the formula from a `grid` object. 
#' 
#' Extracts the \link[stats]{formula} from a `rtsFit` object stored in a grid object. Only returns the top level formula. For region models 
#' this is the formula at the region level, otherwise the grid-level formula is returned. No random effects 
#' specifications are included in the returned formula.
#' @param x A `grid` object.
#' @param ... Further arguments passed from other methods
#' @return A \link[stats]{formula} object.
#' @method formula grid
#' @export
formula.grid <- function(x,...){
  return(formula(x$model_fit()))
}