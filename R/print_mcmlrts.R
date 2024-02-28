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
#' @param object an object of class "`rtsFitSummary`" 
#' @param ... Further arguments passed from other methods
#' @details 
#' `print.rtsFitSummary` prints the summary of an rtsFit, see \link[rts2]{summary.rtsFit}
#' @return No return value, called for side effects.
#' @method print rtsFitSummary
#' @export
print.rtsFitSummary <- function(object,...){
  cat("\nAn rts model fit summary\n")
  digits <- 4
  cat(ifelse(!object$ml,
             "Bayesian Log Gaussian Cox Process Model\n",
             "Maximum Likelihood Log Gaussian Cox Process Model\n"))
  cat("\nUsing: ",object$algo)
  if(object$ml){
    if(grepl("MCMC Maximum Likelihood Expectation Maximisation",object$algo)){
      cat("\nNumber of Monte Carlo simulations per iteration: ",object$m," with tolerance ",object$tol,"\n\n")
    } else {
      cat("\nFinal umber of Monte Carlo simulations: ",object$m)
    }
  } else {
    cat("\nMCMC sample size: ",object$m," with ",object$m/object$iter," chains")
  }
  cat("\nApproximation: ",object$approx)
  cat("\n\nRandom effects: \n")
  print(object$covpars)
  cat("\nFixed effects: \n")
  print(object$fixef)
  if(object$ml){
    cat("\ncAIC: ",round(object$aic,digits))
    cat("\nApproximate R-squared: Conditional: ",round(object$Rsq[1],digits)," Marginal: ",round(object$Rsq[2],digits))
    cat("\nLog-likelihood: ",round(object$logl,digits))
  }
  cat("\n\nModel predictions: \n")
  nT <- length(object$preds)
  for(t in 1:nT){
    cat("\n\U2BC8 Time period ",t,"\n     \U2BA1 Relative risk\n")
    print(summary(exp(rowMeans(object$preds[[t]]$rr))))
    cat("     \U2BA1 Predicted incidence\n")
    print(summary(rowMeans(object$preds[[t]]$y)))
  }
}


#' Extracts coefficients from a mcml object
#' 
#' Extracts the coefficients from an `mcml` object returned from a call of `MCML` or `LA` in the \link[glmmrBase]{Model} class.
#' @param object An `mcml` model fit.
#' @param ... Further arguments passed from other methods
#' @return A data frame summarising the parameters including the random effects.
#' @method coef rtsFit
#' @export
coef.rtsFit <- function(object,...){
  return(object$coefficients)
}

#' Extracts the log-likelihood from an mcml object
#' 
#' Extracts the final log-likelihood value from an mcml object returned from call of `MCML` or `LA` in the \link[glmmrBase]{Model} class. The fitting algorithm estimates
#' the fixed effects, random effects, and covariance parameters all separately. The log-likelihood is separable in the fixed and covariance parameters, so one can return 
#' the log-likelihood for either component, or the overall log-likelihood.
#' @param object An `mcml` model fit.
#' @param fixed Logical whether to include the log-likelihood value from the fixed effects.
#' @param covaraince Logical whether to include the log-likelihood value from the covariance parameters.
#' @param ... Further arguments passed from other methods
#' @return A numeric value. If both `fixed` and `covariance` are FALSE then it returns NA.
#' @method logLik rtsFit
#' @export
logLik.rtsFit <- function(object, fixed = TRUE, covariance = TRUE, ...){
  if(object$method%in%c("mcmc","vb"))stop("Not a maximum likelihood model fit.")
  ll <- 0
  if(fixed) ll <- ll + object$logl
  if(covariance) ll <- ll + object$logl_theta
  if(!fixed & !covariance) ll <- NA
  return(NA)
}

