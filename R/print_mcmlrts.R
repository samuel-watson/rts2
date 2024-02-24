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
  cat("\nUsing: ",algo)
  if(ml){
    cat("\nNumber of Monte Carlo simulations per iteration: ",x$m," with tolerance ",x$tol,"\n\n")
  } else {
    cat("\nMCMC sample size: ",x$m," with ",x$m/x$iter," chains")
  }
  approx <- switch(x$approx,
                   "none" = "None",
                   "hsgp" = "Hilbert Space Gaussian Process",
                   "nngp" = "Nearest Neighbour Gaussian Process")
  cat("\nApproximation: ",approx)
  if(ml){
    colrange <- 2:7
  } else {
    colrange <- c(2,3,6,7)
  }
  pars <- x$coefficients[,colrange]
  if(ml){
    colnames(pars) <- c("Estimate","Std. Err.","z value","p value","2.5% CI","97.5% CI")
  } else {
    colnames(pars) <- c("Posterior\nmean","Posterior\nStd. Dev.","2.5% CrI","97.5% CrI")
  }
  rownames(pars) <- x$coefficients$par
  pars <- apply(pars,2,round,digits = digits)
  cat("\nRandom effects: \n")
  print(pars[(x$P+1):nrow(pars),1])
  cat("\nFixed effects: \n")
  print(pars[1:x$P,])
  if(ml){
    cat("\ncAIC: ",round(x$aic,digits))
    cat("\nApproximate R-squared: Conditional: ",round(x$Rsq[1],digits)," Marginal: ",round(x$Rsq[2],digits))
    cat("\nLog-likelihood: ",round(x$logl,digits))
  }
  return(invisible(pars))
}

#' Summarises an mcmlrts fit output
#' 
#' Summary method for class "`rtsFit`"
#' 
#' @param object an object of class "`rtsFit`" as a result of a call to fit_ml
#' @param ... Further arguments passed from other methods
#' @details 
#' `print.rtsFit` tries to replicate the output of other regression functions, such
#' as `lm` and `lmer` reporting parameters, standard errors, and z- and p- statistics.
#' The z- and p- statistics should be interpreted cautiously however, as generalised
#' linear miobjected models can suffer from severe small sample biases where the effective
#' sample size relates more to the higher levels of clustering than individual observations.
#' @return A list with random effect names and a data frame with random effect mean and credible intervals
#' @method summary mcmlrts
#' @export
summary.rtsFit <- function(object,...){
  digits <- 2
  pars <- print(object)
  ## summarise random effects
  dfre <- data.frame(Mean = round(apply(object$re.samps,2,mean),digits = digits), 
                     lower = round(apply(object$re.samps,2,function(i)stats::quantile(i,0.025)),digits = digits),
                     upper = round(apply(object$re.samps,2,function(i)stats::quantile(i,0.975)),digits = digits))
  colnames(dfre) <- c("Estimate","2.5% CI","97.5% CI")
  cat("\nRandom effects estimates\n")
  print(dfre)
  ## add in model fit statistics
  return(invisible(list(coefficients = pars,re.terms = dfre)))
}


#' Extracts coefficients from a mcml object
#' 
#' Extracts the coefficients from an `mcml` object returned from a call of `MCML` or `LA` in the \link[glmmrBase]{Model} class.
#' @param object An `mcml` model fit.
#' @param ... Further arguments passed from other methods
#' @return A data frame summarising the parameters including the random effects.
#' @method coef mcml
#' @export
coef.rtsFit <- function(object,...){
  return(object$coefficients)
}

#' Extracts the computed AIC from an mcml object
#' 
#' Extracts the conditional Akaike Information Criterion from an mcml object returned from call of `MCML` or `LA` in the \link[glmmrBase]{Model} class.
#' @param object An `mcml` model fit.
#' @param ... Further arguments passed from other methods
#' @return A numeric value.
#' @method extractAIC mcml
#' @export
extractAIC.rtsFit <- function(object,...){
  if(object$method%in%c("mcmc","vb"))stop("Not a maximum likelihood model fit.")
  return(object$aic)
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
#' @method logLik mcml
#' @export
logLik.rtsFit <- function(object, fixed = TRUE, covariance = TRUE, ...){
  if(object$method%in%c("mcmc","vb"))stop("Not a maximum likelihood model fit.")
  ll <- 0
  if(fixed) ll <- ll + object$logl
  if(covariance) ll <- ll + object$logl_theta
  if(!fixed & !covariance) ll <- NA
  return(NA)
}

#' Extracts the fixed effect estimates
#' 
#' Extracts the fixed effect estimates from an mcml object returned from call of `MCML` or `LA` in the \link[glmmrBase]{Model} class.
#' @param object An `mcml` model fit.
#' @param ... Further arguments passed from other methods
#' @return A named, numeric vector of fixed-effects estimates.
#' @method ranef mcml
#' @export
fixef.rtsFit <- function(object,...){
  fixed <- object$coefficients$est[1:object$P]
  names(fixed) <- object$coefficients$par[1:object$P]
  return(fixed)
}

#' Extracts the random effect estimates
#' 
#' Extracts the random effect estimates or samples from an mcml object returned from call of `MCML` or `LA` in the \link[glmmrBase]{Model} class.
#' @param object An `mcml` model fit.
#' @param ... Further arguments passed from other methods
#' @return A matrix of dimension (number of fixed effects ) x (number of MCMC samples). For Laplace approximation, the number of "samples" equals one.
#' @method ranef mcml
#' @export
ranef.rtsFit <- function(object,...){
  return(object$re.samps)
}