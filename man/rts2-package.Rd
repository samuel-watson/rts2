\name{rts2-package}
\alias{rts2-package}
\alias{rts2}
\docType{package}
\title{
\packageTitle{rts2}
}
\description{
\packageDescription{rts2}
\code{rts2} provides several estimators for the Log Gaussian Cox Process Model (LGCP). The LGCP is a stochastic Poisson model used for modelling case counts of phenomena of interest, and is particularly useful for predicting risk across an area of interest, such as in disease surveillance applications.
}
\section{Workflow}{
Most of the functionality of the \pkg{rts2} package is provided by the \code{grid} class. The computational strategy for the LGCP is to divide up the area of interest into a regular grid and aggregate case counts within cells. For models with count data aggregated to an irregular set of polygons, such as census tracts, the latent surface is also modelled as a regular grid. A typical workflow using this package would be:
\enumerate{
\item Create a new \code{grid} object, e.g. \code{g1 <- grid$new(poly, cellsize = 0.1)}. The class is initialized with either a single polygon describing the area of interest or a collection
of polygons if spatially aggregated data are used. The \pkg{sf} package is used for all spatial data.
\item If the location (and times) of cases are available (i.e. the data are not spatially aggregated), then we map the points to the computational grid. The function \link[rts2]{create_points} can generate point data in the correct \pkg{sf} format. The member function \code{points_to_grid} will then map these data to the grid. Counts can also be manually added to grid data. For region data, since the counts are assumed to be already aggregated, these must be manually provided by the user. The case counts must appear in columns with specific names. If there is only a single time period then the counts
must be in a column named \code{y}. If there are multiple time periods then the counts must be in columns names \code{t1}, \code{t2}, \code{t3},... Associated columns labelled \code{date1}, \code{date2}, etc. will permit use of some functionality regarding specific time intervals.
\item If any covariates are to be used for the modelling, then these can be mapped to the compuational grid using the function \code{add_covariates()}. Other functions, \code{add_time_indicators()} and \code{get_dow()} will also generate relevant temporal indicators where required. At a minimum we would recommend including
a measure of population density.
\item Fit a model. There are multiple methods for model fitting, which are available through the member functions \code{lgcp_ml()} and \code{lgcp_bayes()} for maximum likelihood and Bayesian approaches, respectively. The results are stored internally and optionally returned as a \code{rtsFit} object.
\item Summarise the output. The main functions for summarising the output are \code{extract_preds()}, which will generate predictions of relative risk, incidence rate ratios, and predicted incidence, and \code{hotspots()}, which will estimate probabilities that these statistics exceed given thresholds. For spatially-aggregated data models, the relative risk applies to the grid, whereas rate ratios and predicted incidence applies to the areas.
\item Predictions can be visualised or aggregated to relevant geographies with the \code{plot()} and \code{aggregate()} functions.
}
}
\section{Estimation methods and model specification}{
The \pkg{rts2} package provide several methods for estimating the model and drawing samples from the latent surface.
\itemize{
\item Maximum Likelihood. We include stochastic maximum likelihood estimation methods including both Markov Chain Monte Carlo (MCMC) Maximum Likelihood and Stochastic Approximation Expectation Maximisation (SAEM). MCMC-ML can use Newton-Raphson, quasi-Newton, or derivative free methods to estimate the model parameters. Both algorithms have three steps: 1. Sample the random effects using MCMC; 2. Estimate the fixed effect parameters conditional on the sampled random effects; 3. Estimate the covariance parameters. The process is iterated until convergence. Stochastic maximum likelihood estimators are provided by the function \code{grid$lgcp_ml()}.
\item Bayesian. We also include Bayesian estimation of the model using Stan via either \pkg{rstan} or \pkg{cmdstanr}, and allow both MCMC and Variational Bayes methods. 
}

The LGCP can be computationally complex and scales poorly with sample size (number of grid cells and time periods), due to the large covariance matrix that must be inverted to estimate the covariance parameters. We offer several strategies and approximations for efficient model fitting:
\itemize{
\item Gaussian Process Approximations. The package includes both Hilbert Space Gaussian Process (see Solin and Särkkä (2020) <doi:10.1007/s11222-019-09886-w> and Riutort-Mayol et al (2020) <arXiv:2004.11408>) and the Nearest Neighbour Gaussian Process (Datta et al (2016) <doi:10.1080/01621459.2015.1044091>).
\item For spatio-temporal models we use a "spatial innovation" formulation of the spatio-temporal Gaussian process, for which the computational complexity is linear in the number of time periods.
}
} 
\section{Package development}{
The package is still in development and there may still be bugs and errors. While we do not expect the general user interface to change there may be changes to the underlying library as well as new additions and functionality.
}
\author{
\packageAuthor{rts2}

Maintainer: \packageMaintainer{rts2}
}
\keyword{ package }
\references{
Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
}
