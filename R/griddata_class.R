#' R6 class holding sf grid data with data and analysis functions
#'
#' Grid data consists of the computational grid over the area of interest. Outcomes and
#' covariates are projected onto the grid, which can then be sent to the LGCP model.
#' @importFrom R6 R6Class
#' @importFrom stats weighted.mean sd model.matrix
#' @importFrom utils stack capture.output
#' @export
grid <- R6::R6Class("grid",
                         public = list(
                           #' @field grid_data sf object specifying the computational grid for the analysis
                           grid_data = NULL,
                           #' @field region_data sf object specifying an irregular lattice, such as census areas, 
                           #' within which case counts are aggregated. Only used if polygon data are provided on
                           #' class initialisation.
                           region_data = NULL,
                           #' @field priors list of prior distributions for the analysis
                           priors = NULL,
                           #' @field boundary sf object showing the boundary of the area of interest
                           boundary = NULL,
                           #' @description
                           #' Create a new griddata object
                           #'
                           #' Produces a regular grid over an area of interest as an sf object
                           #'
                           #' Given a contiguous boundary describing an area of interest, which is stored as an sf
                           #' object of a regular grid within the limits of the boundary at `$grid_data`. The boundary
                           #' is also stored in the object as `$boundary`
                           #'
                           #' @param poly An sf object containing either one polygon describing the area of interest or multiple polygons
                           #' representing survey or census regions in which the case data counts are aggregated
                           #' @param cellsize The dimension of the grid cells
                           #' @param verbose Logical indicating whether to provide feedback to the console.
                           #' @return NULL
                           #' @examples
                           #' b1 = sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
                           #' g1 <- grid$new(b1,0.5)
                           initialize = function(poly,
                                                 cellsize,
                                                 verbose = TRUE){

                             if(!is(poly,"sf"))stop("data not sf")
                             if(!is(cellsize,"numeric"))stop("cellsize not numeric")
                             if(nrow(poly)>1 & verbose)message("Multiple polygons in data. Assuming analysis uses counts aggregated to an irregular lattice and not point data.")
                             #if(nrow(boundary)!=1)stop("boundary should only contain one polygon")
                             
                             if(nrow(poly)==1){
                               bgrid <- sf::st_make_grid(poly,cellsize = cellsize)
                               bgrid <- sf::st_sf(bgrid)
                               idx1<- sf::st_contains(y=bgrid,x=poly)
                               idx2<- sf::st_intersects(y=bgrid,x=poly)
                               idx <- c(unlist( sapply( idx1, `[`) ),unlist( sapply( idx2, `[`) ))
                               bgrid <- bgrid[idx,]
                               
                               #class(bgrid) <- c(class(bgrid),"rts_grid")
                               self$boundary <- poly
                               bgrid <- bgrid[!duplicated(sf::st_coordinates(sf::st_centroid(bgrid))),]
                               self$grid_data <- bgrid
                             } else if(nrow(poly)>1){
                               boundary <- sf::st_union(poly)
                               boundary <- sf::st_sfc(boundary)
                               self$boundary <- sf::st_sf(boundary)
                               bboundary<<-boundary
                               bgrid <- sf::st_make_grid(self$boundary,cellsize = cellsize)
                               bgrid <- sf::st_sf(bgrid)
                               idx1<- sf::st_contains(y=bgrid,x=self$boundary)
                               idx2<- sf::st_intersects(y=bgrid,x=self$boundary)
                               idx <- c(unlist( sapply( idx1, `[`) ),unlist( sapply( idx2, `[`) ))
                               bgrid <- bgrid[idx,]
                               bgrid <- bgrid[!duplicated(sf::st_coordinates(sf::st_centroid(bgrid))),]
                               self$grid_data <- bgrid
                               
                               sf::st_agr(poly) = "constant"
                               sf::st_agr(self$grid_data) = "constant"
                               self$grid_data$grid_id <- 1:nrow(self$grid_data)
                               poly$region_id <- 1:nrow(poly)
                               tmp <- suppressWarnings(sf::st_intersection(self$grid_data[,"grid_id"],poly[,"region_id"]))
                               n_Q <- nrow(tmp)
                               rID <- poly$region_id
                               tmp$area <- as.numeric(sf::st_area(tmp))
                               a1 <-rep(aggregate(tmp$area,list(tmp$region_id),sum)$x,unname(table(tmp$region_id)))
                               tmp$w <- tmp$area/a1
                               
                               self$region_data <- poly
                               private$intersection_data <- tmp
                             }
                             if(verbose){
                               msg <- paste0("Created grid with ",nrow(self$grid_data)," cells.")
                               if(!is.null(self$region_data)){
                                 msg <- paste0(msg," Region data has ",nrow(self$region_data)," areas.",
                                               "There are ",nrow(private$intersection_data)," intersection areas.")
                               }
                               message(msg)
                             }
                           },
                           #' @description
                           #' Prints the $grid_data sf object
                           print = function(){
                             print(self$grid_data)
                           },
                           #' @description
                           #' Plots the grid data
                           #'
                           #' @details
                           #' If `zcol` is not specified then only the geometry is plotted, otherwise the covariates specified will be plotted.
                           #' The user can also use sf plotting functions on grid$grid_data directly.
                           #' @param zcol Vector of strings specifying names of columns of `grid_data` to plot
                           #' @return A plot
                           #' @examples
                           #' b1 = sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
                           #' g1 <- grid$new(b1,0.5)
                           #' g1$plot()
                           plot = function(zcol){
                             if(missing(zcol)){
                               plot(sf::st_geometry(self$grid_data))
                             } else {
                               plot(self$grid_data[,zcol])
                             }
                           },
                           #' @description
                           #' Generates case counts of points over the grid
                           #'
                           #' Counts the number of cases in each time period in each grid cell
                           #'
                           #' Given the sf object with the point locations and date output from
                           #' `create_points()`, the functions will add columns to `grid_data` indicating
                           #' the case count in each cell in each time period.
                           #' @details
                           #' Case counts are generated for each grid cell for each time period. The user
                           #' can specify the length of each time period; currently `day`, `week`, and `month`
                           #' are supported.
                           #'
                           #' The user must also specify the number of time periods to include with the
                           #' `laglength` argument. The total number of time periods is the specified lag
                           #' length counting back from the most recent case. The columns in the output
                           #' will be named `t1`, `t2`,... up to the lag length, where the highest number
                           #' is the most recent period.
                           #' @param point_data sf object describing the point location of cases with a column
                           #' `t` of the date of the case in YYYY-MM-DD format. See \link[rts2]{create_points}
                           #' @param t_win character string. One of "day", "week", or "month" indicating the
                           #' length of the time windows in which to count cases
                           #' @param laglength integer The number of time periods to include counting back from the most
                           #' recent time period
                           #' @param verbose Logical indicating whether to report detailed output
                           #' @return NULL
                           #' @seealso \link[rts2]{create_points}
                           #' @examples
                           #' b1 <- sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
                           #' g1 <- grid$new(b1,0.5)
                           #' dp <- data.frame(y=runif(10,0,3),x=runif(10,0,3),date=paste0("2021-01-",11:20))
                           #' dp <- create_points(dp,pos_vars = c('y','x'),t_var='date')
                           #' g1$points_to_grid(dp, laglength=5)
                           points_to_grid = function(point_data,
                                                     t_win = c("day"),
                                                     laglength = 14,
                                                     verbose = TRUE){

                             if(!is(point_data,"sf"))stop("points not sf")


                             if(sf::st_crs(point_data)!=sf::st_crs(self$grid_data)){
                               warning("CRS not equal. Setting st_crs(point_data)==st_crs(self$grid_data)")
                               sf::st_crs(point_data) <- sf::st_crs(self$grid_data)
                             }

                             if("t"%in%colnames(point_data)){

                               if(!t_win%in%c("day","week","month"))stop("t_win not day, week, or month")

                               #get unique time values to summarise over
                               tvals <- c(as.Date(min(point_data$t)),as.Date(max(point_data$t)))
                               yvals <- lubridate::year(tvals[1]):lubridate::year(tvals[2])
                               tuniq <- tvals[1]:tvals[2]
                               tuniq <- as.Date(tuniq,origin=as.Date("1970-01-01"))

                               if(t_win=="day"){
                                 tdat <- paste0(lubridate::yday(point_data$t),".",lubridate::year(point_data$t))
                                 tuniq <- paste0(lubridate::yday(tuniq),".",lubridate::year(tuniq))
                               }
                               if(t_win=="week"){
                                 tdat <- paste0(lubridate::week(point_data$t),".",lubridate::year(point_data$t))
                                 tuniq <- paste0(lubridate::week(tuniq),".",lubridate::year(tuniq))
                               }
                               if(t_win=="month"){
                                 tdat <- paste0(lubridate::month(point_data$t),".",lubridate::year(point_data$t))
                                 tuniq <- paste0(lubridate::month(tuniq),".",lubridate::year(tuniq))
                               }

                               tuniq <- tuniq[(length(tuniq)-laglength+1):length(tuniq)]



                               for(i in 1:length(tuniq))
                               {
                                 self$grid_data$y <-  lengths(sf::st_intersects(self$grid_data,
                                                                           point_data[tdat==tuniq[i],]))
                                 if(length(tuniq)>1)colnames(self$grid_data)[length(colnames(self$grid_data))] <- paste0("t",i)
                                 self$grid_data$d <- min(point_data[tdat==tuniq[i],]$t)
                                 colnames(self$grid_data)[length(colnames(self$grid_data))] <- paste0("date",i)
                               }
                               
                             } else {
                               self$grid_data$y <-  lengths(sf::st_intersects(self$grid_data,
                                                                         point_data))
                             }

                             if(verbose)message("added points data to grid data")
                           },
                           #' @description
                           #' Adds covariate data to the grid
                           #'
                           #' Maps spatial, temporal, or spatio-temporal covariate data onto the grid
                           #'
                           #' @details
                           #' *ADDING COVARIATES*
                           #' *Spatially-varying data only* 
                           #' `cov_data` is an sf object describing covariate
                           #' values for a set of polygons over the area of interest. The values are mapped
                           #' onto `grid_data`. For each grid cell in `grid_data` a weighted
                           #' average of each covariate listed in `zcols` is generated with weights either
                           #' equal to the area of intersection of the grid cell and the polygons in
                           #' `cov_data` (`weight_type="area"`), or this area multiplied by the population
                           #' density of the polygon for population weighted (`weight_type="pop"`). Columns
                           #' with the names in `zcols` are added to the output.
                           #'
                           #' *Temporally-varying only data* 
                           #' `cov_data` is a data frame with number of rows
                           #' equal to the number of time periods. One of the columns must be called `t` and
                           #' have values from 1 to the number of time periods. The other columns of the data
                           #' frame have the values of the covariates for each time period. See
                           #' `get_dow()` for day of week data. A total of
                           #' length(zcols)*(number of time periods) columns are added to the output: for each
                           #' covariate there will be columns appended with each time period number. For example,
                           #' `dayMon1`, `dayMon2`, etc.
                           #'
                           #' *Spatially and temporally varying data* 
                           #' There are two ways to add data that
                           #' vary both spatially and temporally. The final output for use in analysis must
                           #' have a column for each covariate and each time period with the same name appended
                           #' by the time period number, e.g. `covariateA1`,`covariateA2`,... If the covariate
                           #' values for different time periods are in separate sf objects, one can follow
                           #' the method for spatially-varying only data above and append the time period number
                           #' using the argument `t_label`. If the values for different time periods are in the same
                           #' sf object then they should be named as described above and then can be added
                           #' as for spatially-varying covariates, e.g. `zcols=c("covariateA1","covariateA2")`.
                           #'
                           #' @param cov_data sf object or data.frame. See details.
                           #' @param zcols vector of character strings with the names of the columns of `cov_data`
                           #' to include
                           #' @param weight_type character string. Either "area" for area-weighted average or "pop"
                           #' for population-weighted average
                           #' @param popdens character string. The name of the column in `cov_data` with the
                           #' population density. Required if weight_type="pop"
                           #' @param verbose logical. Whether to provide a progress bar
                           #' @param t_label integer. If adding spatio-temporally varying data by time period,
                           #' this time label should be appended to the column name. See details.
                           #' @return NULL
                           #' @examples
                           #' b1 <-  sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
                           #' g1 <- grid$new(b1,0.5)
                           #' cov1 <- grid$new(b1,0.8)
                           #' cov1$grid_data$cov <- runif(nrow(cov1$grid_data))
                           #' g1$add_covariates(cov1$grid_data,
                           #'                   zcols="cov",
                           #'                   verbose = FALSE)
                           add_covariates = function(cov_data,
                                                     zcols,
                                                     weight_type="area",
                                                     popdens=NULL,
                                                     verbose=TRUE,
                                                     t_label=NULL){

                             if(!weight_type%in%c("area","pop"))stop("Type must be area or pop.")
                             if((weight_type=="pop"&is.null(popdens)))stop("Pop. dens. variable not found.")
                             if(weight_type=="pop"&!is.null(popdens)){
                               if(!popdens%in%colnames(cov_data))stop("Pop. dens. variable not found.")
                             }
                             if(any(!zcols%in%colnames(cov_data)))stop("variable names not in cov_data")
                             if(!is(cov_data,"sf"))if(!"t"%in%colnames(cov_data))stop("not column named t in cov_data")

                             if(is(cov_data,"sf")){
                               sf::st_agr(cov_data) = "constant"
                               sf::st_agr(self$grid_data) = "constant"
                               self$grid_data$grid_id <- 1:nrow(self$grid_data)
                               cov_data$region_id <- 1:nrow(cov_data)
                               if(weight_type=="pop"){
                                 tmp <- suppressWarnings(sf::st_intersection(g1$grid_data[,"grid_id"],cov_data[,c("region_id",popdens)]))
                               } else {
                                 tmp <- suppressWarnings(sf::st_intersection(self$grid_data[,"grid_id"],cov_data[,"region_id"]))
                               }
                               n_Q <- nrow(tmp)
                               tmp$area <- as.numeric(sf::st_area(tmp))
                               tmp <- tmp[order(tmp$grid_id),]
                               a1 <-rep(aggregate(tmp$area,list(tmp$grid_id),sum)$x,unname(table(tmp$grid_id)))
                               tmp$w <- tmp$area/a1
                               if(weight_type=="pop"){
                                 tmp$w <- tmp$w*tmp[,popdens,drop=TRUE]
                                 a1 <-rep(aggregate(tmp$w,list(tmp$grid_id),sum)$x,unname(table(tmp$grid_id)))
                                 tmp$w <- tmp$w/a1
                               }
                               
                               vals <- matrix(nrow=nrow(self$grid_data),ncol=length(zcols))

                               if(verbose)cat("Overlaying geographies\n")

                               for(i in 1:nrow(self$grid_data)){
                                 for(j in 1:length(zcols)){
                                   vals[i,j] <- sum(cov_data[tmp[tmp$grid_id==i,]$region_id,zcols[j],drop=TRUE]*tmp[tmp$grid_id==i,]$w)
                                 }

                                 if(verbose)cat("\r",progress_bar(i,nrow(self$grid_data)))
                               }

                               for(j in 1:length(zcols)){
                                 self$grid_data$x <- vals[,j]
                                 if(is.null(t_label)){
                                   colnames(self$grid_data)[length(colnames(self$grid_data))] <- zcols[j]
                                 } else {
                                   colnames(self$grid_data)[length(colnames(self$grid_data))] <- paste0(zcols[j],t_label)
                                 }

                               }
                             } else {
                               nT <- max(cov_data$t)
                               for(j in zcols){
                                 for(t in 1:nT){
                                   self$grid_data$x <- cov_data[cov_data$t==t,j]
                                   colnames(self$grid_data)[length(colnames(self$grid_data))] <- paste0(j,t)
                                 }
                               }
                             }

                             if(verbose)message(paste0("added covariates ",zcols))

                           },
                           #' @description
                           #' Generate day of week data
                           #'
                           #' Create data frame with day of week indicators
                           #'
                           #' Generates a data frame with indicator
                           #' variables for each day of the week for use in the `add_covariates()` function.
                           #'@return data.frame with columns `t`, `day`, and `dayMon` to `daySun`
                           #'@examples
                           #' b1 <- sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
                           #' g1 <- grid$new(b1,0.5)
                           #' dp <- data.frame(y=runif(10,0,3),x=runif(10,0,3),date=paste0("2021-01-",11:20))
                           #' dp <- create_points(dp,pos_vars = c('y','x'),t_var='date')
                           #' g1$points_to_grid(dp, laglength=5)
                           #' g1$get_dow()
                           get_dow = function(){
                             nT <- length(colnames(self$grid_data)[grepl("t[0-9]",colnames(self$grid_data))])
                             dw <- data.frame(t=1:nT,day=NA)
                             for(i in 1:nT){
                               dw$day[i] <- as.character(lubridate::wday(as.data.frame(self$grid_data)[1,paste0("date",i)],
                                                                         label = TRUE))
                             }
                             dx <- model.matrix(~day-1,data=dw)
                             dw <- cbind(dw,as.data.frame(dx))

                             return(dw)
                           },
                           #' @description
                           #' Fit an (approximate) log-Gaussian Cox Process model
                           #'
                           #' @details
                           #' *BAYESIAN MODEL FITTING*
                           #' The grid data must contain columns `t*`, giving the case
                           #' count in each time period (see `points_to_grid`), as well as any covariates to include in the model
                           #' (see `add_covariates`) and the population density.
                           #'
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
                           #' The argument `approx` specifies whether to use a full LGCP model (`approx='none'`) or whether
                           #' to use either a nearest neighbour approximation (`approx='nngp'`) or a "Hilbert space" approximation
                           #' (`approx='hsgp'`). For full details of NNGPs see XX and for Hilbert space approximations see references (1) and (2).
                           #'                           #'
                           #' *Priors*
                           #' For Bayesian model fitting, the priors should be provided as a list to the griddata object:
                           #'```
                           #'griddata$priors <- list(
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
                           #' @param popdens character string. Name of the population density column
                           #' @param covs vector of character string. Base names of the covariates to
                           #' include. For temporally-varying covariates only the stem is required and not
                           #' the individual column names for each time period (e.g. `dayMon` and not `dayMon1`,
                           #' `dayMon2`, etc.)
                           #' @param approx Either "rank" for reduced rank approximation, or "nngp" for nearest 
                           #' neighbour Gaussian process. 
                           #' @param m integer. Number of basis functions for reduced rank approximation, or
                           #' number of nearest neighbours for nearest neighbour Gaussian process. See Details.
                           #' @param L integer. For reduced rank approximation, boundary condition as proportionate extension of area, e.g.
                           #' `L=2` is a doubling of the analysis area. See Details.
                           #' @param model Either "exp" for exponential covariance function or "sqexp" for squared exponential
                           #' covariance function
                           #' @param known_theta An optional vector of two values of the covariance parameters. If these are provided
                           #' then the covariance parameters are assumed to be known and will not be estimated.
                           #' @param dir character string. Directory to save ouptut.
                           #' @param iter_warmup integer. Number of warmup iterations
                           #' @param iter_sampling integer. Number of sampling iterations
                           #' @param chains integer. Number of chains
                           #' @param parallel_chains integer. Number of parallel chains
                           #' @param priors list. See Details
                           #' @param verbose logical. Provide feedback on progress
                           #' @param vb Logical indicating whether to use variational Bayes (TRUE) or full MCMC sampling (FALSE)
                           #' @param use_cmdstanr logical. Defaults to false. If true then cmdstanr will be used
                           #' instead of rstan.
                           #' @param ... additional options to pass to `$sample()``, see \link[cmdstanr]{sample}
                           #' @return A \link[rstan]{stanfit} or a \link[cmdstanr]{CmdStanMCMC} object
                           #' @seealso points_to_grid, add_covariates
                           #' @examples
                           #' \dontrun{
                           #' b1 <- sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
                           #' g1 <- grid$new(b1,0.5)
                           #' dp <- data.frame(y=runif(10,0,3),x=runif(10,0,3),date=paste0("2021-01-",11:20))
                           #' dp <- create_points(dp,pos_vars = c('y','x'),t_var='date')
                           #' cov1 <- grid$new(b1,0.8)
                           #' cov1$grid_data$cov <- runif(nrow(cov1$grid_data))
                           #' g1$add_covariates(cov1,
                           #'                   zcols="cov",
                           #'                   verbose = FALSE)
                           #' g1$points_to_grid(dp, laglength=5)
                           #' g1$priors <- list(
                           #'   prior_lscale=c(0,0.5),
                           #'   prior_var=c(0,0.5),
                           #'   prior_linpred_mean=c(0),
                           #'   prior_linpred_sd=c(5)
                           #'   )
                           #' res <- g1$lgcp_fit(popdens="cov")
                           #' }
                           lgcp_fit = function(popdens,
                                               covs=NULL,
                                               approx = "nngp",
                                               m=10,
                                               L=1.5,
                                               model = "exp",
                                               known_theta = NULL,
                                               dir=NULL,
                                               iter_warmup=500,
                                               iter_sampling=500,
                                               chains=3,
                                               parallel_chains=3,
                                               verbose=TRUE,
                                               vb = FALSE,
                                               use_cmdstanr = FALSE,
                                               ...){

                             if(verbose)if(is.null(dir))message("dir not set, files will be lost after session restart")
                             if(m<1)stop("m must be positive")
                             if(!approx%in%c('nngp','hsgp','none'))stop("approx must be one of nngp, hsgp, or none")
                             if(m >25 & verbose)warning("m is large, sampling may take a long time.")
                             if(!is.null(self$region_data)&verbose)message("Using regional data model.")
                             
                             #prepare data for model fit
                             datlist <- private$prepare_data(m,model,approx,popdens,covs,verbose,TRUE)
                             if(!is.null(known_theta)){
                               if(length(known_theta)!=2)stop("Theta should be of length 2")
                               datlist$known_cov <- 1
                               datlist$sigma_data <- as.array(known_theta[1])
                               if(approx=="hsgp"){
                                 datlist$phi_data <- as.array(c(known_theta[2],known_theta[2]))
                               } else {
                                 datlist$phi_data <- as.array(known_theta[2])
                               }
                             } else {
                               datlist$known_cov <- 0
                               datlist$sigma_data <- c()
                               datlist$phi_data <- c()
                             }

                             if(approx == "hsgp"){
                               if(!is.null(self$region_data)){
                                 filen <- "approxlgcp_region_cmd.stan"
                                 fname <- "approxlgcp_region"
                               } else {
                                 filen <- "approxlgcp_cmd.stan"
                                 fname <- "approxlgcp"
                               }
                             } else if(approx == "nngp"){
                               if(!is.null(self$region_data)){
                                 filen <- "approxlgcp_nngp_region_cmd.stan"
                                 fname <- "approxlgcp_nngp_region"
                               } else {
                                 filen <- "approxlgcp_nngp_cmd.stan"
                                 fname <- "approxlgcp_nngp"
                               }
                             } else {
                               if(!is.null(self$region_data)){
                                 filen <- "lgcp_region.stan"
                                 fname <- "lgcp_region"
                               } else {
                                 filen <- "lgcp.stan"
                                 fname <- "lgcp"
                               }
                             }

                             if(use_cmdstanr){
                               if(!requireNamespace("cmdstanr"))stop("cmdstanr not available.")
                               model_file <- system.file("cmdstan",
                                                         filen,
                                                         package = "rts2",
                                                         mustWork = TRUE)
                               #used for debugging without installing package
                               #model_file <- paste0("D:/Documents/R/rts2/inst/cmdstan/",filen)
                               model <- cmdstanr::cmdstan_model(model_file)
                               if(!vb){
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
                                 res <- model$variational(
                                   data=datlist,
                                   output_dir=dir,
                                   ...
                                 )
                               }
                               
                             } else {
                               if(!verbose){
                                 if(!vb){
                                   capture.output(suppressWarnings( res <- rstan::sampling(stanmodels[[fname]],
                                                                                           data=datlist,
                                                                                           chains=chains,
                                                                                           iter = iter_warmup+iter_sampling,
                                                                                           warmup = iter_warmup,
                                                                                           cores = parallel_chains,
                                                                                           refresh = 0)), file=tempfile())
                                 } else {
                                   capture.output(suppressWarnings( res <- rstan::vb(stanmodels[[fname]],
                                                                                           data=datlist)), file=tempfile())
                                 }
                                 
                               } else {
                                 if(!vb){
                                   res <- rstan::sampling(stanmodels[[fname]],
                                                          data=datlist,
                                                          chains=chains,
                                                          iter = iter_warmup+iter_sampling,
                                                          warmup = iter_warmup,
                                                          cores = parallel_chains)
                                 } else {
                                   res <- rstan::vb(stanmodels[[fname]],
                                                    data=datlist)
                                 }
                                 
                               }

                             }


                             return(res)
                           },
                           #' @description
                           #' Extract predictions
                           #'
                           #' Extract incidence and relative risk predictions
                           #'
                           #' @param fit A \link[rstan]{stanfit}, \link[cmdstanr]{CmdStanMCMC}, \link[cmdstanr]{CmdStanVB} object.
                           #' Output of `lgcp_fit()` or the output of `lgcp_fit_ml()` or `lgcp_fit_la()`
                           #' @param type Vector of character strings. Any combination of "pred", "rr", and "irr", which are,
                           #' posterior mean incidence (overall and population standardised), relative risk,
                           #' and incidence rate ratio, respectively.
                           #' @param irr.lag integer. If "irr" is requested as `type` then the number of time
                           #' periods lag previous the ratio is in comparison to
                           #' @param t.lag integer. Extract predictions for previous time periods.
                           #' @param popdens character string. Name of the column in `grid_data` with the
                           #' population density data
                           #' @param verbose Logical indicating whether to print messages to the console
                           #' @return NULL
                           #' @details
                           #' *EXTRACTING PREDICTIONS*
                           #' Three outputs can be extracted from the model fit, which will be added as columns to `grid_data`:
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
                           #' g1 <- grid$new(b1,0.5)
                           #' dp <- data.frame(y=runif(10,0,3),x=runif(10,0,3),date=paste0("2021-01-",11:20))
                           #' dp <- create_points(dp,pos_vars = c('y','x'),t_var='date')
                           #' cov1 <- grid$new(b1,0.8)
                           #' cov1$grid_data$cov <- runif(nrow(cov1$grid_data))
                           #' g1$add_covariates(cov1,
                           #'                   zcols="cov",
                           #'                   verbose = FALSE)
                           #' g1$points_to_grid(dp, laglength=5)
                           #' g1$priors <- list(
                           #'   prior_lscale=c(0,0.5),
                           #'   prior_var=c(0,0.5),
                           #'   prior_linpred_mean=c(0),
                           #'   prior_linpred_sd=c(5)
                           #'   )
                           #' res <- g1$lgcp_fit(popdens="cov")
                           #' g1$extract_preds(res,
                           #'                  type=c("pred","rr"),
                           #'                  popdens="cov")
                           #' }
                           #' @importFrom stats sd
                           extract_preds = function(fit,
                                                    type=c("pred","rr","irr"),
                                                    irr.lag=NULL,
                                                    t.lag=0,
                                                    popdens=NULL,
                                                    verbose = TRUE){

                             if("irr"%in%type&is.null(irr.lag))stop("For irr set irr.lag")
                             if(!(is(fit,"CmdStanMCMC")|is(fit,"stanfit")|
                                  is(fit,"CmdStanVB")|is(fit,"lgcp_mcmcml")|
                                  is(fit,"lgcp_la")))stop("stan fit or MCMCML fit required")
                             if("pred"%in%type&is.null(popdens))stop("set popdens for pred")

                             nCells <- nrow(self$grid_data)
                             if(!is.null(self$region_data))nRegion <- nrow(self$region_data)
                             if(is(fit,"stanfit")){
                               ypred <- rstan::extract(fit,"y_grid_predict")
                               ypred <- ypred$y_grid_predict
                               f <- rstan::extract(fit,"f")
                               f <- f$f
                               nT <- dim(ypred)[2]/nCells
                               cmdst <- FALSE
                             } else if(is(fit,"CmdStanMCMC")|is(fit,"CmdStanVB")){
                               if(requireNamespace("cmdstanr")){
                                 ypred <- fit$draws("y_grid_predict")
                                 f <- fit$draws("f")
                                 if(length(dim(ypred))==2){
                                   #to convert to 3d if VB is used
                                   ypred <- array(drop(ypred),dim = c(1,dim(ypred)))
                                   f <- array(drop(f),dim = c(1,dim(f)))
                                 }
                                 nT <- dim(ypred)[3]/nCells
                                 cmdst <- TRUE
                               } else {
                                 stop("No cmdstanr package")
                               }
                             } else if(is(fit,"lgcp_mcmcml")|is(fit,"lgcp_la")){
                               ypred <- t(fit$y_grid_predict)
                               f <- t(fit$u)
                               nT <- ncol(f)/nCells
                               cmdst <- FALSE
                             }

                             #print(nT)
                             #print(nCells)

                             if(nT>1){
                               if("pred"%in%type){
                                 if(is.null(self$region_data)){
                                   popd <- as.data.frame(self$grid_data)[,popdens]
                                   if(!cmdst){
                                     fmu <- ypred[,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells),drop=FALSE]/popd
                                     self$grid_data$pred_mean_total <- apply(ypred[,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells),drop=FALSE],2,mean)
                                     self$grid_data$pred_mean_total_sd <- apply(ypred[,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells),drop=FALSE],2,sd)
                                     self$grid_data$pred_mean_pp <- apply(fmu,2,mean)
                                     self$grid_data$pred_mean_pp_sd <- apply(fmu,2,sd)
                                   } else {
                                     fmu <- ypred[,,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells),drop=FALSE]/popd
                                     self$grid_data$pred_mean_total <- apply(ypred[,,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells),drop=FALSE],3,mean)
                                     self$grid_data$pred_mean_total_sd <- apply(ypred[,,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells),drop=FALSE],3,sd)
                                     self$grid_data$pred_mean_pp <- apply(fmu,3,mean)
                                     self$grid_data$pred_mean_pp_sd <- apply(fmu,3,sd)
                                   }
                                 } else {
                                   if(verbose)message("Predicted rates are added to region_data, rr and irr are added to grid_data")
                                   popd <- as.data.frame(self$region_data)[,popdens]
                                   if(!cmdst){
                                     fmu <- ypred[,((nT-1-t.lag)*nRegion+1):((nT-t.lag)*nRegion),drop=FALSE]/popd
                                     self$region_data$pred_mean_total <- apply(ypred[,((nT-1-t.lag)*nRegion+1):((nT-t.lag)*nRegion),drop=FALSE],2,mean)
                                     self$region_data$pred_mean_total_sd <- apply(ypred[,((nT-1-t.lag)*nRegion+1):((nT-t.lag)*nRegion),drop=FALSE],2,sd)
                                     self$region_data$pred_mean_pp <- apply(fmu,2,mean)
                                     self$region_data$pred_mean_pp_sd <- apply(fmu,2,sd)
                                   } else {
                                     fmu <- ypred[,,((nT-1-t.lag)*nRegion+1):((nT-t.lag)*nRegion),drop=FALSE]/popd
                                     self$region_data$pred_mean_total <- apply(ypred[,,((nT-1-t.lag)*nRegion+1):((nT-t.lag)*nRegion),drop=FALSE],3,mean)
                                     self$region_data$pred_mean_total_sd <- apply(ypred[,,((nT-1-t.lag)*nRegion+1):((nT-t.lag)*nRegion),drop=FALSE],3,sd)
                                     self$region_data$pred_mean_pp <- apply(fmu,3,mean)
                                     self$region_data$pred_mean_pp_sd <- apply(fmu,3,sd)
                                   }
                                 }
                                 

                               }

                               if("rr"%in%type){
                                 if(!cmdst){
                                   self$grid_data$rr <- exp(apply(f[,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells),drop=FALSE],2,mean))
                                   self$grid_data$rr_sd <- exp(apply(f[,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells),drop=FALSE],2,sd))
                                 } else {
                                   self$grid_data$rr <- exp(apply(f[,,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells),drop=FALSE],3,mean))
                                   self$grid_data$rr_sd <- exp(apply(f[,,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells),drop=FALSE],3,sd))
                                 }

                               }

                               if("irr"%in%type){
                                 if(!cmdst){
                                   self$grid_data$irr <- apply(ypred[,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells),drop=FALSE]/
                                                            ypred[,((nT-irr.lag-t.lag)*nCells+1):(((nT-t.lag)-irr.lag+1)*nCells),drop=FALSE],2,mean)
                                   self$grid_data$irr_sd <- apply(ypred[,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells),drop=FALSE]/
                                                               ypred[,((nT-irr.lag-t.lag)*nCells+1):(((nT-t.lag)-irr.lag+1)*nCells),drop=FALSE],2,sd)
                                 } else {
                                   self$grid_data$irr <- apply(ypred[,,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells),drop=FALSE]/
                                                            ypred[,,((nT-irr.lag-t.lag)*nCells+1):(((nT-t.lag)-irr.lag+1)*nCells),drop=FALSE],3,mean)
                                   self$grid_data$irr_sd <- apply(ypred[,,((nT-1-t.lag)*nCells+1):((nT-t.lag)*nCells),drop=FALSE]/
                                                               ypred[,,((nT-irr.lag-t.lag)*nCells+1):(((nT-t.lag)-irr.lag+1)*nCells),drop=FALSE],3,sd)
                                 }

                               }
                             } else {

                               if("irr"%in%type)stop("cannot estimate irr as only one time period")

                               if("pred"%in%type){
                                 if(is.null(self$region_data)){
                                   if(!cmdst){
                                     fmu <- ypred/as.data.frame(self$grid_data)[,popdens]
                                     self$grid_data$pred_mean_total <- apply(ypred,2,mean)
                                     self$grid_data$pred_mean_total_sd <- apply(ypred,2,sd)
                                     self$grid_data$pred_mean_pp <- apply(fmu,2,mean)
                                     self$grid_data$pred_mean_pp_sd <- apply(fmu,2,sd)
                                   } else {
                                     fmu <- ypred/as.data.frame(self$grid_data)[,popdens]
                                     self$grid_data$pred_mean_total <- apply(ypred,3,mean)
                                     self$grid_data$pred_mean_total_sd <- apply(ypred,3,sd)
                                     self$grid_data$pred_mean_pp <- apply(fmu,3,mean)
                                     self$grid_data$pred_mean_pp_sd <- apply(fmu,3,sd)
                                   }
                                 } else {
                                   if(!cmdst){
                                     fmu <- ypred/as.data.frame(self$region_data)[,popdens]
                                     self$region_data$pred_mean_total <- apply(ypred,2,mean)
                                     self$region_data$pred_mean_total_sd <- apply(ypred,2,sd)
                                     self$region_data$pred_mean_pp <- apply(fmu,2,mean)
                                     self$region_data$pred_mean_pp_sd <- apply(fmu,2,sd)
                                   } else {
                                     fmu <- ypred/as.data.frame(self$region_data)[,popdens]
                                     self$region_data$pred_mean_total <- apply(ypred,3,mean)
                                     self$region_data$pred_mean_total_sd <- apply(ypred,3,sd)
                                     self$region_data$pred_mean_pp <- apply(fmu,3,mean)
                                     self$region_data$pred_mean_pp_sd <- apply(fmu,3,sd)
                                   }
                                 }
                                 

                               }

                               if("rr"%in%type){
                                 if(!cmdst){
                                   self$grid_data$rr <- exp(apply(f,2,mean))
                                   self$grid_data$rr_sd <- exp(apply(f,2,sd))
                                 } else {
                                   self$grid_data$rr <- exp(apply(f,3,mean))
                                   self$grid_data$rr_sd <- exp(apply(f,3,sd))
                                 }

                               }


                             }


                           },
                           #' @description
                           #' Hotspots
                           #'
                           #' Generate hotspot probabilities
                           #'
                           #' Given a definition of a hotspot in terms of threshold(s) for incidence,
                           #' relative risk, and/or incidence rate ratio, returns the probabilities
                           #' each area is a "hotspot". See Details of `extract_preds`. Columns
                           #' will be added to `grid_data`. Note that for incidence threshold, the threshold should
                           #' be specified as the per individual incidence.
                           #'
                           #' @param stan_fit A \link[rstan]{stanfit} or \link[cmdstanr]{CmdStanMCMC} object.
                           #' Output of `lgcp_fit()`
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
                           #' @return NULL
                           #' @examples
                           #' \dontrun{
                           #' b1 <- sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
                           #' g1 <- grid$new(b1,0.5)
                           #' dp <- data.frame(y=runif(10,0,3),x=runif(10,0,3),date=paste0("2021-01-",11:20))
                           #' dp <- create_points(dp,pos_vars = c('y','x'),t_var='date')
                           #' cov1 <- grid$new(b1,0.8)
                           #' cov1$grid_data$cov <- runif(nrow(cov1$grid_data))
                           #' g1$add_covariates(cov1,
                           #'                   zcols="cov",
                           #'                   verbose = FALSE)
                           #' g1$points_to_grid(dp, laglength=5)
                           #' g1$priors <- list(
                           #'   prior_lscale=c(0,0.5),
                           #'   prior_var=c(0,0.5),
                           #'   prior_linpred_mean=c(0),
                           #'   prior_linpred_sd=c(5)
                           #'   )
                           #' res <- g1$lgcp_fit(popdens="cov")
                           #' g1$hotspots(res,
                           #'             incidence.threshold=1,
                           #'             popdens="cov")
                           #' }
                           hotspots = function(stan_fit,
                                               incidence.threshold=NULL,
                                               irr.threshold=NULL,
                                               irr.lag=NULL,
                                               rr.threshold=NULL,
                                               popdens,
                                               col_label=NULL){

                             if(all(is.null(incidence.threshold),is.null(irr.threshold),is.null(rr.threshold)))stop("At least one criterion required.")
                             if(!is.null(irr.threshold)&is.null(irr.lag))stop("irr.lag must be set")
                             if(!(is(stan_fit,"CmdStanMCMC")|is(stan_fit,"stanfit")))stop("stan fit required")

                             nCells <- nrow(self$grid_data)
                             nT <- sum(grepl("\\bt[0-9]",colnames(self$grid_data)))
                             if(nT == 0)nT <- 1
                             if(is(stan_fit,"stanfit")){
                               ypred <- rstan::extract(stan_fit,"y_grid_predict")$y_grid_predict
                               f <- rstan::extract(stan_fit,"f")$f
                               f <- f[,((nT-1)*nCells+1):(nT*nCells),drop=FALSE]
                             } else if(is(stan_fit,"CmdStanMCMC")){
                               if(requireNamespace("cmdstanr")){
                                 ypred <- stan_fit$draws("y_grid_predict")
                                 ypred <- matrix(ypred, prod(dim(ypred)[1:2]), dim(ypred)[3])
                                 f <- stan_fit$draws("f")
                                 f <- f[,,((nT-1)*nCells+1):(nT*nCells),drop=FALSE]
                                 f <- matrix(f, prod(dim(f)[1:2]), dim(f)[3])
                               }
                             }
                             
                             #nT <- dim(ypred)[3]/nCells

                             nCr <- sum(c(!is.null(incidence.threshold),
                                          !is.null(irr.threshold),
                                          !is.null(rr.threshold)))


                             inc1 <- matrix(0,nrow=nrow(f),ncol=nCells)

                             if(!is.null(incidence.threshold)){
                               fmu <- matrix(0,nrow=nrow(f),ncol=nCells)
                               for(i in 1:nCells){
                                 fmu[,i] <- ypred[,((nT-1)*nCells+i)]/as.data.frame(self$grid_data)[i,popdens]
                               }
                               inc1 <- inc1 + I(fmu > incidence.threshold)*1

                             }

                             if(!is.null(irr.threshold)){
                               if(nT==1)stop("cannot estimate irr as only one time period") else {
                                 inc1 <- inc1 + I(ypred[,((nT-1)*nCells+1):(nT*nCells),drop=FALSE]/
                                                    ypred[,((nT-irr.lag)*nCells+1):((nT-irr.lag+1)*nCells),drop=FALSE] > irr.threshold)*1
                               }

                             }

                             if(!is.null(rr.threshold)){
                               inc1 <- inc1 + I(exp(f) > rr.threshold)*1
                             }

                             inc1 <- I(inc1 == nCr)*1
                             self$grid_data$hotspot_prob <- apply(inc1,2,mean)

                             if(!is.null(col_label)){
                               colnames(self$grid_data)[length(colnames(self$grid_data))] <- col_label
                             }

                           },
                           #' @description
                           #' Aggregate output
                           #'
                           #' Aggregate `lgcp_fit` output to another geography
                           #'
                           #' @param new_geom sf object. A set of polygons covering the same area as `boundary`
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
                           #' g1 <- grid$new(b1,0.5)
                           #' dp <- data.frame(y=runif(10,0,3),x=runif(10,0,3),date=paste0("2021-01-",11:20))
                           #' dp <- create_points(dp,pos_vars = c('y','x'),t_var='date')
                           #' cov1 <- grid$new(b1,0.8)
                           #' cov1$grid_data$cov <- runif(nrow(cov1$grid_data))
                           #' g1$add_covariates(cov1,
                           #'                   zcols="cov",
                           #'                   verbose = FALSE)
                           #' g1$points_to_grid(dp, laglength=5)
                           #' g1$priors <- list(
                           #'   prior_lscale=c(0,0.5),
                           #'   prior_var=c(0,0.5),
                           #'   prior_linpred_mean=c(0),
                           #'   prior_linpred_sd=c(5)
                           #'   )
                           #' res <- g1$lgcp_fit(popdens="cov")
                           #' g1$extract_preds(res,
                           #'                  type=c("pred","rr"),
                           #'                  popdens="cov")
                           #' new1 <- g1$aggregate_output(cov1,
                           #'                             zcols="rr")
                           #' }
                           aggregate_output = function(new_geom,
                                                       zcols,
                                                       weight_type="area",
                                                       popdens=NULL,
                                                       verbose=TRUE){

                             if(sf::st_crs(self$grid_data)!=sf::st_crs(new_geom)){
                               sf::st_crs(self$grid_data) <- sf::st_crs(new_geom)
                               warning("st_crs(self$grid_data)!=st_crs(new_geom) setting equal")
                             }

                             tmp <- sf::st_intersection(self$grid_data[,zcols],new_geom)
                             tmp_len <- lengths(sf::st_intersects(new_geom,self$grid_data))
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

                           },
                           #' @description
                           #' Returns scale conversion factor
                           #'
                           #' Coordinates are scaled to `[-1,1]` for `lgcp_fit()`. This function
                           #' returns the scaling factor for this conversion.
                           #'
                           #' @return numeric
                           #' @examples
                           #' b1 = sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
                           #' g1 <- grid$new(b1,0.5)
                           #' g1$scale_conversion_factor()
                           scale_conversion_factor = function(){
                             x_grid <- as.data.frame(suppressWarnings(sf::st_coordinates(sf::st_centroid(self$grid_data))))
                             xrange <- range(x_grid[,1])
                             yrange <- range(x_grid[,2])
                             std_val <- max(max(xrange - mean(xrange)),max(yrange - mean(yrange)))
                             return(std_val)
                           },
                           #' @description
                           #' Fit the LGCP using Markov Chain Monte Carlo Maximum likelihood
                           #' 
                           #' Various Markov Chain Monte Carlo Maximum Likelihood algorithms are provided for the full and 
                           #' regional models, and using a full LGCP or nearest neighbour approximation.
                           #' 
                           #' @details
                           #' *MAXIMUM LIKELIHOOD MODEL FITTING*
                           #' The arguments `mcml_options` and `la_options` for the functions `lgcp_fit_ml` and `lgcp_fit_la` are named lists 
                           #' with the options of whether to use NNGP (`useNN`), whether to use 
                           #' Newton-Raphson (`mcnr`), whether to treat the covariance parameters as known (`known_theta`), whether to provide
                           #' more detailed output (`trace`), the number of nearest neighbours if using NNGP (`nNN`), the tolerance for when to 
                           #' terminate the algorithm (`tol`), and the maximum number of algorithm iterations (`maxiter`). The argument 
                           #' `mcmc_options` is a named list with the number of warmup and sampling iterations for the MCMC sampler.
                           #' @param popdens character string. Name of the population density column
                           #' @param covs vector of character string. Base names of the covariates to
                           #' include. For temporally-varying covariates only the stem is required and not
                           #' the individual column names for each time period (e.g. `dayMon` and not `dayMon1`,
                           #' `dayMon2`, etc.)
                           #' @param start Starting values of the model parameteters in the order c(beta, theta, rho). Rho is optional for models with nT>1
                           #' @param model Either "exp" for exponential covariance function or "sqexp" for squared exponential
                           #' covariance function
                           #' @param mcml_options List of options for the algorithm. See details.
                           #' @param mcmc_options List of options for the MCMC sampling. See details.
                           #' @param verbose logical. Provide feedback on progress
                           #' @param use_cmdstanr logical. Defaults to false. If true then cmdstanr will be used
                           #' instead of rstan.
                           #' @return A named list with the estimated parameters, number of iterations, whether the algorithm converged, and 
                           #' the random effect samples.
                           lgcp_fit_ml = function(popdens,
                                                  covs=NULL,
                                                  start,
                                                  model = "exp",
                                                  mcml_options = list(useNN=FALSE,mcnr=FALSE,known_theta=FALSE,trace=1,nNN=10,tol=1e-2, maxiter=10),
                                                  mcmc_options = list(warmup = 100, sampling = 100),
                                                  verbose=TRUE,
                                                  use_cmdstanr = FALSE
                                                  ){
                             
                             if(mcml_options$nNN<1)stop("number of nearest neighbours must be positive")
                             if(!all(c("useNN","mcnr","known_theta","trace","nNN","tol","maxiter")%in%names(mcml_options)))stop("mcml options does not have all options specified")
                             if(!all(c("warmup","sampling")%in%names(mcmc_options)))stop("mcmc options does not have all options specified")
                             
                             if(!is.null(self$region_data)&verbose)message("Using regional data model.")
                             approx <- ifelse(mcml_options$useNN,"nngp","none")
                             datlist <- private$prepare_data(mcml_options$nNN,model,approx,popdens,covs,verbose,FALSE)
                             npar <- ncol(datlist$X) + 2
                             if(datlist$nT > 1)npar = npar + 1
                             if(length(start) != npar)stop("Wrong number of starting values")
                             
                             
                             args <- list(y= datlist$y,
                                          X = datlist$X,
                                          coords = datlist$x_grid,
                                          popdens = log(datlist$popdens),
                                          nT = datlist$nT,
                                          start = start,
                                          mod = model,
                                          mcml_options = mcml_options,
                                          mcmc_options = mcmc_options,
                                          verbose = verbose,
                                          use_cmdstanr = use_cmdstanr)
                             
                             if(!is.null(self$region_data)){
                               regiondat <- list(ncell = datlist$n_cell,
                                                   cell_id = datlist$cell_id,
                                                   q_weights = datlist$q_weights)
                               args <- append(args,list(region_data = regiondat))
                             }
                             
                             out <- do.call("lgcp_mcmcml",args)
                             
                             xb <- datlist$X%*%out$beta
                             zu <- get_ZLu(diag(nrow(datlist$x_grid)),u = out$u,nT = datlist$nT,rho = out$rho)
                             
                             if(is.null(self$region_data)){
                               y <- zu
                               for(i in 1:ncol(y)){
                                 y[,i] <- exp(y[,i] + xb + log(datlist$popdens))
                               } 
                             } else {
                               y <- region_intensity(xb = xb,y = datlist$y,zu = zu,
                                                     offset = log(datlist$popdens),nT = datlist$nT,
                                                     n_cell = regiondat$ncell,cell_id = regiondat$cell_id,
                                                     q_weights = regiondat$q_weights)
                             }
                             
                             out$y_grid_predict <- y
                             
                             
                             class(out) <- "lgcp_mcmcml"
                             return(out)
                           },
                           #' @description 
                           #' Fit the LGCP using Markov Chain Monte Carlo Maximum likelihood
                           #' 
                           #' Various Markov Chain Monte Carlo Maximum Likelihood algorithms are provided for the full and 
                           #' regional models, and using a full LGCP or nearest neighbour approximation.
                           #' 
                           #' @param popdens character string. Name of the population density column
                           #' @param covs vector of character string. Base names of the covariates to
                           #' include. For temporally-varying covariates only the stem is required and not
                           #' the individual column names for each time period (e.g. `dayMon` and not `dayMon1`,
                           #' `dayMon2`, etc.)
                           #' @param start Starting values of the model parameteters in the order c(beta, theta, rho). Rho is optional for models with nT>1
                           #' @param model Either "exp" for exponential covariance function or "sqexp" for squared exponential
                           #' covariance function
                           #' @param la_options List of options for the algorithm. See details.
                           #' @param verbose logical. Provide feedback on progress
                           #' @param use_cmdstanr logical. Defaults to false. If true then cmdstanr will be used
                           #' instead of rstan to sample random effects at the end of the algorithm
                           #' @return A named list with the estimated parameters, number of iterations, whether the algorithm converged, and 
                           #' the random effect samples.
                           lgcp_fit_la = function(popdens,
                                                  covs=NULL,
                                                  start,
                                                  model = "exp",
                                                  la_options = list(useNN=FALSE,known_theta=FALSE,nr=FALSE,trace=1,nNN=10,tol=1e-2, maxiter=10),
                                                  verbose=TRUE,
                                                  use_cmdstanr = TRUE
                           ){
                             
                             if(!is.null(self$region_data)&verbose)stop("Laplace approximation model fitting not available for regional data model.")
                             if(la_options$nNN<1)stop("Number of nearest neighbours must be positive")
                             if(!all(c("useNN","known_theta","nr","trace","nNN","tol","maxiter")%in%names(la_options)))stop("la options does not have all options specified")
                             
                             approx <- ifelse(la_options$useNN,"nngp","none")
                             datlist <- private$prepare_data(la_options$nNN,model,approx,popdens,covs,verbose,FALSE)
                             npar <- ncol(datlist$X) + 2
                             if(datlist$nT > 1)npar = npar + 1
                             if(length(start) != npar)stop("Wrong number of starting values")
                             
                             
                             args <- list(y= datlist$y,
                                          X = datlist$X,
                                          coords = datlist$x_grid,
                                          popdens = log(datlist$popdens),
                                          nT = datlist$nT,
                                          start = start,
                                          mod = model,
                                          la_options = la_options,
                                          verbose = verbose)
                             
                             out <- do.call("lgcp_la",args)
                             ## sample the random effects
                             coords <- as.data.frame(datlist$x_grid)
                             colnames(coords) <- c('X','Y')
                             
                             f1 <- ifelse(model=="exp","~(1|fexp(X,Y))","~(1|sqexp(X,Y))")
                             
                             cov <- glmmrBase::Covariance$new(
                               formula =formula(f1),
                               parameters = out$theta,
                               data=coords
                             )
                             
                             
                             L <- cov$get_chol_D()
                             xb <- datlist$X%*%out$beta + log(datlist$popdens)
                             v <- sample_u(y = datlist$y, xb =xb, coords = coords,
                                           nT = datlist$nT,
                                           start = c(out$theta,out$rho),verbose=verbose,
                                           use_cmdstanr = use_cmdstanr)
                             v <- get_ZLu(L,v,datlist$nT,out$rho)
                             out$u_la <- out$u
                             out$u <- v
                             y <- v
                             for(i in 1:ncol(y)){
                               y[,i] <- exp(y[,i] + xb)
                             } 
                             
                             out$y_grid_predict <- y
                             
                             class(out) <- "lgcp_la"
                             return(out)
                           },
                           #' @description 
                           #' Returns summary data of the region/grid intersections
                           #' 
                           #' Information on the intersection between the region areas and the computational grid
                           #' including the number of cells intersecting each region (`n_cell`), the indexes of the
                           #' cells intersecting each region in order (`cell_id`), and the proportion of each region's 
                           #' area covered by each intersecting grid cell (`q_weights`).
                           #' @return A named list
                           get_region_data = function(){
                             if(is.null(self$region_data))stop("No region data")
                             ncell <- unname(table(private$intersection_data$region_id))
                             ncell <- c(1, cumsum(ncell)+1)
                             datlist <- list(
                               n_region = nrow(self$region_data),
                               n_Q = nrow(private$intersection_data),
                               n_cell = ncell,
                               cell_id = private$intersection_data$grid_id,
                               q_weights = private$intersection_data$w
                             )
                             return(datlist)
                           },
                           #' @description 
                           #' Plots the empirical semi-variogram
                           #' 
                           #' @param popdens String naming the variable in the data specifying the offset. If not 
                           #' provided then no offset is used.
                           #' @param nbins The number of bins in the empirical semivariogram
                           #' @return A ggplot plot is printed and optionally returned
                           variogram = function(popdens,
                                                yvar,
                                                nbins = 20){
                             
                             if(is.null(self$region_data)){
                               
                               if(!missing(popdens)){
                                 if(!popdens%in%colnames(self$grid_data))stop("popdens variable not in grid data")
                                 offs <- as.data.frame(self$grid_data)[,popdens]
                               } else {
                                 offs <- rep(1,nrow(self$grid_data))
                               }
                               
                               if(missing(yvar)){
                                 nT <- sum(grepl("\\bt[0-9]",colnames(self$grid_data)))
                                 if(nT==0){
                                   if("y"%in%colnames(self$grid_data)){
                                     yvar <- "y"
                                   } else {
                                     stop("case counts not defined in data")
                                   }
                                 } else {
                                   yvar <- paste0("t",nT)
                                 }
                               } else {
                                 if(!yvar%in%colnames(self$grid_data))stop("yvar not in grid data")
                               }
                               
                               
                               df <- suppressWarnings(as.data.frame(sf::st_coordinates(sf::st_centroid(self$grid_data))))
                               dfs <- semivariogram(as.matrix(df),
                                                    offs = as.data.frame(self$grid_data)[,popdens],
                                                    y = as.data.frame(self$grid_data)[,yvar],nbins)
                               dfs <- as.data.frame(dfs)
                               colnames(dfs) <- c("bin","val")
                               p <- ggplot2::ggplot(data=dfs,ggplot2::aes(x=bin,y=val))+
                                 ggplot2::geom_point()+
                                 ggplot2::theme_bw()+
                                 ggplot2::theme(panel.grid=ggplot2::element_blank())+
                                 ggplot2::labs(x="Distance",y="Semivariogram function")
                               print(p)
                               return(invisible(p))
                             } else {
                               
                               if(!missing(popdens)){
                                 if(!popdens%in%colnames(self$region_data))stop("popdens variable not in region data")
                                 offs <- as.data.frame(g1$region_data)[,popdens]
                               } else {
                                 offs <- rep(1,nrow(self$region_data))
                               }
                               
                               if(missing(yvar)){
                                 nT <- sum(grepl("\\bt[0-9]",colnames(self$region_data)))
                                 if(nT==0){
                                   if("y"%in%colnames(self$region_data)){
                                     yvar <- "y"
                                   } else {
                                     stop("case counts not defined in data")
                                   }
                                 } else {
                                   yvar <- paste0("t",nT)
                                 }
                               } else {
                                 if(!yvar%in%colnames(self$region_data))stop("yvar not in region data")
                               }
                               
                               message("Using centroids of regions as locations")
                               df <- suppressWarnings(as.data.frame(sf::st_coordinates(sf::st_centroid(self$region_data))))
                               dfs <- semivariogram(as.matrix(df),
                                                    offs = as.data.frame(self$region_data)[,popdens],
                                                    y = as.data.frame(self$region_data)[,yvar],nbins)
                               dfs <- as.data.frame(dfs)
                               colnames(dfs) <- c("bin","val")
                               p <- ggplot2::ggplot(data=dfs,ggplot2::aes(x=bin,y=val))+
                                 ggplot2::geom_point()+
                                 ggplot2::theme_bw()+
                                 ggplot2::theme(panel.grid=ggplot2::element_blank())+
                                 ggplot2::labs(x="Distance",y="Semivariogram function")
                               print(p)
                               return(invisible(p))
                             }
                             
                           },
                           #' @description 
                           #' Re-orders the computational grid
                           #' 
                           #' The quality of the nearest neighbour approximation can depend on the ordering of 
                           #' the grid cells. This function reorders the grid cells.
                           #' @param option Either "y" for order of the y coordinate, "x" for order of the x coordinate,
                           #' "minimax"  in which the next observation in the order is the one which maximises the
                           #'  minimum distance to the previous observations,g1$grid_data <- g1$grid_data[o0,] or "random" which randomly orders them.
                           #'  @return No return, used for effects.
                           reorder = function(option="y"){
                             df <- suppressWarnings(as.data.frame(sf::st_coordinates(sf::st_centroid(self$grid_data))))
                             colnames(df) <- c("x","y")
                             if(option=="y"){
                               o0 <- order(df$y,df$x,decreasing = FALSE)
                               self$grid_data <- self$grid_data[o0,]
                             } else if(option=="x"){
                               o0 <- order(df$x,df$y,decreasing = FALSE)
                               self$grid_data <- self$grid_data[o0,]
                             } else if(option=="minimax"){
                               o1 <- rep(NA,nrow(df))
                               o1[which.min(sqrt((df$x-0.5)^2 + (df$y - 0.5)^2))] <- 1
                               for(i in 2:nrow(df)){
                                 dists <- matrix(0,nrow=nrow(df)-i+1, ncol=(i-1))
                                 for(j in 1:(i-1)){
                                   dists[,j] <- sqrt((df$x[is.na(o1)] - df$x[which(o1 == j)])^2 +
                                                       (df$y[is.na(o1)] - df$y[which(o1 == j)])^2)
                                 }
                                 mindists <- apply(dists,1,min)
                                 o1[is.na(o1)][which.max(mindists)] <- i
                               }
                               self$grid_data <- self$grid_data[o1,]
                             } else if(option=="random"){
                               o2 <- sample(1:nrow(df),nrow(df),replace=FALSE)
                               self$grid_data <- self$grid_data[o2,]
                             }
                           }
                         ),
                    private = list(
                      intersection_data = NULL,
                      prepare_data = function(m,
                                              model,
                                              approx,
                                              popdens,
                                              covs,
                                              verbose,
                                              bayes){
                        
                        mod <- NA
                        if(model == "exp"){
                          mod <- 1
                        } else if(model=="sqexp"){
                          mod <- 0
                        }
                        
                        ind <- as.matrix(expand.grid(1:m,1:m))
                        if(is.null(self$region_data)){
                          nT <- sum(grepl("\\bt[0-9]",colnames(self$grid_data)))
                          if(nT==0){
                            if("y"%in%colnames(self$grid_data)){
                              nT <- 1
                            } else {
                              stop("case counts not defined in data")
                            }
                          }
                        } else {
                          nT <- sum(grepl("\\bt[0-9]",colnames(self$region_data)))
                          if(nT==0){
                            if("y"%in%colnames(self$region_data)){
                              nT <- 1
                            } else {
                              stop("case counts not defined in data")
                            }
                          }
                        }
                        
                        nCell <- nrow(self$grid_data)
                        
                        x_grid <- as.data.frame(suppressWarnings(sf::st_coordinates(
                          sf::st_centroid(self$grid_data))))
                        
                        
                        if(approx=="hsgp"){
                          # scale to -1,1 in all dimensions
                          xrange <- range(x_grid[,1])
                          yrange <- range(x_grid[,2])
                          std_val <- max(max(xrange - mean(xrange)),max(yrange - mean(yrange)))
                          
                          x_grid[,1] <- (x_grid[,1]- mean(xrange))/std_val
                          x_grid[,2] <- (x_grid[,2]- mean(yrange))/std_val
                        }
                        
                        
                        if(is.null(self$region_data)){
                          # outcome data
                          if(nT > 1){
                            y <- stack(as.data.frame(self$grid_data)[,paste0("t",1:nT)])[,1]
                          } else {
                            y <- as.data.frame(self$grid_data)[,"y"]
                          }
                          
                          #population density
                          nColP <- sum(grepl(popdens,colnames(self$grid_data)))
                          if(nColP==1){
                            popd <- rep(as.data.frame(self$grid_data)[,popdens],nT)
                          } else if(nColP==0){
                            stop("popdens variable not found in grid data")
                          } else {
                            if(nT>1){
                              popd <- stack(as.data.frame(self$grid_data)[,paste0(popdens,1:nT)])[,1]
                            } else {
                              popd <- as.data.frame(self$grid_data)[,popdens]
                            }
                          }
                          
                          #add covariates
                          if(!is.null(covs)){
                            nQ <- length(covs)
                            X <- matrix(NA,nrow=length(y),ncol=nQ+1)
                            X[,1] <- 1
                            for(i in 1:nQ){
                              nColV <- sum(grepl(covs[i],colnames(self$grid_data)))
                              if(nColV==1){
                                X[,i+1] <- rep(as.data.frame(self$grid_data)[,covs[i]],nT)
                              } else if(nColV==0){
                                stop(paste0(covs[i]," not found in grid data"))
                              } else {
                                if(nT>1){
                                  X[,i+1] <- stack(as.data.frame(self$grid_data)[,paste0(covs[i],1:nT)])[,1]
                                } else {
                                  X[,i+1] <- as.data.frame(self$grid_data)[,covs[i]]
                                }
                              }
                              Q <- nQ+1
                            }
                          } else {
                            X <- matrix(1,nrow=length(y),ncol=1)
                            Q <- 1
                          }
                        } else {
                          #outcomes
                          if(nT > 1){
                            y <- stack(as.data.frame(self$region_data)[,paste0("t",1:nT)])[,1]
                          } else {
                            y <- as.data.frame(self$region_data)[,"y"]
                          }
                          
                          #population density
                          nColP <- sum(grepl(popdens,colnames(self$region_data)))
                          if(nColP==1){
                            popd <- rep(as.data.frame(self$region_data)[,popdens],nT)
                          } else if(nColP==0){
                            stop("popdens variable not found in region data")
                          } else {
                            if(nT>1){
                              popd <- stack(as.data.frame(self$region_data)[,paste0(popdens,1:nT)])[,1]
                            } else {
                              popd <- as.data.frame(self$region_data)[,popdens]
                            }
                          }
                          
                          #add covariates
                          if(!is.null(covs)){
                            nQ <- length(covs)
                            X <- matrix(NA,nrow=length(y),ncol=nQ+1)
                            X[,1] <- 1
                            for(i in 1:nQ){
                              nColV <- sum(grepl(covs[i],colnames(self$region_data)))
                              if(nColV==1){
                                X[,i+1] <- rep(as.data.frame(self$region_data)[,covs[i]],nT)
                              } else if(nColV==0){
                                stop(paste0(covs[i]," not found in region data"))
                              } else {
                                if(nT>1){
                                  X[,i+1] <- stack(as.data.frame(self$region_data)[,paste0(covs[i],1:nT)])[,1]
                                } else {
                                  X[,i+1] <- as.data.frame(self$region_data)[,covs[i]]
                                }
                              }
                              Q <- nQ+1
                            }
                          } else {
                            X <- matrix(1,nrow=length(y),ncol=1)
                            Q <- 1
                          }
                        }
                        
                        if(bayes){
                          if(!is.null(self$priors)){
                            if(length(self$priors$prior_linpred_mean)!=Q|length(self$priors$prior_linpred_sd)!=Q)
                              stop("Prior mean or sd vector for linear predictior is not equal to number of covariates")
                            if(length(self$priors$prior_lscale)!=2|length(self$priors$prior_var)!=2)
                              stop("prior_lscale or prior_var not of length 2")
                          } else {
                            warning("priors not set, using defaults")
                            self$priors <- list(
                              prior_lscale = c(0,0.5),
                              prior_var = c(0,0.5),
                              prior_linpred_mean = rep(0,Q),
                              prior_linpred_sd = rep(5,Q)
                            )
                          }
                        }
                        
                        
                        if(verbose)message(paste0(nCell," grid cells ",nT," time periods, and ",Q," covariates. Starting sampling..."))
                        if(approx == "hsgp"){
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
                            mod = mod
                          )
                          
                          if(!is.null(self$region_data)){
                            ncell <- unname(table(private$intersection_data$region_id))
                            ncell <- c(1, cumsum(ncell)+1)
                            datlist <- append(datlist,list(
                              n_region = nrow(self$region_data),
                              n_Q = nrow(private$intersection_data),
                              n_cell = ncell,
                              cell_id = private$intersection_data$grid_id,
                              q_weights = private$intersection_data$w
                            ))
                          } 
                          
                          
                        } else if(approx == "nngp"){
                          NN <- suppressWarnings(genNN(sf::st_coordinates(sf::st_centroid(self$grid_data)),m))
                          NN <- NN+1
                          datlist <- list(
                            D = 2,
                            Q = Q,
                            M = m,
                            Nsample = nCell,
                            nT = nT,
                            NN = NN,
                            y = y,
                            x_grid = x_grid[,1:2],
                            popdens = popd,
                            X = X,
                            mod = mod
                          )
                          if(!is.null(self$region_data)){
                            ncell <- unname(table(private$intersection_data$region_id))
                            ncell <- c(1,cumsum(ncell)+1)
                            datlist <- append(datlist,list(
                              n_region = nrow(self$region_data),
                              n_Q = nrow(private$intersection_data),
                              n_cell = ncell,
                              cell_id = private$intersection_data$grid_id,
                              q_weights = private$intersection_data$w
                            ))
                          } 
                          
                          
                        } else {
                          datlist <- list(
                            D = 2,
                            Q = Q,
                            Nsample = nCell,
                            nT = nT,
                            y = y,
                            x_grid = x_grid[,1:2],
                            popdens = popd,
                            X = X,
                            mod = mod
                          )
                          
                          if(!is.null(self$region_data)){
                            ncell <- unname(table(private$intersection_data$region_id))
                            ncell <- c(1,cumsum(ncell)+1)
                            datlist <- append(datlist,list(
                              n_region = nrow(self$region_data),
                              n_Q = nrow(private$intersection_data),
                              n_cell = ncell,
                              cell_id = private$intersection_data$grid_id,
                              q_weights = private$intersection_data$w
                            ))
                          } 
                        }
                        
                        if(bayes){
                          datlist <- append(datlist,list(prior_lscale=self$priors$prior_lscale,
                                                         prior_var=self$priors$prior_var,
                                                         prior_linpred_mean = as.array(self$priors$prior_linpred_mean),
                                                         prior_linpred_sd=as.array(self$priors$prior_linpred_sd)))
                        }
                        
                        return(datlist)
                      }
                    ))

