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
                           #' @field bobyqa_control list of control parameters for the BOBYQA algorithm, must contain named
                           #' elements any or all of `npt`, `rhobeg`, `rhoend`, `covrhobeg`, `covrhoend`. 
                           #' Only has an effect for the HSGP and NNGP approximations. The latter two parameters control the 
                           #' covariance parameter optimisation, while the former control the linear predictor. 
                           bobyqa_control = NULL,
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
                               if(nrow(tmp)==0)stop("Covariate data does not overlap grid")
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
                           #' Adds time period indicators to the data
                           #' 
                           #' Adds indicator variables for each time period to the data. To include
                           #' these in a model fitting procedure use, for example, `covs = c("time1i, time2i,...)`
                           #' @return 
                           #' Nothing. Called for effects.
                           add_time_indicators = function(){
                             if(is.null(self$region_data)){
                               nT <- sum(grepl("\\bt[0-9]",colnames(self$grid_data)))
                               if(nT == 1){
                                 stop("Only one time period")
                               } else {
                                 for(i in 1:nT){
                                   for(j in 1:nT){
                                     if(i == j){
                                       self$grid_data$newvar <- 1
                                     } else {
                                       self$grid_data$newvar <- 0
                                     }
                                     colnames(self$grid_data)[length(colnames(self$grid_data))] <- paste0("time",i,"i",j)
                                   }
                                 }
                               }
                             } else {
                               nT <- sum(grepl("\\bt[0-9]",colnames(self$region_data)))
                               if(nT == 1){
                                 stop("Only one time period")
                               } else {
                                 for(i in 1:nT){
                                   for(j in 1:nT){
                                     if(i == j){
                                       self$region_data$newvar <- 1
                                     } else {
                                       self$region_data$newvar <- 0
                                     }
                                     colnames(self$region_data)[length(colnames(self$region_data))] <- paste0("time",i,"i",j)
                                   }
                                 }
                               }
                             }
                           },
                           #' @description
                           #' Fit an (approximate) log-Gaussian Cox Process model using Bayesian methods
                           #'
                           #' @details
                           #' *BAYESIAN MODEL FITTING*
                           #' The grid data must contain columns `t*`, giving the case
                           #' count in each time period (see `points_to_grid`), as well as any covariates to include in the model
                           #' (see `add_covariates`) and the population density. Otherwise, if the data are regional data, then the outcome
                           #' counts must be in self$region_data
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
                           #' @param covs_grid If using a region model, covariates at the level of the grid can also be specified by providing their
                           #' names to this argument.
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
                           #' res <- g1$lgcp_bayes(popdens="cov")
                           #' }
                           lgcp_bayes = function(popdens,
                                               covs=NULL,
                                               covs_grid = NULL,
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
                                               use_cmdstanr = TRUE,
                                               ...){

                             if(verbose)if(is.null(dir))message("dir not set, files will be lost after session restart")
                             if(!approx%in%c('nngp','hsgp','none'))stop("approx must be one of nngp, hsgp, or none")
                             if(m<=1 & approx %in% c('nngp','hsgp'))stop("m must be greater than one")
                             if(m >25 & verbose)message("m is large, sampling may take a long time.")
                             if(!is.null(self$region_data)&verbose)message("Using regional data model.")
                             #prepare data for model fit
                             datlist <- private$prepare_data(m,model,approx,popdens,covs,covs_grid,verbose,TRUE,L)
                             if(!is.null(known_theta)){
                               if(length(known_theta)!=2)stop("Theta should be of length 2")
                               datlist$known_cov <- 1
                               datlist$sigma_data <- known_theta[1]
                               datlist$phi_data <- known_theta[2]
                             } else {
                               datlist$known_cov <- 0
                               datlist$sigma_data <- c(1)
                               datlist$phi_data <- c(1)
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
                                 filen <- "lgcp_region_cmd.stan"
                                 fname <- "lgcp_region"
                               } else {
                                 filen <- "lgcp_cmd.stan"
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
                           #' Fit an (approximate) log-Gaussian Cox Process model using Maximum Likelihood
                           #'
                           #' @details
                           #' *MAXIMUM LIKELIHOOD MODEL FITTING*
                           #' The grid data must contain columns `t*`, giving the case
                           #' count in each time period (see `points_to_grid`), as well as any covariates to include in the model
                           #' (see `add_covariates`) and the population density. Otherwise, if the data are regional data, then the outcome
                           #' counts must be in self$region_data. See `lgcp_bayes()` for more details on the model.
                           #' 
                           #' The argument `approx` specifies whether to use a full LGCP model (`approx='none'`) or whether
                           #' to use either a nearest neighbour approximation (`approx='nngp'`) 
                           #' 
                           #' Model fitting uses MCMC Maximum likelihood, which has three steps:
                           #'  1. Sample random effects using MCMC. Using cmdstanr is recommended as it is much faster. The arguments 
                           #'    `mcmc_warmup` and `mcmc_sampling` specify the warmup and sampling iterations for this step.
                           #'  2. Fit fixed effect parameters using expectation maximisation.
                           #'  3. Fit covariance parameters using expectation maximisation. This third step is the slowest. The NNGP approximation
                           #'     provides some speed improvements. Otherwise this step can be skipped if the covaraince parameters are "known".   
                           #' 
                           #' @param popdens character string. Name of the population density column
                           #' @param covs vector of strings. Base names of the covariates to
                           #' include. For temporally-varying covariates only the stem is required and not
                           #' the individual column names for each time period (e.g. `dayMon` and not `dayMon1`,
                           #' `dayMon2`, etc.) Alternatively, a formula can be passed to the `formula` arguments below.
                           #' @param covs_grid If using a region model, covariates at the level of the grid can also be specified by providing their
                           #' names to this argument. Alternatively, a formula can be passed to the `formula` arguments below.
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
                           #' @param starting_values An optional list providing starting values of the model parameters. The list can have named elements
                           #' `gamma` for the linear predictor parameters, `theta` for the covariance parameters, and `ar` for the auto-regressive parameter.
                           #' If there are covariates for the grid in a region data model then their parameters are `gamma_g`. The list elements must be a 
                           #' vector of starting values. If this is not provided then the non-intercept linear predictor parameters are initialised randomly
                           #' as N(0,0.1), the covariance parameters as Uniform(0,0.5) and the auto-regressive parameter to 0.1. 
                           #' @param lower_bound Optional. Vector of lower bound values for the fixed effect parameters.
                           #' @param upper_bound Optional. Vector of upper bound values for the fixed effect parameters.
                           #' @param formula_1 Optional. Instead of providing a list of covariates above (to `covs`) a formula can be specified here. For a regional model, this 
                           #' argument specified the regional-level fixed effects model.
                           #' @param formula_2 Optional. Instead of providing a list of covariates above (to `covs_grid`) a formula can be specified here. For a regional model, this 
                           #' argument specified the grid-level fixed effects model.
                           #' @param algo integer. 1 = L-BFGS for beta and non-approximate covariance parameters (default), 2 = BOBYQA for both, 3 = L-BFGS for beta, BOBYQA for covariance parameters.
                           #' @param iter_warmup integer. Number of warmup iterations
                           #' @param iter_sampling integer. Number of sampling iterations
                           #' @param trace Integer. Level of detail of information printed to the console. 0 = none, 1 = some (default), 2 = most.
                           #' @param use_cmdstanr logical. Defaults to false. If true then cmdstanr will be used
                           #' instead of rstan.
                           #' @param ... additional options to pass to `$sample()``, see \link[cmdstanr]{sample}
                           #' @return A `mcmlrts` model fit object
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
                           #' res <- g1$lgcp_ml(popdens="cov")
                           #' }
                           lgcp_ml = function(popdens,
                                              covs=NULL,
                                              covs_grid = NULL,
                                              approx = "nngp",
                                              m=10,
                                              L=1.5,
                                              model = "exp",
                                              known_theta = NULL,
                                              starting_values = NULL,
                                              lower_bound = NULL,
                                              upper_bound = NULL,
                                              formula_1 = NULL,
                                              formula_2 = NULL,
                                              algo = 1,
                                              tol = 1e-2,
                                              max.iter = 30,
                                              iter_warmup=100,
                                              iter_sampling=250,
                                              trace = 1,
                                              use_cmdstanr = TRUE){
                             
                             # some checks at the beginning
                             if(!approx%in%c('nngp','none','hsgp'))stop("approx must be one of nngp, hsgp or none")
                             if(m<=1 & approx == 'nngp')stop("m must be greater than one")
                             if(m >25 & trace >= 1)message("m is large, sampling may take a long time.")
                             if(!is.null(self$region_data)& trace >= 1)message("Using regional data model.")
                             
                             # set up main data and initialise the pointer to the C++ class
                             datlist <- private$update_ptr(m,model,approx,popdens,covs,covs_grid,L,TRUE, formula_1, formula_2)
                             if(!is.null(known_theta)){
                               if((length(known_theta)!=2 & datlist$nT == 1) | (length(known_theta)!=3 & datlist$nT > 1))stop("Theta should be of length 2 (T=1) or 3")
                               datlist$known_cov <- 1
                               datlist$sigma_data <- as.array(known_theta[1])
                               datlist$phi_data <- as.array(known_theta[2])
                               rtsModel__update_theta(private$ptr,known_theta[1:2],private$cov_type,private$lp_type)
                               if(datlist$nT>1)rtsModel__update_rho(private$ptr,known_theta[3],private$cov_type,private$lp_type)
                             } else {
                               datlist$known_cov <- 0
                               datlist$sigma_data <- c()
                               datlist$phi_data <- c()
                             }
                             rtsModel__set_trace(private$ptr,trace,private$cov_type,private$lp_type)
                             
                             ## deal with starting values and initialise parameters
                             if(!is.null(starting_values)){
                               if("gamma"%in%names(starting_values)){
                                 gamma_current <- rtsModel__get_beta(private$ptr,private$cov_type,private$lp_type)
                                 gamma_start <- starting_values[["gamma"]]
                                 if("gamma_g"%in%names(starting_values)){
                                   if(is.null(self$region_data)){
                                     message("gamma_g starting values ignored as no region data")
                                   } else {
                                     gamma_start <- c(gamma_start,starting_values[["gamma_g"]])
                                   }
                                 }
                                 if(length(gamma_start)!=length(gamma_current))stop("Gamma wrong length")
                                 rtsModel__update_beta(private$ptr,gamma_start,private$cov_type,private$lp_type)
                               }
                               if("theta"%in%names(starting_values)){
                                 if(length(starting_values[["theta"]])!=2)stop("length(theta) != 2")
                                 rtsModel__update_theta(private$ptr,starting_values[["theta"]],private$cov_type,private$lp_type)
                               }
                               if("ar"%in%names(starting_values)){
                                 if(datlist$nT == 1){
                                   message("ar ignored as single period")
                                 } else {
                                   rtsModel__update_rho(private$ptr,starting_values[["ar"]],private$cov_type,private$lp_type)
                                 }
                               }
                             }
                             
                             # ## update the BOBYQA control parameters if necessary
                             # if(!is.null(self$bobyqa_control)){
                             #   npt <- ifelse("npt"%in%names(self$bobyqa_control),self$bobyqa_control$npt,0)
                             #   rhobeg <- ifelse("rhobeg"%in%names(self$bobyqa_control),self$bobyqa_control$rhobeg,0)
                             #   rhoend <- ifelse("rhoend"%in%names(self$bobyqa_control),self$bobyqa_control$rhoend,0)
                             #   #if(verbose)cat("\nBOBYQA control parameters: npt(",npt,"), rhobeg(",rhobeg,"), rhoend(",rhoend,")")
                             #   rtsModel__set_bobyqa_control(private$ptr,private$cov_type,private$lp_type,npt,rhobeg,rhoend)
                             #   rhobeg <- ifelse("covrhobeg"%in%names(self$bobyqa_control),self$bobyqa_control$covrhobeg,0)
                             #   rhoend <- ifelse("covrhoend"%in%names(self$bobyqa_control),self$bobyqa_control$covrhoend,0)
                             #   rtsModel__set_cov_bobyqa_control(private$ptr,private$cov_type,private$lp_type,rhobeg,rhoend)
                             # }
                             
                             if(!is.null(lower_bound)){
                               rtsModel__set_bound(private$ptr,private$cov_type,private$lp_type,lower_bound,lower=TRUE)
                             }
                             if(!is.null(upper_bound)){
                               rtsModel__set_bound(private$ptr,private$cov_type,private$lp_type,upper_bound,lower=FALSE)
                             }
                             
                             # initialise the parameters and data on the R side
                             beta <- rtsModel__get_beta(private$ptr,private$cov_type,private$lp_type)
                             theta <- rtsModel__get_theta(private$ptr,private$cov_type,private$lp_type)
                             rho <- rtsModel__get_rho(private$ptr,private$cov_type,private$lp_type)
                             
                             xb <- rtsModel__xb(private$ptr,private$cov_type,private$lp_type)
                             all_pars <- c(beta = beta,theta = theta)
                             if(datlist$nT > 1) all_pars <- c(all_pars, rho = rho)
                             all_pars_new <- rep(1,length(all_pars))
                             beta_new <- beta
                             theta_new <- theta
                             rho_new <- rho
                             if(trace >= 1)cat("\nStart: ",all_pars,"\n")
                             t_start <- Sys.time()
                             L <- rtsModel__L(private$ptr,private$cov_type,private$lp_type)
                             ar_chol <- rtsModel__ar_chol(private$ptr,private$cov_type,private$lp_type)
                             data <- list(
                               N = datlist$Nsample,
                               Q = ifelse(approx=="hsgp", m * m, datlist$Nsample),
                               nT = datlist$nT,
                               Xb = xb,
                               ZL = as.matrix(L),
                               y = datlist$y,
                               rho = rho,
                               ar_chol = ar_chol
                             )
                             
                             ## set up the stan model
                             if(!is.null(self$region_data)){
                               filecmd <- "mcml_poisson_region_cmd.stan"
                               filers <- "mcml_poisson_region"
                               if(private$lp_type == 3){
                                 xb_grid <- rtsModel__region_grid_xb(private$ptr,private$cov_type)
                               } else {
                                 xb_grid <- matrix(0,nrow=datlist$Nsample*datlist$nT,ncol=1)
                               }
                               P <- rtsModel__grid_to_region_multiplier_matrix(private$ptr,private$cov_type,private$lp_type)
                               data$Xb = exp(data$Xb)
                               data <- append(data,list(nRegion = datlist$n_region,
                                                        ssize = length(P$Ai),
                                                        Ap = P$Ap+1,
                                                        Ai = P$Ai+1,
                                                        Ax = P$Ax))
                             } else {
                               filecmd <- "mcml_poisson_cmd.stan"
                               filers <- "mcml_poisson"
                             }
                             if(use_cmdstanr){
                               if(!requireNamespace("cmdstanr")){
                                 stop("cmdstanr is recommended but not installed, see https://mc-stan.org/cmdstanr/ for details on how to install.\n
                                    Set option use_cmdstanr=FALSE to use rstan instead.")
                               } else {
                                 if(trace == 2)message("If this is the first time running this model, it will be compiled by cmdstan.")
                                 model_file <- system.file("cmdstan",
                                                           filecmd,
                                                           package = "rts2",
                                                           mustWork = TRUE)
                                 mod <- suppressMessages(cmdstanr::cmdstan_model(model_file))
                               }
                             }
                             
                             # this is the main algorithm. iterate until convergence
                             iter <- 0
                             while(any(abs(all_pars-all_pars_new)>tol)&iter < max.iter){
                               
                               # step 1. MCMC sampling of random effects
                               all_pars <- all_pars_new
                               theta <- theta_new
                               if(data$nT > 1) rho <- rho_new
                               iter <- iter + 1
                               if(trace >= 1)cat("\nIter: ",iter,"\n",Reduce(paste0,rep("-",40)))
                               if(trace==2)t1 <- Sys.time()
                               
                               # update the data for stan
                               data$Xb <- rtsModel__xb(private$ptr,private$cov_type,private$lp_type)
                               data$ZL <- rtsModel__L(private$ptr,private$cov_type,private$lp_type)
                               data$rho <- rho
                               if(!is.null(self$region_data)){
                                 P <- rtsModel__grid_to_region_multiplier_matrix(private$ptr,private$cov_type,private$lp_type)
                                 data$Ax <- P$Ax
                                 data$Xb <- exp(data$Xb)
                               }
                               if(data$nT > 1){
                                 data$ar_chol <- rtsModel__ar_chol(private$ptr,private$cov_type,private$lp_type)
                               }
                               if(use_cmdstanr){
                                 if(is.null(self$region_data))data$Xb <-  rtsModel__xb(private$ptr,private$cov_type,private$lp_type)
                                 if(trace >= 1)cat("\nStarting MCMC sampling")
                                 capture.output(fit <- mod$sample(data = data,
                                                                  chains = 1,
                                                                  iter_warmup = iter_warmup,
                                                                  iter_sampling = iter_sampling,
                                                                  refresh = 0),
                                                file=tempfile())
                                 dsamps <- fit$draws("gamma",format = "matrix")
                                 class(dsamps) <- "matrix"
                                 rtsModel__update_u(private$ptr,as.matrix(t(dsamps)),private$cov_type,private$lp_type)
                               } else {
                                 capture.output(suppressWarnings( res <- rstan::sampling(stanmodels[[filers]],
                                                                                         data=data,
                                                                                         chains=1,
                                                                                         iter = iter_warmup+iter_sampling,
                                                                                         warmup = iter_warmup,
                                                                                         cores = 1,
                                                                                         refresh = 0)), file=tempfile())
                                 dsamps <- rstan::extract(res,"gamma",permuted = FALSE)
                                 dsamps <- dsamps[,1,]
                                 rtsModel__update_u(private$ptr,as.matrix(t(dsamps)),private$cov_type,private$lp_type)
                               }
                               if(trace==2)t2 <- Sys.time()
                               td1 <- t2-t1
                               if(trace==2)cat("\nMCMC sampling took: ",td1[[1]],attr(td1,"units"))
                               
                               # step 2. fit beta 
                               if(algo %in% c(1,3) & private$lp_type == 1){
                                 tryCatch(rtsModel__ml_beta(private$ptr,2,private$cov_type,private$lp_type),
                                          error = function(e){
                                            message(conditionMessage(e))
                                            cat("\nL-BFGS failed beta: trying BOBYQA")
                                            rtsModel__ml_beta(private$ptr,0,private$cov_type,private$lp_type)
                                          })
                               } else {
                                 rtsModel__ml_beta(private$ptr,0,private$cov_type,private$lp_type)
                               }
                               
                               if(is.null(known_theta)){
                                 if(algo == 1){
                                   tryCatch(rtsModel__ml_theta(private$ptr,2,private$cov_type,private$lp_type),
                                            error = function(e){
                                              message(conditionMessage(e))
                                              cat("\nL-BFGS failed theta: trying BOBYQA")
                                              rtsModel__ml_theta(private$ptr,0,private$cov_type,private$lp_type)
                                            })
                                 } else {
                                   rtsModel__ml_theta(private$ptr,0,private$cov_type,private$lp_type)
                                 }
                                if(datlist$nT > 1){
                                  if(algo == 1){
                                    tryCatch(rtsModel__ml_rho(private$ptr,2,private$cov_type,private$lp_type),
                                             error = function(e){
                                               message(conditionMessage(e))
                                               cat("\nL-BFGS failed rho, trying BOBYQA")
                                               rtsModel__ml_rho(private$ptr,0,private$cov_type,private$lp_type)
                                             })
                                  } else {
                                    rtsModel__ml_rho(private$ptr,0,private$cov_type,private$lp_type)
                                  }
                                }
                               }
                               beta_new <- rtsModel__get_beta(private$ptr,private$cov_type,private$lp_type)
                               
                               # step 3 fit covariance parameters
                               theta_new <- rtsModel__get_theta(private$ptr,private$cov_type,private$lp_type)
                               all_pars_new <- c(beta_new,theta_new)
                               if(datlist$nT > 1){
                                 rho_new <- rtsModel__get_rho(private$ptr,private$cov_type,private$lp_type)
                                 all_pars_new <- c(all_pars_new,rho_new)
                               }
                               
                               # some summary output
                               if(trace==2)t3 <- Sys.time()
                               td1 <- t3-t2
                               if(trace==2)cat("\nModel fitting took: ",td1[[1]],attr(td1,"units"))
                               if(trace==2){
                                 cat("\nPARAMETER VALUES:\n")
                                 cat("\nBeta: ", beta_new)
                                 cat("\nTheta: ", theta_new)
                                 if(datlist$nT > 1) cat(" | Rho: ", rho_new)
                                 cat("\nDifferences:")
                                 cat("\nBeta diff: ", round(max(abs(all_pars-all_pars_new)[1:length(beta)]),5))
                                 cat("\nTheta diff: ", round(max(abs(theta-theta_new)),5))
                                 if(datlist$nT > 1) cat(" | Rho diff: ", round(abs(rho-rho_new),5))
                                 cat("\nMax. diff: ", round(max(abs(all_pars-all_pars_new)),5))
                                 cat("\n",Reduce(paste0,rep("-",40)))
                               }
                             }
                             
                             # end of algorithm. 
                             not_conv <- iter > max.iter|any(abs(all_pars-all_pars_new)>tol)
                             if(not_conv)message(paste0("algorithm not converged. Max. difference between iterations :",round(max(abs(all_pars-all_pars_new)),4)))
                             
                             ## get the standard errors
                             if(trace >= 1)cat("\n\nCalculating standard errors...\n")
                             u <- rtsModel__u(private$ptr, private$cov_type,private$lp_type)
                             if(private$lp_type == 1){
                               M <- rtsModel__information_matrix(private$ptr,private$cov_type,private$lp_type)
                             } else {
                               M <- rtsModel__information_matrix_region(private$ptr,private$cov_type,private$lp_type)
                               # if(private$cov_type == 1 | private$cov_type == 2){
                               #   M <- rtsModel__information_matrix_region(private$ptr,private$cov_type,private$lp_type)
                               # } else {
                               #   M <- rtsModel__hessian_numerical(private$ptr,1e-4,private$cov_type,private$lp_type)
                               # }
                               M <- tryCatch(solve(M),error = function(i)return(diag(nrow(M))))
                             }
                             SE <- sqrt(diag(M))[1:length(beta)]
                             
                             # prepare output
                             beta_names <- rtsModel__beta_parameter_names(private$ptr,private$cov_type,private$lp_type)
                             theta_names <- c("theta_1","theta_2")
                             rho_names <- "rho"
                             mf_pars_names <- c(beta_names, theta_names)
                             SE <- c(SE,NA,NA)
                             if(datlist$nT > 1)SE <- c(SE,NA)
                             SE <- c(SE,rep(NA,nrow(u)))
                             if(datlist$nT > 1) mf_pars_names <- c(mf_pars_names, rho_names)
                             res <- data.frame(par = c(mf_pars_names,paste0("d",1:nrow(u))),
                                               est = c(all_pars_new,rowMeans(u)),
                                               SE=SE,
                                               t = NA,
                                               p = NA,
                                               lower = NA,
                                               upper = NA)
                             res$t <- res$est/res$SE
                             res$p <- 2*(1-stats::pnorm(abs(res$t)))
                             res$lower <- res$est - qnorm(1-0.05/2)*res$SE
                             res$upper <- res$est + qnorm(1-0.05/2)*res$SE
                             aic <- rtsModel__aic(private$ptr,private$cov_type,private$lp_type)
                             xb <- rtsModel__xb(private$ptr,private$cov_type,private$lp_type)
                             zd <- rowMeans(u)
                             wdiag <- Matrix::diag(rtsModel__get_W(private$ptr,private$cov_type,private$lp_type))
                             total_var <- var(Matrix::drop(xb)) + var(Matrix::drop(zd)) + mean(wdiag)
                             condR2 <- (var(Matrix::drop(xb)) + var(Matrix::drop(zd)))/total_var
                             margR2 <- var(Matrix::drop(xb))/total_var
                             # now get predictions
                             if(private$lp_type == 1){
                               ypred <- u
                               for(i in 1:ncol(ypred)){
                                 ypred[,i] <- exp(ypred[,i] + xb)
                               }
                             } else {
                               ypred <- rtsModel__y_pred(private$ptr,private$cov_type,private$lp_type)
                             }
                             t_end <- Sys.time()
                             t_diff <- t_end - t_start
                             if(trace == 2)cat("Total time: ", t_diff[[1]], " ", attr(t_diff,"units"))
                             out <- list(coefficients = res,
                                         converged = !not_conv,
                                         method = "mcem",
                                         m = dim(u)[2],
                                         tol = tol,
                                         sim_lik = FALSE,
                                         aic = aic,
                                         se="gls",
                                         Rsq = c(cond = condR2,marg=margR2),
                                         logl = rtsModel__log_likelihood(private$ptr,private$cov_type,private$lp_type),
                                         mean_form = "",
                                         cov_form = "",
                                         family = "poisson",
                                         link = "log",
                                         re.samps = u,
                                         iter = iter,
                                         time = t_diff,
                                         dof = length(xb),
                                         P = length(beta_new),
                                         Q = 2,
                                         var_par_family = FALSE,
                                         y=datlist$y,
                                         y_predicted  = ypred)
                             class(out) <- "mcmlrts"
                             return(out)
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
                                  is(fit,"CmdStanVB")|is(fit,"mcmlrts")))stop("stan fit or MCMCML fit required")
                             if("pred"%in%type&is.null(popdens))stop("set popdens for pred")

                             nCells <- nrow(self$grid_data)
                             nRegion <- ifelse(is.null(self$region_data),0,nrow(self$region_data))
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
                             } else if(is(fit,"mcmlrts")){
                               ypred <- t(fit$y_predicted)
                               f <- t(fit$re.samps)
                               nT <- ncol(f)/nCells
                               cmdst <- FALSE
                             }

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
                                 if(is.null(self$region_data)){
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
                                 } else {
                                   if(!cmdst){
                                     self$region_data$irr <- apply(ypred[,((nT-1-t.lag)*nRegion+1):((nT-t.lag)*nRegion),drop=FALSE]/
                                                                   ypred[,((nT-irr.lag-t.lag)*nRegion+1):(((nT-t.lag)-irr.lag+1)*nRegion),drop=FALSE],2,mean)
                                     self$region_data$irr_sd <- apply(ypred[,((nT-1-t.lag)*nRegion+1):((nT-t.lag)*nRegion),drop=FALSE]/
                                                                      ypred[,((nT-irr.lag-t.lag)*nRegion+1):(((nT-t.lag)-irr.lag+1)*nRegion),drop=FALSE],2,sd)
                                   } else {
                                     self$region_data$irr <- apply(ypred[,,((nT-1-t.lag)*nRegion+1):((nT-t.lag)*nRegion),drop=FALSE]/
                                                                   ypred[,,((nT-irr.lag-t.lag)*nRegion+1):(((nT-t.lag)-irr.lag+1)*nRegion),drop=FALSE],3,mean)
                                     self$region_data$irr_sd <- apply(ypred[,,((nT-1-t.lag)*nRegion+1):((nT-t.lag)*nRegion),drop=FALSE]/
                                                                      ypred[,,((nT-irr.lag-t.lag)*nRegion+1):(((nT-t.lag)-irr.lag+1)*nRegion),drop=FALSE],3,sd)
                                   }
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
                           #' @param fit A \link[rstan]{stanfit}, \link[cmdstanr]{CmdStanMCMC}, or `mcmlrts` object.
                           #' Output of `lgcp_bayes()` or `lgcp_ml()`
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
                           #' @return None, called for effects. Columns are added to grid or region data.
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
                           hotspots = function(fit,
                                               incidence.threshold=NULL,
                                               irr.threshold=NULL,
                                               irr.lag=NULL,
                                               rr.threshold=NULL,
                                               popdens,
                                               col_label=NULL){

                             if(all(is.null(incidence.threshold),is.null(irr.threshold),is.null(rr.threshold)))stop("At least one criterion required.")
                             if(!is.null(irr.threshold)&is.null(irr.lag))stop("irr.lag must be set")
                             if(!(is(fit,"CmdStanMCMC")|is(fit,"CmdStanVB")|is(fit,"stanfit")|is(fit,"mcmlrts")))stop("model fit required")
                             if(!is.null(self$region_data) & !is.null(rr.threshold) & (!is.null(irr.threshold) | !is.null(incidence.threshold)))stop("Cannot combine region-level measures (IRR/incidence) with grid relative risk.")
                             
                             useRegion <- FALSE
                             if(!is.null(self$region_data) & is.null(rr.threshold))useRegion <- TRUE
                             
                             nCells <- nrow(self$grid_data)
                             nRegion <- ifelse(is.null(self$region_data),0,nrow(self$region_data))
                             if(is(fit,"stanfit")){
                               ypred <- rstan::extract(fit,"y_grid_predict")$y_grid_predict
                               f <- rstan::extract(fit,"f")$f
                               nT <- dim(ypred)[2]/nCells
                               f <- f[,((nT-1)*nCells+1):(nT*nCells),drop=FALSE]
                             } else if(is(fit,"CmdStanMCMC")){
                               if(requireNamespace("cmdstanr")){
                                 ypred <- fit$draws("y_grid_predict")
                                 ypred <- matrix(ypred, prod(dim(ypred)[1:2]), dim(ypred)[3])
                                 nT <- dim(ypred)[3]/nCells
                                 f <- fit$draws("f")
                                 f <- f[,,((nT-1)*nCells+1):(nT*nCells),drop=FALSE]
                                 f <- matrix(f, prod(dim(f)[1:2]), dim(f)[3])
                               }
                             } else if(is(fit,"CmdStanVB")){
                               if(requireNamespace("cmdstanr")){
                                 ypred <- fit$draws("y_grid_predict")
                                 nT <- dim(ypred)[3]/nCells
                                 f <- fit$draws("f")
                                 f <- f[,((nT-1)*nCells+1):(nT*nCells),drop=FALSE]
                               }
                             } else if(is(fit,"mcmlrts")){
                               ypred <- t(fit$y_predicted)
                               f <- t(fit$re.samps)
                               nT <- ncol(f)/nCells
                               f <- f[,((nT-1)*nCells+1):(nT*nCells),drop=FALSE]
                             }

                             nCr <- sum(c(!is.null(incidence.threshold),
                                          !is.null(irr.threshold),
                                          !is.null(rr.threshold)))
                             inc1 <- matrix(0,nrow=nrow(f),ncol=ifelse(useRegion,nRegion,nCells))

                             if(!is.null(incidence.threshold)){
                               if(!useRegion){
                                 fmu <- matrix(0,nrow=nrow(f),ncol=nCells)
                                 for(i in 1:nCells){
                                   fmu[,i] <- ypred[,((nT-1)*nCells+i)]/as.data.frame(self$grid_data)[i,popdens]
                                 }
                                 inc1 <- inc1 + I(fmu > incidence.threshold)*1
                               } else {
                                 fmu <- matrix(0,nrow=nrow(f),ncol=nRegion)
                                 for(i in 1:nRegion){
                                   fmu[,i] <- ypred[,((nT-1)*nRegion+i)]/as.data.frame(self$region_data)[i,popdens]
                                 }
                                 inc1 <- inc1 + I(fmu > incidence.threshold)*1
                               }
                             }

                             if(!is.null(irr.threshold)){
                               if(!useRegion){
                                 if(nT==1)stop("cannot estimate irr as only one time period") else {
                                   inc1 <- inc1 + I(ypred[,((nT-1)*nCells+1):(nT*nCells),drop=FALSE]/
                                                      ypred[,((nT-irr.lag)*nCells+1):((nT-irr.lag+1)*nCells),drop=FALSE] > irr.threshold)*1
                                 }
                               } else {
                                 if(nT==1)stop("cannot estimate irr as only one time period") else {
                                   inc1 <- inc1 + I(ypred[,((nT-1)*nRegion+1):(nT*nRegion),drop=FALSE]/
                                                      ypred[,((nT-irr.lag)*nRegion+1):((nT-irr.lag+1)*nRegion),drop=FALSE] > irr.threshold)*1
                                 }
                               }
                             }

                             if(!is.null(rr.threshold)){
                               inc1 <- inc1 + I(exp(f) > rr.threshold)*1
                             }

                             inc1 <- I(inc1 == nCr)*1
                             if(useRegion){
                               self$region_data$hotspot_prob <- apply(inc1,2,mean)
                               if(!is.null(col_label)){
                                 colnames(self$region_data)[length(colnames(self$region_data))] <- col_label
                               }
                             } else {
                               self$grid_data$hotspot_prob <- apply(inc1,2,mean)
                               if(!is.null(col_label)){
                                 colnames(self$grid_data)[length(colnames(self$grid_data))] <- col_label
                               }
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
                           #' @param yvar String naming the outcome variable to calculate the variogram for. Optional, if
                           #' not provided then the outcome count data will be used.
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
                                   yvar <- paste0("t",1:nT)
                                 }
                               } else {
                                 if(!yvar%in%colnames(self$grid_data))stop("yvar not in grid data")
                               }
                               
                               df <- suppressWarnings(as.data.frame(sf::st_coordinates(sf::st_centroid(self$grid_data))))
                               dfs <- semivariogram(as.matrix(df),
                                                    offs = offs,
                                                    y = stack(as.data.frame(self$grid_data)[,yvar])$values,
                                                    nbins,
                                                    nT = nT)
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
                                 offs <- as.data.frame(self$region_data)[,popdens]
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
                                   yvar <- paste0("t",1:nT)
                                 }
                               } else {
                                 if(!yvar%in%colnames(self$region_data))stop("yvar not in region data")
                               }
                               
                               message("Using centroids of regions as locations")
                               df <- suppressWarnings(as.data.frame(sf::st_coordinates(sf::st_centroid(self$region_data))))
                               dfs <- semivariogram(as.matrix(df),
                                                    offs = offs,
                                                    y = stack(as.data.frame(self$region_data)[,yvar])$values,
                                                    nbins,
                                                    nT)
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
                           #' the grid cells. This function reorders the grid cells. If this is a region data model,
                           #' then the intersections are recomputed.
                           #' @param option Either "y" for order of the y coordinate, "x" for order of the x coordinate,
                           #' "minimax"  in which the next observation in the order is the one which maximises the
                           #'  minimum distance to the previous observations,g1$grid_data <- g1$grid_data[o0,] or "random" which randomly orders them.
                           #'  @return No return, used for effects.
                           reorder = function(option="y", verbose = TRUE){
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
                               xm <- min(df$x) + diff(range(df$x))
                               ym <- min(df$y) + diff(range(df$y))
                               o1[which.min((df$x-xm)^2 + (df$y - ym)^2)] <- 1
                               for(i in 2:nrow(df)){
                                 dists <- matrix(0,nrow=nrow(df)-i+1, ncol=(i-1))
                                 for(j in 1:(i-1)){
                                   dists[,j] <- (df$x[is.na(o1)] - df$x[which(o1 == j)])^2 +
                                                       (df$y[is.na(o1)] - df$y[which(o1 == j)])^2
                                 }
                                 mindists <- apply(dists,1,min)
                                 o1[is.na(o1)][which.max(mindists)] <- i
                                 if(verbose)cat("\r",progress_bar(i,nrow(df)))
                               }
                               self$grid_data <- self$grid_data[o1,]
                             } else if(option=="random"){
                               o2 <- sample(1:nrow(df),nrow(df),replace=FALSE)
                               self$grid_data <- self$grid_data[o2,]
                             }
                             self$grid_data$grid_id <- 1:nrow(self$grid_data)
                             if(!is.null(self$region_data)){
                               tmp <- suppressWarnings(sf::st_intersection(self$grid_data[,"grid_id"],self$region_data[,"region_id"]))
                               n_Q <- nrow(tmp)
                               rID <- self$region_data$region_id
                               tmp$area <- as.numeric(sf::st_area(tmp))
                               a1 <-rep(aggregate(tmp$area,list(tmp$region_id),sum)$x,unname(table(tmp$region_id)))
                               tmp$w <- tmp$area/a1
                               private$intersection_data <- tmp
                             }
                           },
                           #' @description 
                           #' A list of prepared data
                           #' 
                           #' The class prepares data for use in the in-built estimation functions. The same data could be used 
                           #' for alternative models. This is a utility function to facilitate model fitting for custom models.
                           #' @param m The number of nearest neighbours or basis functions. 
                           #' @param popdens String naming the variable in the data specifying the offset. If not 
                           #' provided then no offset is used.
                           #' @param covs An optional vector of covariate names. For regional data models, this is specifically for the region-level covariates.
                           #' @param covs_grid An optional vector of covariate names for region data models, identifying the covariates at the grid level.
                           #' @return A named list of data items used in model fitting
                           data = function(m,
                                           approx,
                                           popdens,
                                           covs,
                                           covs_grid){
                             return(private$prepare_data(m,
                                                         model = "exp",
                                                         approx,
                                                         popdens,
                                                         covs,
                                                         covs_grid,
                                                         FALSE,
                                                         FALSE))
                           },
                           #' @description 
                           #' Returns the random effects stored in the object (if any) after using MCMCML fitting. For example, 
                           #' if a fitting procedure is stopped, the random effects can still be returned.
                           #' @return A matrix of random effects samples if a MCMCML model has been initialised, otherwise returns FALSE
                           get_random_effects = function(){
                             if(!is.null(private$ptr)){
                               u <- rtsModel__u(private$ptr,data$y,private$cov_type,private$lp_type)
                               return(u)
                             } else {
                               return(FALSE)
                             }
                           }
                         ),
                    private = list(
                      intersection_data = NULL,
                      ptr = NULL,
                      grid_ptr = NULL,
                      region_ptr = NULL,
                      cov_type = 1,
                      lp_type = 1,
                      update_ptr = function(m,
                                            model,
                                            approx,
                                            popdens,
                                            covs,
                                            covs_grid,
                                            L = 1.5,
                                            update = TRUE,
                                            formula_1 = NULL,
                                            formula_2 = NULL){

                        data <- private$prepare_data(m,
                                                     model,
                                                     approx,
                                                     popdens,
                                                     covs,
                                                     covs_grid,
                                                     FALSE,
                                                     FALSE,
                                                     L)
                        
                        if(is.null(private$grid_ptr)){
                          private$grid_ptr <- GridData__new(as.matrix(data$x_grid),data$nT)
                        }

                        if(!is.null(self$region_data) & (is.null(private$region_ptr) | update)){
                          rdata <- self$get_region_data()
                          private$region_ptr <- RegionData__new(rdata$n_cell-1,rdata$cell_id-1,rdata$q_weights, nrow(data$x_grid), data$nT)
                        }

                        if(is.null(private$ptr) || update){
                          # build formulae
                          if(!is.null(formula_1)){
                            if(!is(formula_1,"formula"))stop("Not a formula")
                            f1 <- as.character(formula_1[[2]])
                            f1 <- gsub(" ","",f1)
                          } else {
                            f1 <- "1"
                            if(length(covs) > 0){
                              for(i in 1:length(covs)){
                                f1 <- paste0(f1,"+",covs[[i]])
                              }
                            }
                          }
                          

                          if(length(covs_grid)>0 || !is.null(formula_2)){
                            if(!is.null(formula_2)){
                              if(!is(formula_2,"formula"))stop("Not a formula")
                              f2 <- as.character(formula_2[[2]])
                              f2 <- gsub(" ","",f2)
                            } else {
                              f2 <- "1"
                              for(i in 1:length(covs)){
                                f2 <- paste0(f2,"+",covs_grid[[i]])
                              }
                            }
                          } else {
                            f2 <- "1"
                          }
                          
                          # add random effects structure
                          if(model=="exp"){
                            f1 <- paste0(f1,"+(1|fexp(X,Y))")
                            if(length(covs_grid)>0)f2 <- paste0(f2,"+(1|fexp(X,Y))")
                          } else if(model=="sqexp"){
                            f1 <- paste0(f1,"+(1|sqexp(X,Y))")
                            if(length(covs_grid)>0)f2 <- paste0(f2,"+(1|sqexp(X,Y))")
                          } else {
                            stop("Only exp and spexp for now.")
                          }

                          if(approx == "nngp"){
                            private$cov_type <- 2
                          } else if(approx == "lgcp" || approx == "none"){
                            private$cov_type <- 1
                          } else {
                            private$cov_type <- 3
                          }
                          if(!is.null(self$region_data)){
                            if(length(covs_grid) > 0){
                              private$lp_type <- 3
                            } else {
                              private$lp_type <- 2
                            }
                          } else {
                            private$lp_type <- 1
                          }
                          
                          # use random starting values with sensible intercept
                          P <- length(covs) + length(covs_grid)
                          beta <- c(mean(log(mean(data$y)) - log(data$popdens)))
                          if(P>0)beta <- c(beta, rnorm(P,0,0.1))
                          # set sensible starting value for theta 1
                          theta1 <- abs((log(max(data$y)) - beta[1])/2)/2
                          theta <- c(theta1,runif(1,0.1,0.5))
                          if(private$cov_type == 2 & private$lp_type == 1){
                            private$ptr <- Model_nngp_lp__new(f1,
                                                              as.matrix(data$X),
                                                              as.matrix(data$x_grid),
                                                              c("intercept",covs),
                                                              beta,
                                                              theta,
                                                              data$nT,
                                                              m,
                                                              private$grid_ptr)
                          } else if(private$cov_type == 1 & private$lp_type == 1){
                            private$ptr <- Model_ar_lp__new(f1,
                                                            as.matrix(data$X),
                                                            as.matrix(data$x_grid),
                                                            c("intercept",covs),
                                                            beta,
                                                            theta,
                                                            data$nT)
                          } else if(private$cov_type == 3 & private$lp_type == 1){
                            private$ptr <- Model_hsgp_lp__new(f1,
                                                              as.matrix(data$X),
                                                              as.matrix(data$x_grid),
                                                            c("intercept",covs),
                                                            beta,
                                                            theta,
                                                            data$nT,
                                                            m,
                                                            data$L)
                          } else if(private$cov_type == 1 & private$lp_type == 2){
                            private$ptr <- Model_ar_region__new(f1,
                                                                as.matrix(data$X),
                                                                as.matrix(data$x_grid),
                                                                c("intercept",covs),
                                                                beta,
                                                                theta,
                                                                data$nT,
                                                                private$region_ptr)
                          } else if(private$cov_type == 2 & private$lp_type == 2){
                            private$ptr <- Model_nngp_region__new(f1,
                                                              as.matrix(data$X),
                                                              as.matrix(data$x_grid),
                                                              c("intercept",covs),
                                                              beta,
                                                              theta,
                                                              data$nT,
                                                              m,
                                                              private$region_ptr,
                                                              private$grid_ptr)
                          } else if(private$cov_type == 3 & private$lp_type == 2){
                            private$ptr <- Model_hsgp_region__new(f1,
                                                                  as.matrix(data$X),
                                                                  as.matrix(data$x_grid),
                                                                  c("intercept",covs),
                                                                  beta,
                                                                  theta,
                                                                  data$nT,
                                                                  m,
                                                                  private$region_ptr,
                                                                  data$L)
                          } else if(private$cov_type == 1 & private$lp_type == 3){

                            private$ptr <- Model_ar_region_grid__new(f1,
                                                                f2,
                                                                as.matrix(data$X),
                                                                as.matrix(data$x_grid),
                                                                c("intercept",covs),
                                                                colnames(data$x_grid),
                                                                beta,
                                                                theta,
                                                                private$region_ptr,
                                                                data$nT)
                          } else if(private$cov_type == 2 & private$lp_type == 3){
                            private$ptr <- Model_nngp_region_grid__new(f1,
                                                              f2,
                                                              as.matrix(data$X),
                                                              as.matrix(data$x_grid),
                                                              c("intercept",covs),
                                                              colnames(data$x_grid),
                                                              beta,
                                                              theta,
                                                              private$region_ptr,
                                                              private$grid_ptr,
                                                              data$nT,
                                                              m)
                          } else if(private$cov_type == 3 & private$lp_type == 3){
                            private$ptr <- Model_hsgp_region_grid__new(f1,
                                                                       f2,
                                                                       as.matrix(data$X),
                                                                       as.matrix(data$x_grid),
                                                                       c("intercept",covs),
                                                                       colnames(data$x_grid),
                                                                       beta,
                                                                       theta,
                                                                       private$region_ptr,
                                                                       data$nT,
                                                                       m,
                                                                       data$L)
                          }

                          rtsModel__set_offset(private$ptr,log(data$popdens),private$cov_type,private$lp_type)
                          rtsModel__set_y(private$ptr,data$y,private$cov_type,private$lp_type)
                          return(data)
                        }
                      },
                      prepare_data = function(m,
                                              model,
                                              approx,
                                              popdens,
                                              covs,
                                              covs_grid,
                                              verbose,
                                              bayes,
                                              L = 1.5){
                        
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
                          diffs <- c((xrange[1]+xrange[2])/2, (yrange[1]+yrange[2])/2)
                          # x_grid[,1] <- (x_grid[,1]- mean(xrange))/std_val
                          # x_grid[,2] <- (x_grid[,2]- mean(yrange))/std_val
                          x_grid[,1] <- (x_grid[,1]- diffs[1])
                          x_grid[,2] <- (x_grid[,2]- diffs[2])
                          L_boundary <- c(L*(xrange[2]-xrange[1])/2, L*(yrange[2]-yrange[1])/2)
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
                            }
                            Q <- nQ+1
                          } else {
                            X <- matrix(1,nrow=length(y),ncol=1)
                            Q <- 1
                          }
                          if(length(covs_grid)>0){
                            nG <- length(covs_grid)
                            X_g <- matrix(NA,nrow=nrow(self$grid_data)*nT,ncol=nG)
                            for(i in 1:nG){
                              nColV <- sum(grepl(covs_grid[i],colnames(self$grid_data)))
                              if(nColV==1){
                                X_g[,i] <- rep(as.data.frame(self$grid_data)[,covs[i]],nT)
                              } else if(nColV==0){
                                stop(paste0(covs_grid[i]," not found in grid data"))
                              } else {
                                if(nT>1){
                                  X_g[,i] <- stack(as.data.frame(self$grid_data)[,paste0(covs_grid[i],1:nT)])[,1]
                                } else {
                                  X_g[,i] <- as.data.frame(self$grid_data)[,covs_grid[i]]
                                }
                              }
                            }
                          } else {
                            nG <- 0
                            X_g <- matrix(0,nrow=nrow(self$grid_data)*nT,ncol=1)
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
                            L = L_boundary,
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
                              Q_g = nG,
                              X_g = X_g,
                              n_region = nrow(self$region_data),
                              n_Q = nrow(private$intersection_data),
                              n_cell = ncell,
                              cell_id = private$intersection_data$grid_id,
                              q_weights = private$intersection_data$w
                            ))
                            if(length(covs_grid)>0)datlist$x_grid <- cbind(datlist$x_grid,as.data.frame(self$grid_data)[,covs_grid])
                          } 
                          
                          
                        } else if(approx == "nngp"){
                          if(is.null(private$grid_ptr)){
                            private$grid_ptr <- GridData__new(as.matrix(x_grid[,1:2]),nT)
                          }
                          GridData__gen_NN(private$grid_ptr,m)
                          NN <- GridData__NN(private$grid_ptr)
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
                              Q_g = nG,
                              X_g = X_g,
                              n_region = nrow(self$region_data),
                              n_Q = nrow(private$intersection_data),
                              n_cell = ncell,
                              cell_id = private$intersection_data$grid_id,
                              q_weights = private$intersection_data$w
                            ))
                            if(length(covs_grid)>0)datlist$x_grid <- cbind(datlist$x_grid,as.data.frame(self$grid_data)[,covs_grid])
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
                              Q_g = nG,
                              X_g = X_g,
                              n_region = nrow(self$region_data),
                              n_Q = nrow(private$intersection_data),
                              n_cell = ncell,
                              cell_id = private$intersection_data$grid_id,
                              q_weights = private$intersection_data$w
                            ))
                            if(length(covs_grid)>0)datlist$x_grid <- cbind(datlist$x_grid,as.data.frame(self$grid_data)[,covs_grid])
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

