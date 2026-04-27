#' An rts grid object
#' 
#' An rts grid object is an R6 class holding the spatial data with data, model fitting, and analysis functions. 
#' @details
#' **INTRODUCTION**
#' 
#' The various methods of the class include examples and details of their implementation. The \pkg{sf} package is used for 
#' all spatial data. A typical workflow with this class would 
#' be:
#' \enumerate{
#'  \item Create a new grid object. The class is initialized with either a single polygon describing the area of interest or a collection 
#'  of polygons if spatially aggregated data are used. 
#'  \item If the location (and times) of cases are available (i.e. the data are not spatially aggregated), then we map the points to the computational 
#'  grid. The function \link[rts2]{create_points} can generate point data in the correct \pkg{sf} format. The member function `points_to_grid` will then 
#'  map these data to the grid. Counts can also be manually added to grid data. For region data, since the counts are assumed to be already aggregated, these 
#'  must be manually provided by the user. The case counts must appear in columns with specific names. If there is only a single time period then the counts 
#'  must be in a column named `y`. If there are multiple time periods then the counts must be in columns names `t1`, `t2`, `t3`,... Associated columns labelled 
#'  `date1`, `date2`, etc. will permit use of some functionality regarding specific time intervals.
#'  \item If any covariates are to be used for the modelling, then these can be mapped to the compuational grid using the function `add_covariates()`. Other 
#'  functions, `add_time_indicators()` and `get_dow()` will also generate relevant temporal indicators where required. At a minimum we would recommend including 
#'  a measure of population density.
#'  \item Fit a model. There are multiple methods for model fitting, which are available through the member functions `lgcp_ml()` and `lgcp_bayes()` for maximum likelihood 
#'  and Bayesian approaches, respectively. The results are stored internally and optionally returned as a `rtsFit` object. 
#'  \item Summarise the output. The main functions for summarising the output are `extract_preds()`, which will generate predictions of relative risk, incidence rate 
#'  ratios, and predicted incidence, and `hotspots()`, which will estimate probabilities that these statistics exceed given thresholds. For spatially-aggregated data models, 
#'  the relative risk applies to the grid, whereas rate ratios and predicted incidence applies to the areas.
#'  \item Predictions can be visualised or aggregated to relevant geographies with the `plot()` and `aggregate()` functions.
#' }
#' Specific details of the implementation of each of these functions along with examples appear below. 
#' 
#' @importFrom R6 R6Class
#' @importFrom stats weighted.mean sd model.matrix
#' @importFrom utils stack capture.output
#' @export
grid <- R6::R6Class("grid",
                         public = list(
                           #' @field grid_data sf object specifying the computational grid for the analysis
                           grid_data = NULL,
                           #' @field point_data sf object with the event locations
                           point_data = NULL,
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
                           #' Create a new grid object
                           #'
                           #' Produces a regular grid over an area of interest as an sf object, see details for information on initialisation.
                           #'
                           #' @param poly An sf object containing either one polygon describing the area of interest or multiple polygons
                           #' representing survey or census regions in which the case data counts are aggregated
                           #' @param cellsize The dimension of the grid cells
                           #' @param verbose Logical indicating whether to provide feedback to the console.
                           #' @return NULL
                           #' @examples
                           #' # a simple example with a square and a small number of cells
                           #' # this same running example is used for the other functions 
                           #' b1 = sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
                           #' g1 <- grid$new(b1,0.5)
                           #' 
                           #' # an example with multiple polygons
                           #' data("birmingham_crime")
                           #' g2 <- grid$new(birmingham_crime,cellsize = 1000)
                           initialize = function(poly,
                                                 cellsize,
                                                 verbose = TRUE){
                             if(!is(poly,"sf"))stop("data not sf")
                             if(!is(cellsize,"numeric"))stop("cellsize not numeric")
                             if(nrow(poly)>1 & verbose)message("Multiple polygons in data. Assuming analysis uses counts aggregated to an irregular lattice and not point data.")
                             if(nrow(poly)==1){
                               bgrid <- sf::st_make_grid(poly,cellsize = cellsize)
                               bgrid <- sf::st_sf(bgrid)
                               idx1<- sf::st_contains(y=bgrid,x=poly)
                               idx2<- sf::st_intersects(y=bgrid,x=poly)
                               idx <- c(unlist( sapply( idx1, `[`) ),unlist( sapply( idx2, `[`) ))
                               bgrid <- bgrid[idx,]
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
                               #a1 <-rep(aggregate(tmp$area,list(tmp$region_id),sum)$x,unname(table(tmp$region_id)))
                               cellsize <- sf::st_area(self$grid_data[1,])
                               tmp$w <- tmp$area/cellsize
                               self$region_data <- poly
                               tmp <- tmp[order(tmp$region_id),]
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
                           #' Prints this object
                           #' @return None. called for effects.
                           print = function(){
                             is_setup <- FALSE
                             cat("\nAn rts2 grid object\n")
                             cat("\n Computational grid\n     \U2BA1 ",nrow(self$grid_data)," cells")
                             if(is.null(self$region_data)){
                               nT <- sum(grepl("\\bt[0-9]",colnames(self$grid_data)))
                               if(nT == 0) {
                                 nT <- 1
                                 if("y"%in%colnames(self$grid_data)){
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
                             } else {
                               cat("\n \U2BC8 Spatially aggregated count data:\n     \U2BA1 ",nrow(self$region_data)," regions")
                               nT <- sum(grepl("\\bt[0-9]",colnames(self$region_data)))
                               if(nT == 0){
                                 nT <- 1
                                 if("y"%in%colnames(self$region_data)){
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
                             }
                             cat("\n")
                           },
                           #' @description
                          #' Reproject coordinates of spatial data in the object to a new CRS
                          #' 
                          #' @param crs target coordinate reference system passed to `st_transform`
                          #' @param scale_factor A scale factor to multiply existing geometry by.
                          #' @return None, called for effects
                           transform = function(crs, scale_factor){
                             if(!missing(crs) & !missing(scale_factor))stop("set only one of crs or scale factor")
                             if(!missing(crs)){
                               if(!is.null(self$boundary))self$boundary <- sf::st_transform(self$boundary, crs = crs)
                               if(!is.null(self$point_data)) self$point_data <- sf::st_transform(self$point_data, crs = crs)
                               if(!is.null(self$grid_data)) self$grid_data <- sf::st_transform(self$grid_data, crs = crs)
                               if(!is.null(self$region_data)) self$region_data <- sf::st_transform(self$region_data, crs = crs)
                             } else {
                               if(!is.null(self$boundary))sf::st_geometry(self$boundary) <- sf::st_geometry(self$boundary) * scale_factor 
                               if(!is.null(self$point_data)) sf::st_geometry(self$point_data) <- sf::st_geometry(self$point_data) * scale_factor 
                               if(!is.null(self$grid_data)) sf::st_geometry(self$grid_data) <- sf::st_geometry(self$grid_data) * scale_factor 
                               if(!is.null(self$region_data)) sf::st_geometry(self$region_data) <- sf::st_geometry(self$region_data) * scale_factor 
                             }
                             
                           },
                           #' @description
                           #' Plots the grid data
                           #'
                           #' @details
                           #' **PLOTTING**
                           #' 
                           #' If `zcol` is not specified then only the geometry is plotted, otherwise the covariates specified will be plotted.
                           #' The user can also use sf plotting functions on self$grid_data and self$region_data directly.
                           #' @param zcol Vector of strings specifying names of columns of `grid_data` to plot
                           #' @return A plot
                           #' @examples
                           #' b1 = sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
                           #' g1 <- grid$new(b1,0.5)
                           #' g1$plot()
                           #' 
                           #' # a plot with covariates - we simulate covariates first
                           #' g1$grid_data$cov <- stats::rnorm(nrow(g1$grid_data))
                           #' g1$plot("cov")
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
                           #' @details
                           #' **POINTS TO GRID**
                           #' 
                           #' Given the sf object with the point locations and date output from
                           #' `create_points()`, the functions will add columns to `grid_data` indicating
                           #' the case count in each cell in each time period.
                           #' 
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
                           #' @param date_format String describing the format of the date in the data as a combination of "d" days, "m" months, 
                           #' and "y" years, either "dmy", "myd", "ymd", "ydm", "dym" "mdy" as used by the lubridate package.
                           #' @param verbose Logical indicating whether to report detailed output
                           #' @return NULL
                           #' @seealso \link[rts2]{create_points}
                           #' @examples
                           #' b1 <- sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
                           #' g1 <- grid$new(b1,0.5)
                           #' # simulate some points
                           #' dp <- data.frame(y=runif(10,0,3),x=runif(10,0,3),date=paste0("2021-01-",11:20)) 
                           #' dp <- create_points(dp,pos_vars = c('y','x'),t_var='date')
                           #' g1$points_to_grid(dp, laglength=5)
                           points_to_grid = function(point_data,
                                                     t_win = c("day"),
                                                     laglength = 14,
                                                     date_format = "ymd",
                                                     verbose = TRUE){

                             if(!is(point_data,"sf"))stop("points not sf")
                             if(sf::st_crs(point_data)!=sf::st_crs(self$grid_data)){
                               warning("CRS not equal. Setting st_crs(point_data)==st_crs(self$grid_data)")
                               sf::st_crs(point_data) <- sf::st_crs(self$grid_data)
                             }
                             self$point_data <- point_data
                             if("y"%in%colnames(self$grid_data))self$grid_data <- self$grid_data[,-which(colnames(self$grid_data)=="y")]
                             nT <- sum(grepl("\\bt[0-9]",colnames(self$grid_data)))
                             if(nT > 0){
                               tidx <- which(grepl("\\bt[0-9]",colnames(self$grid_data)))
                               self$grid_data <- self$grid_data[,-tidx]
                               tidx <- which(grepl("\\bdate[0-9]",colnames(self$grid_data)))
                               self$grid_data <- self$grid_data[,-tidx]
                             }
                             
                             if("t"%in%colnames(point_data)){
                               if(!t_win%in%c("day","week","month","year"))stop("t_win not day, week, month, or year")
                               if(!date_format%in%c("dmy","myd","ymd","ydm","dym","mdy"))stop("date format not recognised")
                               
                               tvals <- do.call(eval(parse(text = paste0("lubridate::",date_format))),list(point_data$t))
                               tuniq <- sort(tvals)
                               tdat <- lubridate::floor_date(tvals, t_win, week_start = 1)
                               tuniq <- lubridate::floor_date(tuniq, t_win, week_start = 1)
                               tuniq <- unique(tuniq)
                               difft <- c()
                               for(i in 1:(length(tuniq)-1))difft <- c(difft, lubridate::interval(tuniq[i+1],tuniq[i])%/% months(1))
                               if(verbose)message(paste("There are "),length(tuniq)," unique time periods in the data")
                               if(laglength > length(tuniq))stop("laglength exceeds unique time periods")
                               if(length(unique(difft))>1)message("Note that the data includes non-consecutive or differing length periods. The data will be extracted for the specified number of time periods, however the model will assume equal time differences.")
                               tuniq <- tuniq[(length(tuniq)-laglength+1):length(tuniq)]
                               
                               for(i in 1:length(tuniq))
                               {
                                 self$grid_data$y <-  lengths(sf::st_intersects(self$grid_data, point_data[tdat==tuniq[i],]))
                                 if(length(tuniq)>1)colnames(self$grid_data)[length(colnames(self$grid_data))] <- paste0("t",i)
                                 self$grid_data$d <- min(point_data[tdat==tuniq[i],]$t)
                                 colnames(self$grid_data)[length(colnames(self$grid_data))] <- paste0("date",i)
                               }
                             } else {
                               self$grid_data$y <-  lengths(sf::st_intersects(self$grid_data,point_data))
                             }
                             if(verbose)message("added points data to grid data")
                             return(invisible(self))
                           },
                           #' @description 
                           #' Returns a data frame with the grid data and coordinates
                           #' 
                           #' Returns a standard data frame with the grid data and coordinates, which may be useful to 
                           #' run models in another package.
                           model_data_frame = function(){
                             if(is.null(self$region_data)){
                               df <- as.data.frame(self$grid_data)[,2:ncol(self$grid_data)]
                               df <- cbind(df, st_coordinates(st_centroid(self$grid_data)))
                               return(df)
                             } else {
                               stop("Not yet implemented for region data models")
                             }
                           },
                           #' @description
                           #' Adds covariate data to the grid
                           #'
                           #' Maps spatial, temporal, or spatio-temporal covariate data onto the grid using different methods.
                           #'
                           #' @details
                           #' **ADDING COVARIATES**
                           #' 
                           #' *Spatially-varying data only* 
                           #' 
                           #' `cov_data` is an object describing covariate over the area of interest. 
                           #' sf, RasterLayer and SpatRaster objects are supported, with rasters converted internally to sf.
                           #' The mapping can use a spatially-smoothed method (pynchophylatic) or a variance or entropy minimising 
                           #' method or a simple "flat" mapping using only the weighted average of intersections between the grid and 
                           #' covariate polygons. For non-"flat" mapping,  `lambda_smooth` controls the degree of spatial smoothing - if set to zero then no spatial smoothing 
                           #' is used. The argument `lambda_e` adds a small amount to reduce numerical instability. One can also map strictly positive covariates (e.g. population density)
                           #' by setting the `is_positive` argument to true. In this case, `lambda_e` is used to add an entropy minimising criterion (instead of, or in addition to)
                           #' spatial smoothing crtierion. If population density information, then this can be accounted for in the smoothing by setting 
                           #' `weight_type` to `pop` and specifying the name of covariate to `popdens`, which should be on the grid.
                           #'
                           #' *Temporally-varying only data* 
                           #' 
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
                           #' 
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
                           #' @param cov_data sf, RasterLayer, SpatRaster object or a data.frame. See details.
                           #' @param zcols vector of character strings with the names of the columns of `cov_data`
                           #' to include
                           #' @param weight_type character string. Either "area" for area-weighted average or "pop"
                           #' for population-weighted average
                           #' @param popdens character string. The name of the column in `cov_data` with the
                           #' population density. Required if weight_type="pop"
                           #' @param verbose logical. Whether to provide a progress bar
                           #' @param t_label integer. If adding spatio-temporally varying data by time period,
                           #' this time label should be appended to the column name. See details.
                           #' @param flat Logical indicating if the disaggregation should be flat and just a weighted average over intersections. Cannot be strictly positive.
                           #' @param lambda_smooth weight on spatial smoothness, used if `flat` is FALSE
                           #' @param lambda_e small ridge for numerical stability (needed because L is singular), or for 
                           #' strictly positive covariates, the weight on entropy in the minimisation. 
                           #' @param is_positive Logical. Should the disaggregation be strictly positive.
                           #' @return NULL
                           #' @examples
                           #' b1 <-  sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
                           #' g1 <- grid$new(b1,0.5)
                           #' cov1 <- grid$new(b1,0.8)
                           #' cov1$grid_data$cov <- runif(nrow(cov1$grid_data))
                           #' g1$add_covariates(cov1$grid_data,
                           #'                   zcols="cov")
                           #' 
                           #' \donttest{
                           #' # mapping population data from some other polygons
                           #' data("boundary")
                           #' data("birmingham_crime")
                           #' g2 <- grid$new(boundary,cellsize=0.008)
                           #' msoa <- sf::st_transform(birmingham_crime,crs = 4326)
                           #' suppressWarnings(sf::st_crs(msoa) <- sf::st_crs(g2$grid_data)) # ensure crs matches
                           #' g2$add_covariates(msoa,
                           #'                   zcols="pop",
                           #'                   weight_type="area")
                           #'                   
                           #' # add a case count                  
                           #' g2$add_covariates(msoa,
                           #'                   zcols=c("t1"),
                           #'                   weight_type="area",
                           #'                   is_positive = TRUE,
                           #'                   lambda_smooth = 0,
                           #'                   lambda_e = 1e-6)
                           #'
                           #' g2$grid_data$t1 <- round(g2$grid_data$t1,0)
                           #' g2$plot("pop")
                           #' }
                           #' @importFrom spdep nb2mat
                           add_covariates = function(cov_data,
                                                     zcols,
                                                     weight_type="area",
                                                     popdens=NULL,
                                                     t_label=NULL,
                                                     flat = TRUE,
                                                     lambda_smooth = 1,
                                                     lambda_e = 1e-6,
                                                     is_positive = FALSE){
                              
                             if(flat & is_positive)message("is_positive is ignored as flat=TRUE")
                             if(!weight_type%in%c("area","pop"))stop("Type must be area or pop.")
                             if((weight_type=="pop"&is.null(popdens)))stop("Pop. dens. variable not found.")
                             if(weight_type=="pop"&!is.null(popdens)){
                               if(!popdens%in%colnames(cov_data))stop("Pop. dens. variable not found.")
                             }
                             
                             if(is(cov_data,"RasterLayer")){
                               if(any(!zcols%in%names(cov_data)))stop("variable names not in cov_data")
                               cnames <- names(cov_data)
                               fname <- tempfile(fileext = ".tif")
                               raster::writeRaster(cov_data,fname)
                               x <- stars::read_stars(fname)
                               cov_data <- sf::st_as_sf(x)
                               colnames(cov_data)[1:(ncol(cov_data)-1)] <- cnames
                             } else if(is(cov_data,"SpatRaster")){
                               if(any(!zcols%in%names(cov_data)))stop("variable names not in cov_data")
                               cnames <- names(cov_data)
                               fname <- tempfile(fileext = ".tif")
                               raster::writeRaster(raster::raster(cov_data),fname)
                               x <- stars::read_stars(fname)
                               cov_data <- sf::st_as_sf(x)
                               colnames(cov_data)[1:(ncol(cov_data)-1)] <- cnames
                             } else if(is(cov_data,"sf")){
                               if(any(!zcols%in%colnames(cov_data)))stop("variable names not in cov_data")
                             }
                             
                             for(j in 1:length(zcols)){
                               if(zcols[j] %in% colnames(self$grid_data)){
                                 self$grid_data <- self$grid_data[,-which(colnames(self$grid_data) == zcols[j])]
                               }
                             }
                             
                             if(is(cov_data,"sf")){
                               sf::st_agr(cov_data) = "constant"
                               sf::st_agr(self$grid_data) = "constant"
                               self$grid_data$grid_id <- 1:nrow(self$grid_data)
                               cov_data$region_id <- 1:nrow(cov_data)
                               tmp <- suppressWarnings(sf::st_intersection(self$grid_data[,"grid_id"],cov_data[,"region_id"]))
                               n_Q <- nrow(tmp)
                               tmp$area <- as.numeric(sf::st_area(tmp))
                               if(nrow(tmp)==0)stop("Covariate data does not overlap grid")
                               a1 <-rep(aggregate(tmp$area,list(tmp$grid_id),sum)$x,unname(table(tmp$grid_id)))
                               tmp$w <- tmp$area/a1
                               ncell <- unname(table(tmp$region_id))
                               ncell <- c(1, cumsum(ncell)+1)
                               W <- as.matrix(Matrix::sparseMatrix(j = tmp$grid_id, p = ncell-1, x = as.vector(tmp$w)))
                               if(flat){
                                 for(j in 1:length(zcols)){
                                   out <- flat_disaggregate(as.data.frame(cov_data[unique(tmp$region_id),zcols[j]])[,1], W)
                                   self$grid_data$x <- out
                                   if(is.null(t_label)){
                                     colnames(self$grid_data)[length(colnames(self$grid_data))] <- zcols[j]
                                   } else {
                                     colnames(self$grid_data)[length(colnames(self$grid_data))] <- paste0(zcols[j],t_label)
                                   }
                                 }
                               } else {
                                 if(is.null(private$grid_adj_matrix)){
                                   private$grid_adj_matrix <- spdep::nb2mat(spdep::poly2nb(self$grid_data), style = "B")
                                 }
                                 if(weight_type == "pop"){
                                   pop <- self$grid_data[,popdens]
                                 } else {
                                   pop <- rep(1, nrow(self$grid_data))
                                 }
                                 for(j in 1:length(zcols)){
                                   if(is_positive){
                                     out <- disaggregate_positive(W, as.data.frame(cov_data[unique(tmp$region_id),zcols[j]])[,1], pop, private$grid_adj_matrix, lambda_smooth, lambda_e)
                                     self$grid_data$x <- out$y_s
                                     cat("Mapped positive covariate",zcols[j],"Morans's I:",round(out$moran_I,2),"Residual:",out$residual,"\n")
                                   } else {
                                     out <- disaggregate_covariate(as.data.frame(cov_data[unique(tmp$region_id),zcols[j]])[,1], W, private$grid_adj_matrix, pop, lambda_smooth, lambda_e)
                                     self$grid_data$x <- out$x_s
                                     cat("Mapped covariate",zcols[j],"Morans's I:",round(out$moran_I,2),"Mean Error:",out$mean_error,"\n")
                                   }
                                   if(is.null(t_label)){
                                     colnames(self$grid_data)[length(colnames(self$grid_data))] <- zcols[j]
                                   } else {
                                     colnames(self$grid_data)[length(colnames(self$grid_data))] <- paste0(zcols[j],t_label)
                                   }
                                 }
                               }
                             } else if(is(cov_data,"data.frame")){
                               if(!"t"%in%colnames(cov_data))stop("not column named t in cov_data")
                               nT <- max(cov_data$t)
                               for(j in zcols){
                                 for(t in 1:nT){
                                   self$grid_data$x <- cov_data[cov_data$t==t,j]
                                   colnames(self$grid_data)[length(colnames(self$grid_data))] <- paste0(j,t)
                                 }
                               }
                             } else {
                               stop("Cov_data type not supported")
                             }
                             return(invisible(self))
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
                           #' dow <- g1$get_dow()
                           #' g1$add_covariates(dow,zcols = colnames(dow)[3:ncol(dow)])
                           get_dow = function(){
                             if(is.null(self$region_data)){
                               nT <- sum(grepl("\\bt[0-9]",colnames(self$grid_data)))
                             } else {
                               nT <- sum(grepl("\\bt[0-9]",colnames(self$region_data)))
                             }
                             if(nT <= 1)stop("No time periods to generate day of week indicators for")
                             dw <- data.frame(t=1:nT,day=NA)
                             if(is.null(self$region_data)){
                               for(i in 1:nT) dw$day[i] <- as.character(lubridate::wday(as.data.frame(self$grid_data)[1,paste0("date",i)],label = TRUE))
                             } else {
                               for(i in 1:nT) dw$day[i] <- as.character(lubridate::wday(as.data.frame(self$region_data)[1,paste0("date",i)],label = TRUE))
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
                             } else {
                               nT <- sum(grepl("\\bt[0-9]",colnames(self$region_data)))
                             }
                             
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
                             return(invisible(self))
                           },
                           #' @description
                           #' Fit an (approximate) log-Gaussian Cox Process model using Bayesian methods
                           #'
                           #' @details
                           #' **BAYESIAN MODEL FITTING**
                           #' 
                           #' The grid data must contain columns `t*`, giving the case
                           #' count in each time period (see `points_to_grid`), or a column `y` in single time period cases, as well as any covariates to include in the model
                           #' (see `add_covariates`). If the population density is not provided it is set to one. If the data are regional data, then the outcome
                           #' counts must be in self$region_data
                           #'
                           #' The statistical model is a Log Gaussian cox process,
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
                           #' 
                           #' *Priors*
                           #' 
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
                           #' probabilistic programming. Stat Comput. 2023;33:17.
                           #' doi:10.1007/s11222-022-10167-2.
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
                           #' @param iter_warmup integer. Number of warmup iterations
                           #' @param iter_sampling integer. Number of sampling iterations
                           #' @param chains integer. Number of chains
                           #' @param parallel_chains integer. Number of parallel chains
                           #' @param priors list. See Details
                           #' @param verbose logical. Provide feedback on progress
                           #' @param vb Logical indicating whether to use variational Bayes (TRUE) or full MCMC sampling (FALSE)
                           #' @param return_stan_fit logical. The results of the model fit are stored internally as an `rstFit` object and 
                           #' returned in that format. If this argument is set to TRUE, then the fitted stan object will instead be returned, 
                           #' but the `rtsFit` object will still be saved. 
                           #' @param ... additional options to pass to `$sample()``.
                           #' @return A \link[rstan]{stanfit} or a `CmdStanMCMC` object
                           #' @seealso points_to_grid, add_covariates
                           #' @examples
                           #' # the data are just random simulated points 
                           #' b1 <- sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
                           #' g1 <- grid$new(b1,0.5)
                           #' dp <- data.frame(y=runif(10,0,3),x=runif(10,0,3),date=paste0("2021-01-",11:20))
                           #' dp <- create_points(dp,pos_vars = c('y','x'),t_var='date')
                           #' cov1 <- grid$new(b1,0.8)
                           #' cov1$grid_data$cov <- runif(nrow(cov1$grid_data))
                           #' g1$add_covariates(cov1$grid_data,
                           #'                   zcols="cov")
                           #' g1$points_to_grid(dp, laglength=5)
                           #' g1$priors <- list(
                           #'   prior_lscale=c(0,0.5),
                           #'   prior_var=c(0,0.5),
                           #'   prior_linpred_mean=c(0),
                           #'   prior_linpred_sd=c(5)
                           #'   )
                           #' \donttest{
                           #' g1$lgcp_bayes(popdens="cov", approx = "hsgp", parallel_chains = 0)
                           #' g1$model_fit()
                           #' # we can extract predictions
                           #' g1$extract_preds("rr")
                           #' g1$plot("rr")
                           #' g1$hotspots(threshold = 2, stat = "rr")
                           #' 
                           #'  # this example uses real aggregated data but will take a relatively long time to run
                           #'  data("birmingham_crime")
                           #'  example_data <- birmingham_crime[,c(1:8,21)]
                           #'  example_data$y <- birmingham_crime$t12
                           #'  g2 <- grid$new(example_data,cellsize=1000)
                           #'  g2$priors <- list(
                           #'   prior_lscale=c(0,0.5),
                           #'   prior_var=c(0,0.5),
                           #'   prior_linpred_mean=c(-3),
                           #'   prior_linpred_sd=c(5)
                           #' )
                           #' g2$lgcp_bayes(popdens="pop", approx = "hsgp", parallel_chains = 0)
                           #' g2$model_fit()
                           #' g2$extract_preds("rr")
                           #' g2$plot("rr")
                           #' g2$hotspots(threshold = 2, stat = "rr")
                           #' }
                           lgcp_bayes = function(popdens=NULL,
                                               covs=NULL,
                                               covs_grid = NULL,
                                               approx = "nngp",
                                               m=10,
                                               L=1.5,
                                               model = "exp",
                                               known_theta = NULL,
                                               iter_warmup=500,
                                               iter_sampling=500,
                                               chains=3,
                                               parallel_chains=3,
                                               verbose=TRUE,
                                               vb = FALSE,
                                               return_stan_fit = FALSE,
                                               ...){
                             if(!approx%in%c('nngp','hsgp','none'))stop("approx must be one of nngp, hsgp, or none")
                             if(m<=1 & approx %in% c('nngp','hsgp'))stop("m must be greater than one")
                             if(!is.null(self$region_data)&verbose)message("Using regional data model.")
                             if(is.null(popdens)){
                               self$grid_data$intercept <- 1
                               popdens <- "intercept"
                             }
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
                             if(!is.null(self$region_data)){
                               ## now do regional model variant
                               W1 <- self$get_region_data()
                               W <- Matrix::sparseMatrix(j = W1$cell_id, p = W1$n_cell-1, x = as.vector(W1$q_weights))
                               add_offset <- FALSE
                               if(!is.null(popdens)){
                                 if(popdens %in% colnames(self$grid_data)){
                                   W <- W %*% Matrix::Diagonal(x = as.data.frame(self$grid_data[,popdens])[,1])
                                   datlist$popdens <- rep(1, datlist$nT*nrow(self$grid_data))
                                   datlist$q_weights <- W@x
                                 } else if(popdens %in% colnames(self$region_data)){
                                   W <- Matrix::Diagonal(x = as.data.frame(self$region_data[,popdens])[,1]) %*% W
                                   datlist$popdens <- rep(1, datlist$nT*nrow(self$grid_data))
                                   datlist$q_weights <- W@x
                                 } 
                               }
                             }
                             if(!verbose){
                               if(!vb){
                                 capture.output(suppressWarnings( res <- rstan::sampling(stanmodels$rtsbayes,
                                                                                         data=datlist,
                                                                                         chains=chains,
                                                                                         iter = iter_warmup+iter_sampling,
                                                                                         warmup = iter_warmup,
                                                                                         cores = parallel_chains,
                                                                                         refresh = 0)), file=tempfile())
                               } else {
                                 capture.output(suppressWarnings( res <- rstan::vb(stanmodels$rtsbayes,
                                                                                   data=datlist)), file=tempfile())
                               }
                               
                             } else {
                               if(!vb){
                                 res <- rstan::sampling(stanmodels$rtsbayes,
                                                        data=datlist,
                                                        chains=chains,
                                                        iter = iter_warmup+iter_sampling,
                                                        warmup = iter_warmup,
                                                        cores = parallel_chains)
                               } else {
                                 res <- rstan::vb(stanmodels$rtsbayes,
                                                  data=datlist)
                               }
                               if(!is.null(self$region_data)){
                                 ypred <- rstan::extract(res,"region_predict")$region_predict
                               } else {
                                 ypred <- rstan::extract(res,"y_grid_predict")$y_grid_predict
                               }
                               
                               f <- rstan::extract(res,"f")$f
                               gamma <- rstan::extract(res,"gamma")$gamma
                               if(is.null(known_theta)){
                                 phi <- rstan::extract(res,"phi_param")$phi_param
                                 sigma <- rstan::extract(res,"sigma_param")$sigma_param
                               } 
                               #if(!is.null(self$region_data)) gamma_g <- rstan::extract(res,"gamma_g")$gamma_g
                               if(datlist$nT > 1) ar <- rstan::extract(res,"ar")$ar
                             }
                             
                             pars <- c("(Intercept)",covs)
                             if(!is.null(self$region_data) & length(covs_grid)>0) pars <- c(pars, covs_grid)
                             pars <- c(pars,"sigma","phi")
                             if(datlist$nT > 1)pars <- c(pars, "rho")
                             ests <- colMeans(gamma)
                             #if(!is.null(self$region_data)& length(covs_grid)>0) ests <- c(ests, colMeans(gamma_g))
                             ests <- c(ests, colMeans(sigma),colMeans(phi))
                             if(datlist$nT > 1) ests <- c(ests, colMeans(ar))
                             sds <- apply(gamma,2,sd)
                             #if(!is.null(self$region_data)& length(covs_grid)>0) sds <- apply(gamma_g,2,sd)
                             sds <- c(sds,apply(sigma,2,sd), apply(phi,2,sd))
                             if(datlist$nT > 1) sds <- c(sds, apply(ar,2,sd))
                             lower <- apply(gamma,2,function(x)quantile(x,0.025))
                             #if(!is.null(self$region_data)& length(covs_grid)>0) lower <- c(lower, apply(gamma_g,2,function(x)quantile(x,0.025)))
                             lower <- c(lower,apply(sigma,2,function(x)quantile(x,0.025)), apply(phi,2,function(x)quantile(x,0.025)))
                             if(datlist$nT > 1) lower <- c(lower, apply(ar,2,function(x)quantile(x,0.025)))
                             upper <- apply(gamma,2,function(x)quantile(x,0.975))
                             #if(!is.null(self$region_data)& length(covs_grid)>0) upper <- c(upper, apply(gamma_g,2,function(x)quantile(x,0.975)))
                             upper <- c(upper,apply(sigma,2,function(x)quantile(x,0.975)), apply(phi,2,function(x)quantile(x,0.975)))
                             if(datlist$nT > 1) upper <- c(upper, apply(ar,2,function(x)quantile(x,0.975)))
                               
                             results <- data.frame(par = pars,
                                               est = ests,
                                               SE= sds,
                                               t = NA,
                                               p = NA,
                                               lower = lower,
                                               upper = upper)
                             
                             out <- list(coefficients = results,
                                         converged = NA,
                                         approx = approx,
                                         method = ifelse(vb,"vb","mcmc"),
                                         m = iter_sampling*chains,
                                         tol = NA,
                                         aic = NA,
                                         se=NA,
                                         Rsq = NA,
                                         logl = NA,
                                         re.samps = t(f),
                                         iter = iter_sampling,
                                         time = NA,
                                         P = 1 + length(covs) + length(covs_grid),
                                         region = !is.null(self$region_data),
                                         covs = covs,
                                         y=datlist$y,
                                         y_predicted  = t(ypred),
                                         se_pred = NULL,
                                         nT = datlist$nT,
                                         conv_criterion = NA, 
                                         weights = rep(1/nrow(f),nrow(f)))
                             class(out) <- "rtsFit"
                             private$last_model_fit <- out
                             if(return_stan_fit){
                               return(res)
                             } else {
                               return(invisible(out))
                             }
                           },
                           #' @description
                           #' Fit an (approximate) log-Gaussian Cox Process model using Maximum Likelihood
                           #'
                           #' @details
                           #' **MAXIMUM LIKELIHOOD MODEL FITTING**
                           #' 
                           #' The grid or region data must contain columns `t*`, giving the case
                           #' count in each time period (see `points_to_grid`), or `y` if single time period, as well as any covariates to include in the model
                           #' (see `add_covariates`). If a population density variable is not provided it is set to one. If the data are regional data then the outcome
                           #' counts must be in self$region_data. See `lgcp_bayes()` for Bayesian approaches to model fitting and more details on the model.
                           #' 
                           #' Model fitting uses a fast stochastic maximum likelihood algorithms, which have three steps:
                           #'  1. Sample random effects using MCMC. The argument  
                           #'     `mcmc_sampling` specifies the iterations for this step. 
                           #'  2. Fit fixed effect parameters using a Newton-Raphson step.
                           #'  3. Fit covariance parameters using a Newton-Raphson step. 
                           #'     
                           #'  The algorithm uses Bayes Factor termination criterion to measure the evidence of convergence. The algorithm terminates when
                           #'  the Bayes Factor is greater than `tol` (default is 10). The prior is based on the expected number of iterations until convergence, 
                           #'  given by `k0` (defaul is 10).
                           #' 
                           #' @param popdens character string. Name of the population density column in grid data or region data
                           #' @param covs vector of strings. Base names of the covariates to
                           #' include. For temporally-varying covariates only the stem is required and not
                           #' the individual column names for each time period (e.g. `dayMon` and not `dayMon1`,
                           #' `dayMon2`, etc.) For regional models, covariates should be mapped to the grid currently (see add_covariates)
                           #' @param model Either "fexp", "sqexp", "matern1", or "matern2" for exponential, squared exponential, and matern with shape of 1 or 2, respectively. Add `hsgp_` for hsgp approximation. Other functions 
                           #' from glmmrBase will work, but may be less relevant to spatial models.
                           #' @param iter_sampling integer. Number of random effects samples to draw on each iteration.
                           #' @param max_iter Integer. Maximum number of iterations.
                           #' @param tol Scalar. The tolerance for the Bayes Factor convergence criterion.
                           #' @param hist Integer. The length of the running mean for the convergence criterion for non-Gaussian models.
                           #' @param k0 Integer. The expected nunb2mber of iterations until convergence.
                           #' @param trace Integer. Level of detail of information printed to the console. 0 = none, 1 = some (default), 2 = most.
                           #' @param start_theta Optional. Starting values for the covariance parameters (log(tau sq), log(lambda), rho), with rho only 
                           #' required if more than one time period.
                           #' @param m Number of basis functions per dimension for HSGP approximation. If spatio-temporal model then three values are required.
                           #' @param L Multiplicative boundary extension for HSGP approximation.
                           #' @param se Either "average" for a Monte Carlo averaged standard error, or "point" for standard errors evaluated at the posterior mode of the random effects.
                           #' @param mesh_args If using spde, then a list with any of `max_edge`, `cutoff`, and `offset` specified, see help for `build_spde_data` function.
                           #' @return Optionally, an `rtsFit` model fit object. This fit is stored internally and can be retrieved with `model_fit()`
                           #' @seealso points_to_grid, add_covariates
                           #' @examples
                           #' # a simple example with completely random points
                           #' b1 <- sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
                           #' g1 <- grid$new(b1,0.5)
                           #' dp <- data.frame(y=runif(10,0,3),x=runif(10,0,3))
                           #' dp <- create_points(dp,pos_vars = c('y','x'))
                           #' cov1 <- grid$new(b1,0.8)
                           #' cov1$grid_data$cov <- runif(nrow(cov1$grid_data))
                           #' g1$add_covariates(cov1$grid_data,
                           #'                   zcols="cov")
                           #' g1$points_to_grid(dp)
                           #' 
                           #' # an example using real data
                           #' 
                           #' \donttest{
                           #' g1$lgcp_ml(popdens="cov",iter_sampling = 50)
                           #' g1$model_fit()
                           #' g1$extract_preds("rr")
                           #' g1$plot("rr")
                           #' g1$hotspots(threshold = 2, stat = "rr")
                           #' 
                           #' # this example uses real aggregated spatial data
                           #' # note that the full dataset has 12 time periods
                           #' # and can be used as a spatio-temporal example by removing
                           #' # the lines marked # spatial 
                           #'  data("birmingham_crime")
                           #'  example_data <- birmingham_crime[,c(1:8,21)] # spatial
                           #'  example_data$y <- birmingham_crime$t12 # spatial
                           #'  example_data <- sf::st_transform(example_data, crs = 4326)
                           #'  g2 <- grid$new(example_data,cellsize=0.006)
                           #'  g2$lgcp_ml(popdens = "pop", model = "hsgp_matern1", L = 1.5, m = c(15,15))
                           #'  g2$extract_preds("rr")
                           #'  g2$plot("rr")
                           #' }
                           #' 
                           lgcp_ml = function(popdens=NULL,
                                              covs=NULL,
                                              model = "fexp",
                                              max_iter = 30,
                                              iter_sampling=200,
                                              tol = ifelse(iter_sampling == 1,1,10), 
                                              hist = 5, 
                                              k0 = 10,
                                              start_theta = NULL,
                                              m = c(10,10),
                                              L = 1.2,
                                              se = "average",
                                              mesh_args = NULL,
                                              trace = 1){
                             
                             # some checks at the beginning
                             if(!is.null(self$region_data)& trace >= 1)message("Using regional data model.")
                             if(is.null(popdens)){
                               self$grid_data$intercept <- 1
                               popdens <- "intercept"
                             }
                             
                             if(!is.null(mesh_args)){
                               if(!is("mesh_args",list))stop("mesh_args should be a list if using SPDE")
                               if(! "max_edge" %in% names(mesh_args)) mesh_args <- append(mesh_args, list(max_edge = NULL))
                               if(! "cutoff" %in% names(mesh_args)) mesh_args <- append(mesh_args, list(cutoff = NULL))
                               if(! "offset" %in% names(mesh_args)) mesh_args <- append(mesh_args, list(offset = NULL))
                             } else {
                               mesh_args <- list(max_edge = NULL, cutoff = NULL, offset = NULL)
                             }
                             
                             bb <- sf::st_bbox(self$boundary)
                             span <- min(bb[3] - bb[1], bb[4] - bb[2])
                             if(is.null(start_theta)) start_theta <- log(c(rnorm(1,1,0.05), span/5))
                             
                             data <- private$prepare_data(1,
                                      model,
                                      "none",
                                      popdens,
                                      covs,
                                      NULL,
                                      TRUE,
                                      FALSE,
                                      1)
                             
                             form <- "~"
                             if(!is.null(covs)){
                               for(i in 1:length(covs)){
                                 form <- paste0(form,covs[i],"+")
                               }
                             }
                             is_hsgp <- grepl("hsgp",model)
                             is_spde <- grepl("spde",model)
                             if(packageVersion(pkg = "glmmrBase") >= "1.2.0"){
                               if(data$nT > 1 & !is_hsgp){
                                 if(is.null(self$region_data)){
                                   form <- paste0(form,"(1|ar","_",model,"log(x,y,t=",data$nT,"))")
                                 } else {
                                   form <- paste0(form,"(1|",model,"log(x,y))")
                                 }
                               } else {
                                 form <- paste0(form,"(1|",model,"log(x,y))")
                               }
                             } else {
                               form <- paste0(form,"(1|",model,"(x,y))")
                               print(form)
                               cat("If you're seeing this message then you need to update glmmrBase\n")
                             }
                             
                             
                             if(is.null(self$region_data)){
                               # ── Build the data frame (existing code, unchanged) ─────────────────────────
                               dat <- as.data.frame(data$X)
                               dat <- cbind(dat, data$x_grid)
                               dat <- dat[, 2:ncol(dat)]
                               rownames(dat) <- 1:nrow(dat)
                               colnames(dat) <- c(covs, "x", "y")
                               if (is_hsgp & data$nT > 1) {
                                 dat <- cbind(dat, rep(1:data$nT, each = nrow(data$x_grid)))
                                 colnames(dat) <- c(covs, "x", "y", "t")
                                 form <- gsub("\\)", ",t)", form)
                               }
                               
                               # ── Model construction (NEW SPDE branch + existing default) ─────────────────
                               spde <- NULL
                               if (is_spde) {
                                 if (data$nT > 1) stop("Spatio-temporal SPDE not yet implemented for non-aggregated LGCP")
                                 
                                 events <- sf::st_coordinates(self$point_data)[, 1:2, drop = FALSE]
                                 spde   <- do.call(self$build_spde_data, mesh_args)
                                 #spde   <- self$build_spde_data()
                                 n_v    <- spde$n_v
                                 n_e    <- nrow(events)
                                 vertex_coords <- spde$mesh$loc[, 1:2]
                                 
                                 aug <- data.frame(
                                   x     = c(events[, 1], vertex_coords[, 1]),
                                   y     = c(events[, 2], vertex_coords[, 2]),
                                   count = c(rep(1L, n_e), rep(0L, n_v)),
                                   log_w = c(rep(0,   n_e), log(spde$C))
                                 )
                                 if (!is.null(popdens)) {
                                   pop_aug <- private$evaluate_covariates_at(aug[, c("x", "y")], covs = popdens)[, popdens]
                                   aug$log_w <- aug$log_w + log(pop_aug)
                                 }
                                 cov_aug <- private$evaluate_covariates_at(aug[, c("x", "y")], covs = covs)
                                 for (cv in covs) aug[[cv]] <- cov_aug[, cv]
                                 
                                 if(!is.null(covs)){
                                   spde_form <- as.formula(paste0("count ~ ", paste(covs, collapse = " + "),
                                                                  " + (1 | spde_matern1log(x, y))"))
                                 } else {
                                   spde_form <- as.formula("count ~ (1 | spde_matern1log(x, y))")
                                 }
                                
                                 start_cov <- start_theta
                                 
                                 mod <- glmmrBase::Model$new(
                                   formula    = spde_form,
                                   data       = aug,
                                   family     = poisson(),
                                   mesh       = list(A_loc = spde$A_loc, C = spde$C, G = spde$G),
                                   covariance = start_cov,
                                   mean       = c(log(n_e / sum(spde$C)), rep(0,length(covs))),
                                   offset     = aug$log_w
                                 )$set_trace(1)
                               } else {
                                 mod <- glmmrBase::Model$new(
                                   formula = as.formula(form),
                                   data    = dat,
                                   family  = poisson(),
                                   offset  = log(data$popdens)
                                 )$set_trace(1)
                                 
                                 if (is.null(start_theta)) {
                                   start_cov <- if (is_hsgp & data$nT > 1) log(runif(4)) else log(runif(2 + I(data$nT > 1) * 1))
                                 } else {
                                   start_cov <- start_theta
                                 }
                                 
                                 if(length(start_theta) == 2 & data$nT > 1){
                                   if(is_hsgp){
                                     start_theta <- c(start_theta, start_theta[2], log(diff(range(data$nT))/4), log(runif(1)))
                                   } else {
                                     start_theta <- c(start_theta, log(runif(1)))
                                   }
                                 }
                               }
                               
                               # ── Fit (existing logic, with type-aware y) ─────────────────────────────────
                               mod$update_y(if (is_spde) aug$count else data$y)
                               
                               mod$update_parameters(cov.pars = start_cov)
                               fit <- mod$fit(niter = iter_sampling, max_iter = max_iter, tol = tol,
                                              hist = hist, k0 = k0, se = se)
                               ll <- mod$log_likelihood()
                               X  <- mod$mean$X
                               n_cov_pars <- 2
                               if (data$nT > 1 && !is_hsgp && !is_spde) n_cov_pars <- 3
                               if (data$nT > 1 &&  is_hsgp)             n_cov_pars <- 4
                               cov_par_names <- c("tau_sq (log)", "lambda (log)")
                               if (data$nT > 1 && !is_hsgp && !is_spde) cov_par_names <- c(cov_par_names, "rho")
                               if (data$nT > 1 &&  is_hsgp)             cov_par_names <- c(cov_par_names, "lambda_2 (log)", "lambda_t (log)")
                               fit$coefficients$par[(ncol(X) + 1):(ncol(X) + n_cov_pars)] <- cov_par_names
                               
                               # ── Common extractions ──────────────────────────────────────────────────────
                               ptr   <- mod$.__enclos_env__$private$ptr
                               mtype <- mod$.__enclos_env__$private$model_type()
                               
                               popd  <- private$stack_variable(popdens)
                               beta  <- fit$coefficients$est[1:ncol(X)]
                               u_obs <- mod$u(scaled = TRUE)  
                               if (packageVersion(pkg = "glmmrBase") >= "1.2.0") {
                                 w <- glmmrBase:::Model__get_importance_weights(ptr, mtype)
                               } else {
                                 w <- rep(1 / ncol(u_obs), ncol(u_obs))
                               }
                               sum_w <- sum(w); K <- ncol(u_obs)
                               M <- solve(mod$information_matrix())                # V_β (marginal Schur)
                               if (packageVersion(pkg = "glmmrBase") >= "1.2.0") {
                                 M_theta <- solve(mod$information_matrix(theta = TRUE))
                               } else {
                                 M_theta <- diag(2)
                               }
                               # ── Project u and compute var_rr at the GRID ────────────────────────────────
                               if(is_spde){
                                 X_grid <- data$X
                               } else {
                                 X_grid <- mod$mean$X                                     
                               }
                               n_grid <- nrow(X_grid)
                               XVX <- rowSums((X_grid %*% M) * X_grid)
                               
                               if (is_spde) {
                                 
                                 u_obs <- mod$covariance$Z %*% u_obs
                                 u_v     <- u_obs[(n_e + 1):(n_e + n_v), , drop = FALSE]   # subset to vertices
                                 u_grid  <- as.matrix(spde$A_pred %*% u_v)     
                                 var_rr  <- tryCatch(glmmrBase:::Model_spde__re_var_at(ptr, spde$A_pred),
                                  error = function(e){
                                    print(e)
                                    return(rep(1,nrow(spde$A_pred)))
                                  })
                                 u <-  u_grid  #<- as.matrix(spde$A_pred %*% mod$u(FALSE))
                               } else if (is_hsgp) {
                                 u_grid  <- u_obs                                           # obs = grid
                                 var_rr  <- glmmrBase:::Model_hsgp__re_var(ptr)             # or dense equivalent
                                 u       <- u_grid
                               } else {
                                 u_grid  <- u_obs     
                                 Mv <- mod$information_matrix(include.re = TRUE)
                                 var_rr <- diag(solve(Mv)[(ncol(X)+1):nrow(Mv),(ncol(X)+1):nrow(Mv)])
                                 u      <- u_grid
                               }
                               
                               # ── Grid-level intensity samples and posterior summaries ────────────────────
                               xb  <- drop(X_grid %*% beta)                         # X β only (no offset)
                               
                               mu_pp_samples  <- matrix(0, n_grid, K)
                               mu_tot_samples <- matrix(0, n_grid, K)
                               for (k in 1:K) {
                                 mu_pp_samples[, k]  <- exp(xb + u_grid[, k])
                                 mu_tot_samples[, k] <- exp(xb + u_grid[, k] + log(data$popdens))
                               }
                               
                               if (K == 1) {
                                 mu_pp_mean  <- mu_pp_samples[, 1]
                                 mu_tot_mean <- mu_tot_samples[, 1]
                                 rr_mean     <- u_grid[, 1]
                               } else {
                                 mu_pp_mean  <- rowSums(t(t(mu_pp_samples)  * w)) / sum_w
                                 mu_tot_mean <- rowSums(t(t(mu_tot_samples) * w)) / sum_w
                                 rr_centered <- sweep(u_grid, 2, colMeans(u_grid))
                                 rr_mean     <- drop(rr_centered %*% w) / sum_w
                               }
                               
                               # Delta-method SEs using analytical var_rr (consistent across K)
                               SE_rr  <- sqrt(var_rr)
                               SEeta <- sqrt(var_rr + XVX)         # η-uncertainty — conservative delta
                               SEpp  <- mu_pp_mean  * SEeta
                               SEtot <- mu_tot_mean * SEeta
                               
                               ypred  <- mu_tot_mean
                               mupred <- mu_pp_mean
                             } else {
                               ## now do regional model variant
                               W <- self$get_region_data()$W
                               
                               # ── SPDE-specific setup: replace grid-based W with mesh-based W ──
                               spde_data_obj <- NULL
                               if (is_spde) {
                                 #spde_data_obj <- self$build_spde_data()
                                 spde_data_obj <- do.call(self$build_spde_data, mesh_args)
                                 W <- spde_data_obj$W_mesh   # n_regions × n_v, replaces grid W on the model side
                               }
                               
                               add_offset <- FALSE
                               pop_in_W   <- FALSE
                               
                               if (!is.null(popdens)) {
                                 if (popdens %in% colnames(self$grid_data)) {
                                   if (is_spde) {
                                     pop_at_mesh <- private$evaluate_covariates_at(
                                       spde_data_obj$mesh$loc[, 1:2], popdens
                                     )[, 1]
                                     W <- W %*% Matrix::Diagonal(x = pop_at_mesh)
                                   } else {
                                     W <- W %*% Matrix::Diagonal(x = as.data.frame(self$grid_data[, popdens])[,1])
                                   }
                                   pop_in_W <- TRUE
                                 } else if (popdens %in% colnames(self$region_data)) {
                                   W <- Matrix::Diagonal(x = as.data.frame(self$region_data[, popdens])[,1]) %*% W
                                   pop_in_W <- TRUE
                                 } else if (length(data$popd) == (ncol(W) * data$nT)) {
                                   add_offset <- TRUE
                                 } else {
                                   stop("popdens not in grid/region, or is time varying with SPDE.")
                                 }
                               }
                               
                               if (is_spde) {
                                 X_v <- matrix(1, nrow = nrow(spde_data_obj$mesh$loc), ncol = 1)
                                 if (!is.null(covs)) {
                                   aug_cov <- private$evaluate_covariates_at(spde_data_obj$mesh$loc[, 1:2], covs)
                                   for (i in 1:ncol(aug_cov)) X_v <- cbind(X_v, aug_cov[, i])
                                 }
                                 offset <- rep(0, nrow(spde_data_obj$mesh$loc))   # always zero now
                                 if (add_offset) stop("Time-varying popdens not yet supported with SPDE.")
                               } else {
                                 X_v <- data$X
                                 offset <- if (add_offset) log(data$popd) else rep(0, nrow(X_v))
                               }
                               
                               # parse model class
                               if (is_hsgp) form <- gsub("hsgp_","", form)
                               if (is_spde) form <- gsub("spde_","", form)
                               if (is_hsgp && is_spde) stop("Model cannot be both HSGP and SPDE")
                               form <- gsub("~","", form)
                               
                               if (is_hsgp) {
                                 type <- 2L
                               } else if (is_spde) {
                                 type <- 3L
                               } else {
                                 type <- I(data$nT > 1) * 1
                               }
                               
                               # Covariates evaluated at mesh vertices for SPDE; X stays as the design matrix
                               # the C++ side multiplies against the fitted random effects (n_v × niter). The
                               # user-side responsibility is to supply X_v of dim (n_v × p). For now we evaluate
                               # covariates at vertices using the same column logic as data$X on the grid.
                               if (is_hsgp) {
                                 if (data$nT == 1) {
                                   if(length(m) != 2)stop("m must be length 2")
                                   ptr <- regionModel_hsgp__new(
                                     form, as.matrix(data$x_grid), c('x','y'),
                                     data$X, data$y, iter_sampling, m, L
                                   )
                                   theta_start <- start_theta
                                 } else {
                                   if(length(m) != 3)stop("m must be length 3")
                                   form  <- gsub("\\)", ",t)", form)
                                   gcsize <- nrow(data$x_grid)
                                   gcmat  <- data$x_grid
                                   for (t in 1:(data$nT-1)) gcmat <- rbind(gcmat, data$x_grid)
                                   gcmat  <- cbind(gcmat, data.frame(t = rep(1:data$nT, each = gcsize)))
                                   ptr <- regionModel_hsgp__new(
                                     form, as.matrix(gcmat), c('x','y','t'),
                                     data$X, data$y, iter_sampling, m, L
                                   )
                                   if(length(start_theta) == 2 & data$nT > 1){
                                     theta_start <- c(start_theta, start_theta[2], log(diff(range(data$nT))/4), log(runif(1)))
                                   }
                                 }
                               } else if (is_spde) {
                                 if (data$nT > 1) stop("Spatio-temporal SPDE not yet implemented")
                                 ptr <- regionModel_spde__new(
                                   form,
                                   as.matrix(spde_data_obj$mesh$loc[, 1:2]),
                                   c('x','y'),
                                   X_v,                  # covariates at mesh vertices
                                   data$y,               # region-level counts (length n_regions)
                                   iter_sampling,
                                   spde_data_obj$A_loc,  # n_v × n_v identity (model-side projector)
                                   spde_data_obj$C,
                                   spde_data_obj$G,
                                   2L                    # alpha = 2 → ν = 1 in 2D
                                 )
                                 theta_start <- start_theta
                               } else {
                                 if (data$nT == 1) {
                                   ptr <- regionModel__new(
                                     form, as.matrix(data$x_grid), c('x','y'),
                                     data$X, data$y, iter_sampling, 0
                                   )
                                   theta_start <- start_theta
                                 } else {
                                   ptr <- regionModel_ar__new(
                                     form, as.matrix(data$x_grid), c('x','y'),
                                     data$X, data$y, iter_sampling, data$nT
                                   )
                                   if(length(start_theta) == 2) start_theta <- c(start_theta, log(runif(1)))
                                 }
                               }
                               
                               regionModel__set_theta(ptr, theta_start, type)
                               regionModel__set_weights(ptr, W@i, W@p, W@x, nrow(W), ncol(W), type)
                               regionModel__set_offset(ptr, offset, type)
                               regionModel__fit(ptr, tol, max_iter, hist, k0, type)
                               I_beta <- regionModel__information_matrix(ptr, se == "average", type)
                               V_beta <- solve(I_beta)
                               M      <- I_beta
                               M_theta <- regionModel__information_matrix_theta(ptr, type)
                               V_theta <- solve(M_theta)
                               u_model <- regionModel__u(ptr, TRUE, type)             # n_eval × K  (n_eval = n_cells or n_v)
                               beta    <- regionModel__get_beta(ptr, type)
                               theta   <- regionModel__get_theta(ptr, type)
                               # ── For SPDE, lift to grid for prediction; for grid models, u_model is already at grid ──
                               if (is_spde) {
                                 # Grid centroids — used for both covariate evaluation and offset
                                 grid_centroids <- as.data.frame(
                                   sf::st_coordinates(sf::st_centroid(self$grid_data))
                                 )
                                 
                                 # Project u from mesh to grid via A_pred (n_grid × n_v)
                                 u_s <- regionModel__u(ptr, FALSE, type)              # n_v × K
                                 u   <- as.matrix(spde_data_obj$A_pred %*% u_s)        # n_grid × K
                                 
                                 # Build X_pr at grid cells
                                 X_pr <- matrix(1, nrow = nrow(spde_data_obj$A_pred), ncol = 1)
                                 if (!is.null(covs)) {
                                   aug_cov <- private$evaluate_covariates_at(sf::st_coordinates(sf::st_centroid(self$grid_data)), covs)
                                   X_pr    <- cbind(X_pr, aug_cov)
                                 }
                                 
                                 # Grid-level offset.
                                 # If popdens was absorbed into W_mesh during fitting, we already account
                                 # for it via W_out scaling below; offset stays zero at the grid level.
                                 # If popdens is time-varying, that case needs separate handling and isn't
                                 # supported with SPDE yet.
                                 offset <- rep(0, nrow(X_pr))
                                 
                               } else {
                                 u      <- u_model
                                 X_pr   <- data$X
                                 # offset comes through from upstream as already grid-sized
                               }
                               
                               # From here, the code is SPDE-agnostic ------------------------------------
                               n <- nrow(X_pr); p <- ncol(X_pr)
                               K <- ncol(u)
                               m <- nrow(self$get_region_data()$W)   # or wherever W_out comes from
                               
                               stopifnot(nrow(u) == n)               # cheap guard against the bug class
                               
                               cov_par_names <- c("tau_sq (log)","lambda (log)")
                               if (data$nT > 1 && !is_hsgp && !is_spde) cov_par_names <- c(cov_par_names, "rho")
                               if (data$nT > 1 &&  is_hsgp)             cov_par_names <- c(cov_par_names, "lambda_2 (log)", "lambda_t (log)")
                               
                               res <- data.frame(
                                 par   = c(paste0("beta",1:length(beta)), cov_par_names, paste0("d",1:nrow(u))),
                                 est   = c(beta, theta, rowMeans(u)),
                                 SE    = c(sqrt(diag(V_beta)), sqrt(diag(V_theta)), rep(NA, nrow(u))),
                                 t = NA, p = NA, lower = NA, upper = NA
                               )
                               res$t <- res$est / res$SE
                               res$p <- 2 * (1 - stats::pnorm(abs(res$t)))
                               res$lower <- res$est - qnorm(1 - 0.05/2) * res$SE
                               res$upper <- res$est + qnorm(1 - 0.05/2) * res$SE
                               fit <- list(coefficients = res, converged = TRUE, aic = NA, Rsq = list(NA, NA), iter = NA)
                               # ── Output slice operates on grid-level (X_pr, u) for all model types ──
                               n <- nrow(X_pr); p <- ncol(X_pr); m <- nrow(W); K <- ncol(u)
                               
                               # W for output aggregation (region-from-grid). For SPDE this is the *original*
                               # grid-based W from get_region_data(), not the mesh-based W used in fitting.
                               W_out <- self$get_region_data()$W
                               if (!is.null(popdens)) {
                                 if (popdens %in% colnames(self$grid_data)) {
                                   W_out <- W_out %*% Matrix::Diagonal(x = as.data.frame(self$grid_data[, popdens])[,1])
                                 } else if (popdens %in% colnames(self$region_data)) {
                                   W_out <- Matrix::Diagonal(x = as.data.frame(self$region_data[, popdens])[,1]) %*% W_out
                                 }
                               }
                               w <- regionModel__get_weights(ptr, type)
                               sum_w <- sum(w)
                               xb    <- drop(X_pr %*% beta + offset)
                               XVX   <- rowSums((X_pr %*% V_beta) * X_pr)
                               mu_samples <- matrix(0, n, K)
                               for (k in 1:K) mu_samples[, k] <- exp(xb + u[, k])
                               
                               # ── Output slice (unified MC path) ─────────────────────────────────────────────
                               n <- nrow(X_pr); p <- ncol(X_pr); m <- nrow(W); K <- ncol(u)
                               
                               # W for region-from-grid output (SPDE: original grid-based W; grid models: same as fitting W)
                               W_out <- self$get_region_data()$W
                               # if (is_spde && !is.null(popdens)) {
                               #   if (popdens %in% colnames(self$grid_data)) {
                               #     W_out <- W_out %*% Matrix::Diagonal(x = as.data.frame(self$grid_data[, popdens])[,1])
                               #   } else if (popdens %in% colnames(self$region_data)) {
                               #     W_out <- Matrix::Diagonal(x = as.data.frame(self$region_data[, popdens])[,1]) %*% W_out
                               #   }
                               # }
                               
                               w     <- regionModel__get_weights(ptr, type)
                               sum_w <- sum(w)
                               
                               # ── Grid-level Var(η) — single source of truth for β + u uncertainty ─────────
                               zu_var_grid <- if (is_spde) {
                                 regionModel_spde__zu_variance_full(ptr, spde_data_obj$A_pred, X_pr)
                               } else {
                                 regionModel__zu_variance_full(ptr, type)   # grid/HSGP: ZA is internal projector
                               }
                               # ── MC samples at grid level (Jensen mean-centered per sample) ───────────────
                               # NB: offset must be sized to nrow(X_pr). For SPDE the current code sets
                               # `offset <- rep(0, nrow(X_v))` (mesh-sized) — that's a pre-existing size bug
                               # in the no-offset path; here we treat `offset` as already grid-sized.
                               xb         <- drop(X_pr %*% beta + offset)
                               xb_no_off  <- drop(X_pr %*% beta)
                               mu_samples       <- matrix(0, n, K)
                               mu_samples_noff  <- matrix(0, n, K)
                               for (k in 1:K) {
                                 uk_c <- u[, k] - mean(u[, k])
                                 mu_samples[, k]      <- exp(xb        + uk_c)
                                 mu_samples_noff[, k] <- exp(xb_no_off + uk_c)
                               }
                               
                               # ── Grid-level posterior summaries ───────────────────────────────────────────
                               mu_mean <- rowSums(t(t(mu_samples)      * w)) / sum_w
                               mupred  <- rowSums(t(t(mu_samples_noff) * w)) / sum_w
                               SEpp    <- mu_mean * sqrt(zu_var_grid)              # delta on Var(η)
                               
                               # Log-relative-risk (= η minus its spatial mean across samples)
                               u_centered <- u - matrix(colMeans(u), n, K, byrow = TRUE)
                               rr_mean    <- rowSums(t(t(u_centered) * w)) / sum_w
                               SE_rr      <- sqrt(zu_var_grid)
                               
                               # ── Region-level aggregation (kept as MC; β contribution still needs V_β) ────
                               if (data$nT > 1) {
                                 n_A <- ncol(W_out); n_R <- nrow(W_out)
                                 T_periods <- nrow(mu_samples) / n_A
                                 mu_r_samples <- matrix(0, n_R * T_periods, K)
                                 for (t in 1:T_periods) {
                                   idx_in  <- ((t - 1) * n_A + 1):(t * n_A)
                                   idx_out <- ((t - 1) * n_R + 1):(t * n_R)
                                   mu_r_samples[idx_out, ] <- as.matrix(W_out %*% mu_samples[idx_in, , drop = FALSE])
                                 }
                               } else {
                                 mu_r_samples <- as.matrix(W_out %*% mu_samples)
                               }
                               mu_r_mean    <- rowSums(t(t(mu_r_samples) * w)) / sum_w
                               var_from_u_r <- rowSums(t(t((mu_r_samples - mu_r_mean)^2) * w)) / sum_w
                               
                               var_from_beta_r <- numeric(m)
                               for (k in 1:K) {
                                 for (r in 1:m) {
                                   g_r <- drop(t(X_pr) %*% (W_out[r, ] * mu_samples[, k]))
                                   var_from_beta_r[r] <- var_from_beta_r[r] + w[k] * drop(t(g_r) %*% V_beta %*% g_r)
                                 }
                               }
                               var_from_beta_r <- var_from_beta_r / sum_w
                               SEtot <- sqrt(var_from_beta_r + var_from_u_r)
                               ypred <- mu_r_mean
                               ll <- regionModel__log_likelihood(ptr, type)
                             }

                             out <- list(coefficients = fit$coefficients,
                                         converged = fit$converged,
                                         approx = "none",
                                         method = "mcml",
                                         m = dim(u)[2],
                                         tol = tol,
                                         aic = fit$aic,
                                         se="gls",
                                         Rsq = fit$Rsq,
                                         logl = ll,
                                         re.samps = u,
                                         iter = fit$iter,
                                         time = NA,
                                         region = !is.null(self$region_data),
                                         covs = covs,
                                         vcov = M,
                                         vcov_theta = M_theta,
                                         P = length(covs)+1,
                                         var_par_family = FALSE,
                                         y=data$y,
                                         X = data$X,
                                         y_predicted  = ypred,
                                         mu_predicted = mupred,
                                         rr = rr_mean,
                                         se_pred = list(pp = SEpp, tot = SEtot, rr = SE_rr),
                                         nT = data$nT,
                                         conv_criterion = 0,
                                         weights = w)
                             class(out) <- "rtsFit"
                             private$last_model_fit <- out
                             return(invisible(out))
                           },
                           #' @description
                           #' Extract predictions
                           #'
                           #' Extract incidence and relative risk predictions. The predictions will be extracted from the last model fit. If no previous model fit then use either `lgcp_ml()` or `lgcp_bayes()`, or see 
                           #' `model_fit()` to update the stored model fit.
                           #'
                           #' @param type Vector of character strings. Any combination of "pred", "rr", and "irr", which are,
                           #' posterior mean incidence (overall and population standardised), relative risk,
                           #' and incidence rate ratio, respectively.
                           #' @param irr.lag integer. If "irr" is requested as `type` then the number of time
                           #' periods lag previous the ratio is in comparison to
                           #' @param t.lag integer. Extract predictions for previous time periods.
                           #' @param popdens character string. Name of the column in `grid_data` with the
                           #' population density data
                           #' @return NULL
                           #' @details
                           #' **EXTRACTING PREDICTIONS**
                           #' 
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
                           #' # See examples for lgcp_bayes() and lgcp_ml()
                           #' @importFrom stats sd
                           extract_preds = function(type = c("pred","rr","irr"),
                                                    irr.lag = NULL,
                                                    t.lag = 0,
                                                    popdens = NULL,
                                                    level = 0.95){
                             if ("irr" %in% type & is.null(irr.lag)) stop("For irr set irr.lag")
                             if (is.null(self$region_data) & "pred" %in% type & is.null(popdens)) popdens <- "intercept"
                             
                             nCells  <- nrow(self$grid_data)
                             nRegion <- ifelse(is.null(self$region_data), 0, nrow(self$region_data))
                             nT      <- nrow(private$last_model_fit$re.samps) / nCells
                             alpha   <- 1 - level
                             z_crit  <- qnorm(1 - alpha / 2)
                             
                             is_sampler <- private$last_model_fit$method %in% c("vb", "mcmc")
                             
                             cell_idx <- function(t.offset = 0){
                               t_pos <- nT - t.offset
                               ((t_pos - 1) * nCells + 1):(t_pos * nCells)
                             }
                             region_idx <- function(t.offset = 0){
                               t_pos <- nT - t.offset
                               ((t_pos - 1) * nRegion + 1):(t_pos * nRegion)
                             }
                             
                             # Fast path: pull pre-computed point estimates and SEs from the fit
                             fit      <- private$last_model_fit
                             w        <- fit$weights
                             se_pred  <- fit$se_pred
                             
                             # Helpers ---------------------------------------------------------------
                             # Log-scale mean + SE -> response-scale mean and symmetric log-scale CI
                             set_log <- function(df, prefix, log_mean, se_log) {
                               df[[paste0(prefix)]]          <- exp(log_mean)
                               df[[paste0(prefix, "_se")]]   <- se_log * exp(log_mean)  # delta method
                               df[[paste0(prefix, "_lower")]] <- exp(log_mean - z_crit * se_log)
                               df[[paste0(prefix, "_upper")]] <- exp(log_mean + z_crit * se_log)
                               df
                             }
                             # Response-scale mean + SE -> symmetric Wald CI
                             set_lin <- function(df, prefix, mu, se) {
                               df[[paste0(prefix)]]          <- mu
                               df[[paste0(prefix, "_se")]]   <- se
                               df[[paste0(prefix, "_lower")]] <- mu - z_crit * se
                               df[[paste0(prefix, "_upper")]] <- mu + z_crit * se
                               df
                             }
                             
                             # For the sampler case we still need to summarise draws
                             sampler_summary <- function(mat) {
                               list(mean  = apply(mat, 1, mean),
                                    sd    = apply(mat, 1, sd),
                                    lower = apply(mat, 1, quantile, probs = alpha/2),
                                    upper = apply(mat, 1, quantile, probs = 1 - alpha/2))
                             }
                             
                             # ================================================================
                             # nT > 1
                             # ================================================================
                             if (nT > 1) {
                               ci <- cell_idx(t.lag)
                               ri <- region_idx(t.lag)
                               
                               # ---------- pred ----------
                               if ("pred" %in% type) {
                                 if (is.null(self$region_data)) {
                                   popd <- as.data.frame(self$grid_data)[, popdens]
                                   if (is_sampler) {
                                     y_samps <- fit$y_predicted[ci, , drop = FALSE]
                                     fmu     <- y_samps / popd
                                     st <- sampler_summary(y_samps); sp <- sampler_summary(fmu)
                                     self$grid_data$pred_mean_total       <- st$mean
                                     self$grid_data$pred_mean_total_sd    <- st$sd
                                     self$grid_data$pred_mean_total_lower <- st$lower
                                     self$grid_data$pred_mean_total_upper <- st$upper
                                     self$grid_data$pred_mean_pp          <- sp$mean
                                     self$grid_data$pred_mean_pp_sd       <- sp$sd
                                     self$grid_data$pred_mean_pp_lower    <- sp$lower
                                     self$grid_data$pred_mean_pp_upper    <- sp$upper
                                   } else {
                                     # mu_predicted is the per-cell rate; log it to pair with se_pred$pp (log scale)
                                     self$grid_data <- set_lin(self$grid_data, "pred_mean_pp",
                                                               fit$mu_predicted[ci], se_pred$pp[ci])
                                     self$grid_data <- set_lin(self$grid_data, "pred_mean_total",
                                                               fit$y_predicted[ci], se_pred$tot[ci])
                                   }
                                 } else {
                                   if (is_sampler) {
                                     popd    <- as.data.frame(self$region_data)[, popdens]
                                     y_samps <- fit$y_predicted[ri, , drop = FALSE]
                                     fmu     <- y_samps / popd
                                     st <- sampler_summary(y_samps); sp <- sampler_summary(fmu)
                                     self$region_data$pred_mean_total       <- st$mean
                                     self$region_data$pred_mean_total_sd    <- st$sd
                                     self$region_data$pred_mean_total_lower <- st$lower
                                     self$region_data$pred_mean_total_upper <- st$upper
                                     self$region_data$pred_mean_pp          <- sp$mean
                                     self$region_data$pred_mean_pp_sd       <- sp$sd
                                     self$region_data$pred_mean_pp_lower    <- sp$lower
                                     self$region_data$pred_mean_pp_upper    <- sp$upper
                                   } else {
                                     # Grid-level per-cell rate
                                     self$grid_data <- set_lin(self$grid_data, "pred_mean_pp",
                                                               fit$mu_predicted[ci], se_pred$pp[ci])
                                     # Region-level totals: y_predicted / se_pred$tot are assumed region-indexed
                                     self$region_data <- set_lin(self$region_data, "pred_mean_total",
                                                                 fit$y_predicted[ri], se_pred$tot[ri])
                                   }
                                 }
                               }
                               
                               # ---------- rr ----------
                               if ("rr" %in% type) {
                                 if (is_sampler) {
                                   re_samps <- fit$re.samps[ci, , drop = FALSE]
                                   self$grid_data$rr       <- exp(apply(re_samps, 1, mean))
                                   self$grid_data$rr_sd    <- exp(apply(re_samps, 1, sd))
                                   self$grid_data$rr_lower <- apply(exp(re_samps), 1, quantile, probs = alpha/2)
                                   self$grid_data$rr_upper <- apply(exp(re_samps), 1, quantile, probs = 1 - alpha/2)
                                 } else {
                                   # fit$rr is the log RR mean, se_pred$rr is its SE (also log scale)
                                   self$grid_data <- set_log(self$grid_data, "rr",
                                                             fit$rr[ci], se_pred$rr[ci])
                                 }
                               }
                               
                               # ---------- irr ----------
                               if ("irr" %in% type) {
                                 ci_lag <- cell_idx(t.lag + irr.lag)
                                 ri_lag <- region_idx(t.lag + irr.lag)
                                 
                                 if (is.null(self$region_data)) {
                                   if (is_sampler) {
                                     irr_samps <- fit$y_predicted[ci, , drop = FALSE] /
                                       fit$y_predicted[ci_lag, , drop = FALSE]
                                     s <- sampler_summary(irr_samps)
                                     self$grid_data$irr       <- s$mean
                                     self$grid_data$irr_sd    <- s$sd
                                     self$grid_data$irr_lower <- s$lower
                                     self$grid_data$irr_upper <- s$upper
                                   } else {
                                     # No stored SE for IRR: approximate log-IRR SE from the two log-rate SEs
                                     # (treats them as independent — conservative/approximate)
                                     log_irr <- log(fit$mu_predicted[ci]) - log(fit$mu_predicted[ci_lag])
                                     se_log  <- sqrt(se_pred$pp[ci]^2 + se_pred$pp[ci_lag]^2)
                                     self$grid_data <- set_log(self$grid_data, "irr", log_irr, se_log)
                                   }
                                 } else {
                                   if (is_sampler) {
                                     irr_samps <- fit$y_predicted[ri, , drop = FALSE] /
                                       fit$y_predicted[ri_lag, , drop = FALSE]
                                     s <- sampler_summary(irr_samps)
                                     self$region_data$irr       <- s$mean
                                     self$region_data$irr_sd    <- s$sd
                                     self$region_data$irr_lower <- s$lower
                                     self$region_data$irr_upper <- s$upper
                                   } else {
                                     # No stored region-level log-rate SE; fall back to ratio of totals with
                                     # a delta-method approximation using se_pred$tot
                                     mu_t   <- fit$y_predicted[ri]
                                     mu_l   <- fit$y_predicted[ri_lag]
                                     log_irr <- log(mu_t) - log(mu_l)
                                     se_log  <- sqrt((se_pred$tot[ri]/mu_t)^2 + (se_pred$tot[ri_lag]/mu_l)^2)
                                     self$region_data <- set_log(self$region_data, "irr", log_irr, se_log)
                                   }
                                 }
                               }
                               
                               # ================================================================
                               # nT == 1
                               # ================================================================
                             } else {
                               if ("irr" %in% type) stop("cannot estimate irr as only one time period")
                               
                               if ("pred" %in% type) {
                                 if (is.null(self$region_data)) {
                                   popd <- as.data.frame(self$grid_data)[, popdens]
                                   if (is_sampler) {
                                     fmu <- fit$y_predicted / popd
                                     st <- sampler_summary(fit$y_predicted); sp <- sampler_summary(fmu)
                                     self$grid_data$pred_mean_total       <- st$mean
                                     self$grid_data$pred_mean_total_sd    <- st$sd
                                     self$grid_data$pred_mean_total_lower <- st$lower
                                     self$grid_data$pred_mean_total_upper <- st$upper
                                     self$grid_data$pred_mean_pp          <- sp$mean
                                     self$grid_data$pred_mean_pp_sd       <- sp$sd
                                     self$grid_data$pred_mean_pp_lower    <- sp$lower
                                     self$grid_data$pred_mean_pp_upper    <- sp$upper
                                   } else {
                                     self$grid_data <- set_lin(self$grid_data, "pred_mean_pp",
                                                               fit$mu_predicted, se_pred$pp)
                                     self$grid_data <- set_lin(self$grid_data, "pred_mean_total",
                                                               fit$y_predicted, se_pred$tot)
                                   }
                                 } else {
                                   if (is_sampler) {
                                     popd <- as.data.frame(self$region_data)[, popdens]
                                     fmu  <- fit$y_predicted / popd
                                     st <- sampler_summary(fit$y_predicted); sp <- sampler_summary(fmu)
                                     self$region_data$pred_mean_total       <- st$mean
                                     self$region_data$pred_mean_total_sd    <- st$sd
                                     self$region_data$pred_mean_total_lower <- st$lower
                                     self$region_data$pred_mean_total_upper <- st$upper
                                     self$region_data$pred_mean_pp          <- sp$mean
                                     self$region_data$pred_mean_pp_sd       <- sp$sd
                                     self$region_data$pred_mean_pp_lower    <- sp$lower
                                     self$region_data$pred_mean_pp_upper    <- sp$upper
                                   } else {
                                     self$grid_data   <- set_lin(self$grid_data, "pred_mean_pp",
                                                                 fit$mu_predicted, se_pred$pp)
                                     self$region_data <- set_lin(self$region_data, "pred_mean_total",
                                                                 fit$y_predicted, se_pred$tot)
                                   }
                                 }
                               }
                               
                               if ("rr" %in% type) {
                                 if (is_sampler) {
                                   self$grid_data$rr <- apply(exp(fit$re.samps), 1, mean)
                                   self$grid_data$rr_sd <- apply(exp(fit$re.samps), 1, sd)
                                   self$grid_data$rr_lower <- apply(exp(fit$re.samps), 1, quantile, probs = alpha/2)
                                   self$grid_data$rr_upper <- apply(exp(fit$re.samps), 1, quantile, probs = 1 - alpha/2)
                                 } else {
                                   self$grid_data <- set_log(self$grid_data, "rr", fit$rr, se_pred$rr)
                                 }
                               }
                             }
                             
                             return(invisible(self))
                           },
                           #' @description
                           #' Generate hotspot probabilities
                           #'
                           #' Generate hotspot probabilities. The last model fit will be used to extract
                           #' predictions. If no previous model fit then use either `lgcp_ml()` or `lgcp_bayes()`, or see 
                           #' `model_fit()` to update the stored model fit.
                           #'
                           #' Given a definition of a hotspot in terms of threshold(s) for incidence,
                           #' relative risk, or incidence rate ratio, returns the probabilities
                           #' each area is a "hotspot". See Details of `extract_preds`. Columns
                           #' will be added to `grid_data`. Note that for incidence threshold, the threshold should
                           #' be specified as the per individual incidence.
                           #'
                           #' @param threshold Numeric. Threshold of population standardised incidence
                           #' above which an area is a hotspot
                           #' @param stat Numeric. Threshold of incidence rate ratio
                           #' above which an area is a hotspot.
                           #' @param t.lag integer. Extract predictions for incidence or relative risk for previous time periods.
                           #' @param popdens character string. Name of variable in `grid_data`
                           #' specifying the population density. Needed if `incidence.threshold` is not
                           #' `NULL`
                           #' @return None, called for effects. Columns are added to grid or region data.
                           #' @examples
                           #' \dontrun{
                           #' # See examples for lgcp_bayes() and lgcp_ml()
                           #' }
                           hotspots = function(threshold=NULL,
                                               stat = "rr",
                                               t.lag=0,
                                               popdens = NULL){
                             if((!stat %in% colnames(self$grid_data) & (!is.null(self$region_data) & ! stat %in% colnames(self$region_data)))){
                               self$extract_preds(type = stat, t.lag = t.lag, irr.lag = t.lag, popdens = popdens)
                             }
                             
                             if(stat == "rr"){
                               if(private$last_model_fit$method %in% c("vb","mcmc")){
                                 self$grid_data$hotspot_prob <- 1 - pnorm(log(threshold), log(self$grid_data$rr), self$grid_data$rr_sd)
                               } else {
                                 self$grid_data$hotspot_prob <- 1 - pnorm(log(threshold), log(self$grid_data$rr), self$grid_data$rr_se)
                               }
                             } else if(stat == "irr"){
                               if(private$last_model_fit$region){
                                 if(private$last_model_fit$method %in% c("vb","mcmc")){
                                   self$region_data$hotspot_prob <- 1 - pnorm(log(threshold), log(self$region_data$irr), self$region_data$irr_sd)
                                 } else {
                                   self$region_data$hotspot_prob <- 1 - pnorm(log(threshold), log(self$region_data$irr), self$region_data$irr_se)
                                 }
                               } else {
                                 if(private$last_model_fit$method %in% c("vb","mcmc")){
                                   self$grid_data$hotspot_prob <- 1 - pnorm(log(threshold), log(self$grid_data$irr), self$grid_data$irr_sd)                                
                                 } else {
                                   self$grid_data$hotspot_prob <- 1 - pnorm(log(threshold), log(self$grid_data$irr), self$grid_data$irr_se)                                
                               }
                               }
                             } else if(stat == "pred") {
                               if(private$last_model_fit$region){
                                 if(private$last_model_fit$method %in% c("vb","mcmc")){
                                   self$region_data$hotspot_prob <- 1 - pnorm(threshold, self$region_data$pred_mean_pp, self$region_data$pred_mean_pp_sd)
                                 } else {
                                   self$grid_data$hotspot_prob <- 1 - pnorm(threshold, self$grid_data$pred_mean_pp, self$grid_data$pred_mean_pp_se)
                                 }
                               } else {
                                 if(private$last_model_fit$method %in% c("vb","mcmc")){
                                   self$region_data$hotspot_prob <- 1 - pnorm(threshold, self$grid_data$pred_mean_pp, self$grid_data$pred_mean_pp_sd)
                                 } else {
                                   self$region_data$hotspot_prob <- 1 - pnorm(threshold, self$grid_data$pred_mean_pp, self$grid_data$pred_mean_pp_se)
                                 }
                               }
                             }
                             return(invisible(self))
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
                           #' or population weighted averaging, respectively, or "sum" to take the weighted sum.
                           #' @param popdens character string. If `weight_type` is equal to "pop" then the
                           #' name of the column in `grid_data` with population density data
                           #' @param verbose logical. Whether to provide progress bar.
                           #' @return An `sf` object identical to `new_geom` with additional columns with the
                           #' variables specified in `zcols`
                           #' @examples
                           #' \donttest{
                           #' b1 <- sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
                           #' g1 <- grid$new(b1,0.5)
                           #' dp <- data.frame(y=runif(10,0,3),x=runif(10,0,3),date=paste0("2021-01-",11:20))
                           #' dp <- create_points(dp,pos_vars = c('y','x'),t_var='date')
                           #' cov1 <- grid$new(b1,0.8)
                           #' cov1$grid_data$cov <- runif(nrow(cov1$grid_data))
                           #' g1$add_covariates(cov1$grid_data,
                           #'                   zcols="cov")
                           #' g1$points_to_grid(dp, laglength=5)
                           #' g1$priors <- list(
                           #'   prior_lscale=c(0,0.5),
                           #'   prior_var=c(0,0.5),
                           #'   prior_linpred_mean=c(0),
                           #'   prior_linpred_sd=c(5)
                           #'   )
                           #' res <- g1$lgcp_bayes(popdens="cov", parallel_chains = 1)
                           #' g1$extract_preds(res,
                           #'                  type=c("pred","rr"),
                           #'                  popdens="cov")
                           #' new1 <- g1$aggregate_output(cov1$grid_data,
                           #'                             zcols="rr")
                           #' }
                           aggregate_output = function(new_geom,
                                                       zcols,
                                                       weight_type="area",
                                                       popdens=NULL,
                                                       verbose=TRUE){
                             if(!weight_type %in% c("area","pop","sum"))stop("weight_type must be area, pop, or sum")
                             if(sf::st_crs(self$grid_data)!=sf::st_crs(new_geom)){
                               sf::st_crs(self$grid_data) <- sf::st_crs(new_geom)
                               warning("st_crs(self$grid_data)!=st_crs(new_geom) setting equal")
                             }
                             
                             vals <- matrix(nrow=nrow(new_geom),ncol=length(zcols))
                             if(verbose)cat("Overlaying geographies\n")
                             grid_area <- as.numeric(sf::st_area(self$grid_data[1,]))
                             tmp <- suppressWarnings(sf::st_intersection(new_geom,self$grid_data[,c(popdens,zcols)]))
                             rnames <- rownames(tmp)
                             rnames <- gsub("\\..*","",rnames)
                             gnames <- rownames(new_geom)
                             
                             for(i in 1:nrow(new_geom)){
                               tmp_range <- tmp[which(rnames == gnames[i]),]
                               w <- as.numeric(sf::st_area(tmp_range))/grid_area
                               for(j in 1:length(zcols)){
                                 if(weight_type == "sum"){
                                   vals[i,j] <- sum(as.data.frame(tmp_range)[,zcols[j]][w > 0.5])
                                 } else if(weight_type == "area"){
                                   vals[i,j] <- weighted.mean(as.data.frame(tmp_range)[,zcols[j]],w =w)
                                 } else {
                                   w <- w*as.data.frame(tmp_range)[,popdens]
                                   vals[i,j] <- weighted.mean(as.data.frame(tmp_range)[,zcols[j]],w =w)
                                 }
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
                           #' Coordinates are scaled to `[-1,1]` for LGCP models fit with HSGP. This function
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
                             std_val <- max(diff(xrange),diff(yrange))
                             return(std_val)
                           },
                           #' @description
                          #' Builds the mesh and relevant matrices for SPDE estimation
                          #' 
                          #' @param max_edge The largest allowed triangle edge length. One or two values.
                          #' @param cutoff The minimum allowed distance between points. Point at most as far apart as this are replaced by a single vertex prior to the mesh refinement step.
                          #' @param offset The automatic extension distance. One or two values, for an inner and an optional outer extension. If negative, interpreted as a factor relative to the approximate data diameter 
                            #' @return A list with A, C, G, and W matrices.
                           build_spde_data = function(max_edge = NULL, cutoff = NULL, offset = NULL) {
                             
                             is_aggregated <- !is.null(self$region_data)
                             
                             # ── Mesh seed locations and boundary ───────────────────────────────────────
                             if (is_aggregated) {
                               mesh_loc <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(self$region_data)))
                               bb       <- sf::st_bbox(self$region_data)
                             } else {
                               mesh_loc <- as.matrix(unique(sf::st_coordinates(self$point_data)))
                               bb <- c(xmin = min(mesh_loc[,1]), xmax = max(mesh_loc[,1]),
                                       ymin = min(mesh_loc[,2]), ymax = max(mesh_loc[,2]))
                             }
                             
                             span     <- min(diff(range(mesh_loc[,1])), diff(range(mesh_loc[,2])))   # ~20000
                             if (span < 10) {
                               warning("Domain span (", round(span, 2), ") is small. SPDE precision matrix ",
                                       "may be ill-conditioned. Consider rescaling coordinates by a factor of ",
                                       round(50 / span), " using self$transform(scale_factor = a) for numerical stability.")
                             }
                             # ── Default mesh args ──────────────────────────────────────────────────────
                             if (is.null(max_edge)) {
                               h        <- span / 18                                                    # ~1000 m
                               max_edge <- c(h, 1.5 * h)                                                  # c(1000, 5000)
                             }
                             if (is.null(cutoff)) cutoff   <- h  
                             if (is.null(offset)) offset   <- c(h, 2 * h)
                             
                             mesh <- fmesher::fm_mesh_2d(
                               loc      = mesh_loc,
                               boundary = self$boundary,
                               max.edge = max_edge,
                               cutoff   = cutoff,
                               offset   = offset
                             )
                             fem    <- fmesher::fm_fem(mesh)
                             C_diag <- Matrix::diag(fem$c0)
                             G      <- Matrix::drop0(as(as(fem$g1, "generalMatrix"), "CsparseMatrix"))
                             n_v    <- nrow(mesh$loc)
                             
                             # ── A_loc and W_mesh: branch on aggregation ────────────────────────────────
                             W_mesh <- NULL
                             if (is_aggregated) {
                               # Aggregated: model evaluates at vertices, A_loc = I(n_v)
                               A_loc <- Matrix::sparseMatrix(i = 1:n_v, j = 1:n_v, x = 1,
                                                             dims = c(n_v, n_v))
                               
                               vert_sf <- sf::st_as_sf(
                                 data.frame(x = mesh$loc[,1], y = mesh$loc[,2]),
                                 coords = c("x","y"), crs = sf::st_crs(self$region_data)
                               )
                               rid  <- sf::st_intersects(vert_sf, self$region_data)
                               rid  <- vapply(rid, function(z) if (length(z)) z[1] else NA_integer_, integer(1))
                               keep <- !is.na(rid)
                               W_mesh <- Matrix::sparseMatrix(
                                 i = rid[keep], j = which(keep), x = C_diag[keep],
                                 dims = c(nrow(self$region_data), n_v)
                               )
                             } else {
                               # LGCP: augmented rows = events (barycentric) ++ vertices (identity)
                               aug_loc <- rbind(as.matrix(sf::st_coordinates(self$point_data)), mesh$loc[, 1:2])
                               aug_loc <- unique(aug_loc)
                               A_loc   <- Matrix::drop0(fmesher::fm_basis(mesh, loc = aug_loc))
                             }
                             
                             # ── A_pred: prediction projector to grid centroids ─────────────────────────
                             A_pred <- NULL
                             if (!is.null(self$grid_data)) {
                               grid_xy <- sf::st_coordinates(sf::st_centroid(self$grid_data))
                               A_pred <- Matrix::drop0(fmesher::fm_basis(mesh, loc = grid_xy))
                             }
                             print(c(sum_C = sum(C_diag), domain_area = sum(sf::st_area(self$boundary))))
                             list(mesh = mesh, A_loc = A_loc, C = C_diag, G = G,
                                  W_mesh = W_mesh, A_pred = A_pred, n_v = n_v,
                                  is_aggregated = is_aggregated)
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
                               q_weights = private$intersection_data$w,
                               W = Matrix::sparseMatrix(j = private$intersection_data$grid_id, p = ncell-1, x = as.vector(private$intersection_data$w))
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
                           #'  minimum distance to the previous observations, or "random" which randomly orders them.
                           #' @param verbose Logical indicating whether to print a progress bar (TRUE) or not (FALSE).
                           #' @return No return, used for effects.
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
                           #' @param approx Either "rank" for reduced rank approximation, or "nngp" for nearest 
                           #' neighbour Gaussian process. 
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
                           #' Returns the random effects stored in the object (if any) after using ML fitting. It's main use is
                           #' if a fitting procedure is stopped, the random effects can still be returned.
                           #' @return A matrix of random effects samples if a MCMCML model has been initialised, otherwise returns FALSE
                           get_random_effects = function(){
                             if(!is.null(private$ptr)){
                               u <- rtsModel__u(private$ptr,private$cov_type,private$lp_type)
                               return(u)
                             } else {
                               return(FALSE)
                             }
                           },
                           #' @description
                           #' Either returns the stored last model fit with either `lgcp_ml` or `lgcp_bayes`, or updates 
                           #' the saved model fit if an object is provided.
                           #' @param fit Optional. A previous `rtsFit` object. If provided then the function updates the internally stored model fit.
                           #' @return Either a `rtsFit` object or nothing if no model has been previously fit, or if the fit is updated.
                           model_fit = function(fit = NULL){
                             if(!is.null(fit) & !is(fit,"rtsFit"))stop("fit must be an rtsFit")
                             if(is.null(fit)){
                               if(!is.null(private$last_model_fit)) {
                                 return(private$last_model_fit)
                               } else {
                                 stop("No stored model fit")
                               }
                             } else {
                               private$last_model_fit <- fit
                             }
                           }
                         ),
                    private = list(
                      intersection_data = NULL,
                      grid_ptr = NULL,
                      region_ptr = NULL,
                      cov_type = 1,
                      lp_type = 1,
                      last_model_fit = NULL,
                      grid_adj_matrix = NULL,
                      stack_variable = function(var, use_grid = FALSE){
                        if(is.null(self$region_data)){
                          nT <- sum(grepl("\\bt[0-9]",colnames(self$grid_data)))
                        } else {
                          nT <- sum(grepl("\\bt[0-9]",colnames(self$region_data)))
                        }
                        nT <- ifelse(nT==0,1,nT)
                        if(is.null(self$region_data) | use_grid){
                          #population density
                          nColP <- sum(grepl(var,colnames(self$grid_data)))
                          if(nColP==1){
                            popd <- rep(as.data.frame(self$grid_data)[,var],nT)
                          } else if(nColP==0){
                            stop("variable not found in grid data")
                          } else {
                            if(nT>1){
                              popd <- stack(as.data.frame(self$grid_data)[,paste0(var,1:nT)])[,1]
                            } else {
                              popd <- as.data.frame(self$grid_data)[,var]
                            }
                          }
                        } else {
                          #population density
                          nColP <- sum(grepl(var,colnames(self$region_data)))
                          if(nColP==1){
                            popd <- rep(as.data.frame(self$region_data)[,var],nT)
                          } else if(nColP==0){
                            stop("variable not found in region data")
                          } else {
                            if(nT>1){
                              popd <- stack(as.data.frame(self$region_data)[,paste0(var,1:nT)])[,1]
                            } else {
                              popd <- as.data.frame(self$region_data)[,var]
                            }
                          }
                        }
                        return(popd)
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
                          } else {
                            if("y"%in%colnames(self$grid_data))stop("both t and y case count variables found")
                          }
                        } else {
                          nT <- sum(grepl("\\bt[0-9]",colnames(self$region_data)))
                          if(nT==0){
                            if("y"%in%colnames(self$region_data)){
                              nT <- 1
                            } else {
                              stop("case counts not defined in data")
                            }
                          }  else {
                            if("y"%in%colnames(self$grid_data))stop("both t and y case count variables found")
                          }
                        }
                        nCell <- nrow(self$grid_data)
                        x_grid <- as.data.frame(suppressWarnings(sf::st_coordinates(sf::st_centroid(self$grid_data))))
                        if(approx=="hsgp"){
                          # scale to -1,1 in all dimensions
                          xrange <- range(x_grid[,1])
                          yrange <- range(x_grid[,2])
                          scale_f <- max(diff(xrange), diff(yrange))
                          x_grid[,1] <- -1 + 2*(x_grid[,1] - xrange[1])/scale_f
                          x_grid[,2] <- -1 + 2*(x_grid[,2] - yrange[1])/scale_f
                          L_boundary <- c(L,L)
                        }
                        
                        #add covariates
                        if(!is.null(covs)){
                          nQ <- length(covs)
                          X <- matrix(NA,nrow=nrow(self$grid_data)*nT,ncol=nQ+1)
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
                          X <- matrix(1,nrow=nrow(self$grid_data)*nT,ncol=1)
                          Q <- 1
                        }
                        
                        if(is.null(self$region_data)){
                          # outcome data
                          if(nT > 1){
                            y <- stack(as.data.frame(self$grid_data)[,paste0("t",1:nT)])[,1]
                          } else {
                            y <- as.data.frame(self$grid_data)[,"y"]
                          }
                          popd <- private$stack_variable(popdens)
                          
                        } else {
                          #outcomes
                          if(nT > 1){
                            y <- stack(as.data.frame(self$region_data)[,paste0("t",1:nT)])[,1]
                          } else {
                            y <- as.data.frame(self$region_data)[,"y"]
                          }
                          if(any(grepl(popdens,colnames(self$grid_data))) | any(grepl(popdens,colnames(self$region_data)))){
                            popd <- private$stack_variable(popdens, any(grepl(popdens,colnames(self$grid_data))))
                          } else {
                            popd <- rep(1, nT*nrow(self$grid_data))
                          }
                          nG <- 0
                          X_g <- matrix(0,nrow=nrow(self$grid_data)*nT,ncol=1)
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
                        if(verbose)message(paste0(nCell," grid cells ",nT," time periods, and ",Q," covariates."))
                        
                        datlist <- list(
                          D = 2,
                          Q = Q,
                          Nsample = nCell,
                          nT= nT,
                          y = y,
                          x_grid = x_grid[,1:2],
                          popdens = popd,
                          X= X,
                          mod = mod,
                          approx = 0,
                          is_region = 0,
                          n_region = 1,
                          n_Q = 1,
                          n_cell = as.array(1),
                          cell_id = as.array(1),
                          q_weights = as.array(1),
                          Q_g = 0,
                          X_g = matrix(1,1,1),
                          M = m,
                          M_nD = m^2,
                          indices = matrix(1,1,1),
                          L = c(0,0),
                          NN = matrix(1,1,1)
                        )
                        
                        if(approx == "hsgp"){
                          datlist$approx <- 1
                          datlist$L <- L_boundary
                          datlist$indices <- ind
                          
                        } else if(approx == "nngp"){
                          datlist$approx <- 2
                          gptr <- GridData__new(as.matrix(x_grid[,1:2]), nT)
                          GridData__gen_NN(gptr, m)
                          datlist$NN <- GridData__NN(gptr)+1
                        } 
                        
                        if(!is.null(self$region_data)){
                          ncell <- unname(table(private$intersection_data$region_id))
                          ncell <- c(1, cumsum(ncell)+1)
                          datlist$is_region <- 1
                          datlist$Q_g <- nG
                          datlist$X_g <- X_g
                          datlist$n_region <- nrow(self$region_data)
                          datlist$n_Q <- nrow(private$intersection_data)
                          datlist$n_cell <- ncell
                          datlist$cell_id <- private$intersection_data$grid_id
                          datlist$q_weights <- private$intersection_data$w
                        } 
                        
                        if(bayes){
                          datlist <- append(datlist,list(prior_lscale=self$priors$prior_lscale,
                                                         prior_var=self$priors$prior_var,
                                                         prior_linpred_mean = as.array(self$priors$prior_linpred_mean),
                                                         prior_linpred_sd=as.array(self$priors$prior_linpred_sd)))
                        }
                        return(datlist)
                      },
                      #' @description 
                      #' Interpolates covariate values at coordinates, used typically for SPDE models
                      #' @param coords n x 2 matrix of prediction coordinates. 
                      #' @param covs Vector of strings naming columns in grid data
                      #' @return Matrix of interpolated covariate values.
                      evaluate_covariates_at = function(coords, covs) {
                        
                        if (is.null(covs) || length(covs) == 0) {
                          return(matrix(0, nrow(coords), 0))
                        }
                        
                        grid_cols <- colnames(self$grid_data)
                        missing   <- setdiff(covs, grid_cols)
                        if (length(missing) > 0) {
                          stop("Covariate(s) not found in grid_data: ", paste(missing, collapse = ", "))
                        }
                        
                        query <- if (is.data.frame(coords)) {
                          as.matrix(coords[, c("x", "y")])
                        } else {
                          as.matrix(coords)[, 1:2, drop = FALSE]
                        }
                        
                        grid_xy <- if (inherits(self$grid_data, "sf")) {
                          sf::st_coordinates(sf::st_centroid(sf::st_geometry(self$grid_data)))
                        } else {
                          as.matrix(self$grid_data[, c("x", "y")])
                        }
                        
                        nn_idx <- if (requireNamespace("nabor", quietly = TRUE)) {
                          nabor::knn(grid_xy, query, k = 1)$nn.idx[, 1]
                        } else if (requireNamespace("FNN", quietly = TRUE)) {
                          FNN::get.knnx(grid_xy, query, k = 1)$nn.index[, 1]
                        } else {
                          apply(query, 1, function(p) {
                            which.min((grid_xy[, 1] - p[1])^2 + (grid_xy[, 2] - p[2])^2)
                          })
                        }
                        
                        # Drop sf geometry before subsetting so we get a plain data frame
                        gd  <- if (inherits(self$grid_data, "sf")) sf::st_drop_geometry(self$grid_data) else self$grid_data
                        out <- as.matrix(gd[nn_idx, covs, drop = FALSE])
                        colnames(out) <- covs
                        out
                      }
                    ))

