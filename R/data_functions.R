#' Create sf grid over a boundary
#'
#' Produces a regular grid over an area of interest
#'
#' Given a contiguous boundary describing an area of interest, returns an sf
#' object of a regular grid within the limits of the boundary.
#'
#' @param boundary An sf object containing one polygon describing the area of interest
#' @param cellsize The dimension of the grid cells
#' @return A sf object describing a set of square cells tiled over the area of the boundary
#' @examples
#' b1 = sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
#' create_grid(b1,0.5)
#' @export
create_grid <- function(boundary,
                        cellsize){

  if(!is(boundary,"sf"))stop("boundary not sf")
  if(nrow(boundary)!=1)stop("boundary should only contain one polygon")
  if(!is(cellsize,"numeric"))stop("cellsize not numeric")

  bgrid <- sf::st_make_grid(boundary,cellsize = cellsize)
  bgrid <- sf::st_sf(bgrid)
  idx1<- sf::st_contains(y=bgrid,x=boundary)
  idx2<- sf::st_intersects(y=bgrid,x=boundary)
  idx <- c(unlist( sapply( idx1, `[`) ),unlist( sapply( idx2, `[`) ))
  bgrid <- bgrid[idx,]

  #class(bgrid) <- c(class(bgrid),"rts_grid")

  return(bgrid)
}

#' Create sf object from point location data
#'
#' Produces an sf object with location and time of cases from a data frame
#'
#' Given a data frame containing the point location and date of cases, the function
#' will return an sf object of the points with the date information.
#' @param data data.frame with the x- and y-coordinate of case locations and the date
#' of the case.
#' @param pos_vars vector of length two with the names of the columns
#' containing the y and x coordinates, respectively.
#' @param t_var character string with the name of the column with the date of the case
#' @param format character string with the format of the date specified by t_var. See
#' \link[base]{strptime}
#' @param verbose Logical indicating whether to print information
#' @return An sf object of the same size as `data`
#' @examples
#' dp <- data.frame(y=runif(10,0,3),x=runif(10,0,3),date=paste0("2021-01-",11:20))
#' dp <- create_points(dp,pos_vars = c('y','x'),t_var='date')
#' @export
create_points <- function(data,
                          pos_vars = c('lat','long'),
                          t_var,
                          format="%Y-%m-%d",
                          verbose=TRUE){

  if(verbose)message(paste0("Using ",pos_vars[1]," as y-coordinate and ",pos_vars[2]," as x-coordinate."))
  if(!is(data,"data.frame"))stop("data not data.frame")
  if(any(!pos_vars%in%colnames(data)))stop("pos_vars not in colnames(data)")
  if(!t_var%in%colnames(data))stop("t_var not in colnames(data)")
  if(any(is.na(data[,t_var]))|any(is.na(as.Date(data[,t_var], format=format))))
    stop(paste0(t_var," does not contain date in ",format," format"))

  out <- lapply(1:nrow(data),function(i)sf::st_point(c(data[i,pos_vars[2]],data[i,pos_vars[1]])))
  dp <- sf::st_sfc(out)
  dp <- sf::st_sf(dp, data.frame(t=data[,t_var]))
  dp$t <- as.Date(dp$t, format=format)
  return(dp)
}

#' Generates case counts of points over a grid
#'
#' Counts the number of cases in each time period in each grid cell
#'
#' Given the sf object describing a regular grid output from `create_grid()`
#' and the sf object with the point locations and date output from
#' `create_points()`, returns the grid object with additional columns indicating
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
#' @param grid_data sf object describing a regular grid over the area of interest.
#' See \link[rts2]{create_grid}
#' @param point_data sf object describing the point location of cases with a column
#' `t` of the date of the case in YYYY-MM-DD format. See \link[rts2]{create_points}
#' @param t_win character string. One of "day", "week", or "month" indicating the
#' length of the time windows in which to count cases
#' @param laglength integer The number of time periods to include counting back from the most
#' recent time period
#' @return An sf object identical to `grid_data` but with additional columns with the
#' case count and date of each period
#' @examples
#' b1 <- sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
#' g1 <- create_grid(b1,0.5)
#' dp <- data.frame(y=runif(10,0,3),x=runif(10,0,3),date=paste0("2021-01-",11:20))
#' dp <- create_points(dp,pos_vars = c('y','x'),t_var='date')
#' points_to_grid(g1, dp, laglength=5)
#' @export
points_to_grid <- function(grid_data,
                           point_data,
                           t_win = c("day"),
                           laglength = 14){

  if(!is(grid_data,"sf"))stop("grid not sf")
  if(!is(point_data,"sf"))stop("points not sf")
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

  if(sf::st_crs(point_data)!=sf::st_crs(grid_data)){
    warning("CRS not equal. Setting st_crs(point_data)==st_crs(grid_data)")
    sf::st_crs(point_data) <- sf::st_crs(grid_data)
  }

  for(i in 1:length(tuniq)){
    grid_data$y <-  lengths(sf::st_intersects(grid_data,
                                              point_data[tdat==tuniq[i],]))
    colnames(grid_data)[length(colnames(grid_data))] <- paste0("t",i)
    grid_data$d <- min(point_data[tdat==tuniq[i],]$t)
    colnames(grid_data)[length(colnames(grid_data))] <- paste0("date",i)
  }

  return(grid_data)
}

#' Generate day of week data
#'
#' Create data frame with day of week indicators
#'
#' Takes the output of `points_to_grid()` and generates a data frame with indicator
#' variables for each day of the week for use in the `lgcp_fit()` function.
#'@param grid_data sf object with columns `t[0-9]` and `date[0-9]`. See \link[rts2]{points_to_grid}
#'@return data.frame with columns `t`, `day`, and `dayMon` to `daySun`
#'@examples
#' b1 <- sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
#' g1 <- create_grid(b1,0.5)
#' dp <- data.frame(y=runif(10,0,3),x=runif(10,0,3),date=paste0("2021-01-",11:20))
#' dp <- create_points(dp,pos_vars = c('y','x'),t_var='date')
#' g1 <- points_to_grid(g1, dp, laglength=5)
#' get_dow(g1)
#' @importFrom stats model.matrix
#' @export
get_dow <- function(grid_data){
  if(!is(grid_data,"sf"))stop("grid_data not sf")

  nT <- length(colnames(grid_data)[grepl("t[0-9]",colnames(grid_data))])
  dw <- data.frame(t=1:nT,day=NA)
  for(i in 1:nT){
    dw$day[i] <- as.character(lubridate::wday(as.data.frame(grid_data)[1,paste0("date",i)],
                                   label = TRUE))
  }
  dx <- model.matrix(~day-1,data=dw)
  dw <- cbind(dw,as.data.frame(dx))

  return(dw)
}

#' Adds covariate data to sf object
#'
#' Maps spatial, temporal, or spatio-temporal covariate data onto sf object
#'
#' @details
#' *Spatially-varying data only* `cov_data`` is an sf object describing covariate
#' values for a set of polygons over the area of interest. The values are mapped
#' onto another sf object - the grid `grid_data` covering the same area (see
#' \link[rts2]{create_grid}). For each grid cell in `grid_data` a weighted
#' average of each covariate listed in `zcols` is generated with weights either
#' equal to the area of intersection of the grid cell and the polygons in
#' `cov_data` (`weight_type="area"`), or this area multiplied by the population
#' density of the polygon for population weighted (`weight_type="pop"`). Columns
#' with the names in `zcols` are added to the output.
#'
#' *Temporally-varying only data* `cov_data` is a data frame with number of rows
#' equal to the number of time periods. One of the columns must be called `t` and
#' have values from 1 to the number of time periods. The other columns of the data
#' frame have the values of the covariates for each time period. See
#' \link[rts2]{get_dow} for day of week data. A total of
#' length(zcols)*(number of time periods) columns are added to the output: for each
#' covariate there will be columns appended with each time period number. For example,
#' `dayMon1`, `dayMon2`, etc.
#'
#' *Spatially and temporally varying data* There are two ways to add data that
#' vary both spatially and temporally. The final output for use in analysis must
#' have a column for each covariate and each time period with the same name appended
#' by the time period number, e.g. `covariateA1`,`covariateA2`,... If the covariate
#' values for different time periods are in separate sf objects, one can follow
#' the method for spatially-varying only data above and append the time period number
#' using the argument `t_label`. If the values for different time periods are in the same
#' sf object then they should be named as described above and then can be added
#' as for spatially-varying covariates, e.g. `zcols=c("covariateA1","covariateA2")`.
#'
#' @param grid_data sf object describing a regular grid over the area of interest.
#' See \link[rts2]{create_grid}
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
#' @return sf object identical to grid_data but with additional columns added with covariate data. See details.
#' @examples
#' b1 <-  sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
#' g1 <- create_grid(b1,0.5)
#' cov1 <- create_grid(b1,0.8)
#' cov1$cov <- runif(nrow(cov1))
#' g1 <- add_covariates(g1,
#'                      cov1,
#'                      zcols="cov",
#'                     verbose = FALSE)
#' @importFrom stats weighted.mean
#' @export
add_covariates <- function(grid_data,
                           cov_data,
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
    sf::st_agr(grid_data) = "constant"
    tmp <- sf::st_intersection(cov_data[,zcols],grid_data)
    tmp_len <- lengths(sf::st_intersects(grid_data,cov_data))
    tmp_len <- 1 - tmp_len[1] + cumsum(tmp_len)
    vals <- matrix(nrow=nrow(grid_data),ncol=length(zcols))

    if(verbose)cat("Overlaying geographies\n")

    for(i in 1:nrow(grid_data)){
      if(i < nrow(grid_data)){
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

      if(verbose)cat("\r",progress_bar(i,nrow(grid_data)))
    }

    for(j in 1:length(zcols)){
      grid_data$x <- vals[,j]
      if(is.null(t_label)){
        colnames(grid_data)[length(colnames(grid_data))] <- zcols[j]
      } else {
        colnames(grid_data)[length(colnames(grid_data))] <- paste0(zcols[j],t_label)
      }

    }
  } else {
    nT <- max(cov_data$t)
    for(j in zcols){
      for(t in 1:nT){
        grid_data$x <- cov_data[cov_data$t==t,j]
        colnames(grid_data)[length(colnames(grid_data))] <- paste0(j,t)
      }
    }
  }

  return(grid_data)

}

#' Generates a progress bar
#'
#' Prints a progress bar
#'
#' @param i integer. The current iteration.
#' @param n integer. The total number of interations
#' @param len integer. Length of the progress a number of characters
#' @return A character string
#' @examples
#' progress_bar(10,100)
#' @export
progress_bar <- function(i,n,len=30){
  prop <- floor((i*100/n) / (100/len))
  pt1 <- paste0(rep("=",prop), collapse="")
  pt2 <- paste0(rep(" ",len-prop), collapse = "")
  msg <- paste0("|",pt1,pt2,"| ",round((i*100/n),0),"%")
  return(msg)
}
