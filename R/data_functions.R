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
#' @param t_var character string with the name of the column with the date of the case. If single-period
#' analysis then set t_var to NULL.
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
                          t_var = NULL,
                          format="%Y-%m-%d",
                          verbose=TRUE){

  if(verbose)message(paste0("Using ",pos_vars[1]," as y-coordinate and ",pos_vars[2]," as x-coordinate."))
  if(!is(data,"data.frame"))stop("data not data.frame")
  if(any(!pos_vars%in%colnames(data)))stop("pos_vars not in colnames(data)")

  if(!is.null(t_var))
  {
    if(!t_var%in%colnames(data))stop("t_var not in colnames(data)")
    if(any(is.na(data[,t_var]))|any(is.na(as.Date(data[,t_var], format=format))))
      stop(paste0(t_var," does not contain date in ",format," format"))
  }

  out <- lapply(1:nrow(data),function(i)sf::st_point(c(data[i,pos_vars[2]],data[i,pos_vars[1]])))
  dp <- sf::st_sfc(out)

  if(!is.null(t_var)){
    dp <- sf::st_sf(dp, data.frame(t=data[,t_var]))
    dp$t <- as.Date(dp$t, format=format)
  } else {
    dp <- sf::st_sf(dp)
  }
  
  return(dp)
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

#' Simulated point data for running single-period examples
#' 
#' A set of 261 points simulated within the boundary of the city Birmingham, UK 
#' from a log-Gaussian Cox process.
"example_points"

#' Boundary polygon for Birmingham, UK
#' 
#' A Boundary polygon describing the border of the city of Birmingham, UK.
"boundary"

#' Birmingham crime data
#' 
#' Counts of burglaries for the months of 2022 for the city of Birmingham, UK at the 
#' Middle-Layer Super Output Area.
"birmingham_crime"

