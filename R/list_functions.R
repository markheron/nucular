##' function variations to run on lists of vectors, handling them as one vector
##' 
##' Copies of several functions that are originally for numeric vectors, but in this _list variations work on list of vectors
##' 
##' @name list_functions
##' @author Mark Heron
##' 
NULL


##' mean_list
##'
##' calculates the mean of all list elements together
##' @export
##' @param x list of numeric vectors
##' @return mean of the list
##' @author Mark Heron
mean_list <- function(x, na.rm=TRUE) {
  
  x_sum <- sum(unlist(lapply(x, sum, na.rm=na.rm)), na.rm=na.rm)
  x_len <- sum(unlist(lapply(x, function(tmp) length(na.omit(tmp)))), na.rm=na.rm)
  
  return(x_sum/x_len)
}



##' center_list
##'
##' centers all list elements by their combined mean
##' @export
##' @param x list of numeric vectors
##' @return centered x (mean = 0)
##' @author Mark Heron
center_list <- function(x) {
  return( lapply(x, function(a) a - mean_list(x)) )
}



##' scale_list
##'
##' scales all list elements by their combined mean
##' @export
##' @param x list of numeric vectors
##' @return scaled x (mean = 1)
##' @author Mark Heron
scale_list <- function(x) {
  return( lapply(x, function(a) a / mean_list(x)) )
}


