##' function variations to run on lists of vectors, handling them as one vector
##' 
##' Copies of several functions that are originally for numeric vectors, but in this _list variations work on list of vectors
##' 
##' @name list_functions
##' @author Mark Heron
##' 
NULL


##' sum_list
##'
##' Calculates the sum of all list elements together.
##' 
##' @export
##' @param x list of numeric vectors
##' @param na.rm should NA values be ignored (see \code{\link{sum}})
##' @return (numeric) sum of all values in the list
##' 
sum_list <- function(x, na.rm=TRUE) {
  
  sum <- sum(unlist(lapply(x, sum, na.rm=na.rm)), na.rm=na.rm)
  
  return(sum)
}

##' mean_list
##'
##' Calculates the mean of all list elements together.
##' I.e. the mean of all values contained in list elements, \emph{not} the mean of the mean of all list elements!
##' 
##' @export
##' @param x list of numeric vectors
##' @param na.rm should NA values be ignored (see \code{\link{sum}})
##' @return (numeric) mean of all values in the list
##'
mean_list <- function(x, na.rm=TRUE) {
  
  x_sum <- sum_list(x, na.rm=na.rm)
  x_len <- sum(unlist(lapply(x, function(tmp) length(if(na.rm) {na.omit(tmp)} else {tmp} ))), na.rm=na.rm)
  
  return(x_sum/x_len)
}



##' center_list
##'
##' Centers all list elements by their combined mean.
##' Works with lists of ff vectors.
##' 
##' @export
##' @param x list of numeric vectors (or ff vectors)
##' @return (list of numeric (ff) vectors) centered \code{x} (i.e., mean_list == 0)
##' 
center_list <- function(x) {
  tmp_list_mean <- mean_list(x)
  return( lapply(x, function(a) a - tmp_list_mean) )
}



##' scale_list
##'
##' Scales all list elements by their combined mean.
##' Works with lists of ff vectors.
##' 
##' @export
##' @param x list of numeric vectors (or ff vectors)
##' @return (list of numeric (ff) vectors) scaled \code{x} (i.e., mean_list == 1)
##' 
scale_list <- function(x) {
  tmp_list_mean <- mean_list(x)
  return( lapply(x, function(a) a / tmp_list_mean) )
}


##' intra_scaling_list
##'
##' Scales each vector in the list to mean == 1.
##' Simple wrapper for lapply '/' mean.
##' Works with lists of ff vectors.
##' 
##' @export
##' @param list where each vector should be scaled individually
##' @return intra scaled list
##' 
intra_scaling_list <- function(list) {
  
  return(lapply(list, function (x) x / mean(x, na.rm=TRUE)))
}

