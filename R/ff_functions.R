##' function variations to run on ff vectors
##' 
##' Copies of several functions that are originally for numeric vectors, but in this *_ff variations work on ff vectors
##' 
##' @name ff_functions
##' @author Mark Heron
##' 
##' @import ff
##' @import ffbase
NULL



##' cov_ff
##' 
##' Computest the covariance between two ff vectors.
##' @export
##' @param x a ff vector
##' @param y a second ff vector
##' @param centered boolean, is the mean already = 0 ?
##' @return covariance between x and y
##' @author Mark Heron
cov_ff <- function(x,y,centered=FALSE) {
  
  if(centered) {
    return( sum( x*y , na.rm=TRUE) )
  } else {
    return( sum( (x- mean(x, na.rm=TRUE))*(y-mean(y, na.rm=TRUE)) ,na.rm=TRUE) )
  }
}



##' sd_ff
##' 
##' Computes the standard deviation of a vector.
##' @export
##' @param x a ff vector
##' @param centered boolean, is the mean already = 0 ?
##' @return standard deviation of x
##' @author Mark Heron
sd_ff <- function(x, centered=FALSE) {
  
  if(centered) {
    return( sqrt(sum( x^2 , na.rm=TRUE) ) )
  } else {
    return( sqrt(sum( (x-mean(x, na.rm=TRUE))^2 ,na.rm=TRUE)) )
  }  
}



##' cor_ff
##' 
##' Computes the correlation between two vectors.
##' @export
##' @param x a ff vector
##' @param y a second ff vector
##' @return pearson correlation between x and y
##' @author Mark Heron
cor_ff <- function(x,y) {
  
  centered_x = x - mean(x, na.rm=TRUE)
  centered_y = y - mean(y, na.rm=TRUE)
  return( cov_ff(centered_x, centered_y, centered=TRUE) / (sd_ff(centered_x, centered=TRUE)*sd_ff(centered_y, centered=TRUE)) )  
}
