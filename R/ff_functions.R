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
    return( sum( x*y , na.rm=TRUE) / (sum( !is.na(x) & !is.na(y) )-1) )
  } else {
    return( cov_ff(x- mean(x, na.rm=TRUE)),(y-mean(y, na.rm=TRUE)), centered=TRUE )
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
    return( sqrt(sum( x^2 , na.rm=TRUE) / (sum(!is.na(x))-1) ) )
  } else {
    return( sd_ff(x-mean(x, na.rm=TRUE), centered=TRUE) )
  }  
}



##' cor_ff
##' 
##' Computes the correlation between two vectors.
##' @export
##' @param x a ff vector
##' @param y a second ff vector
##' @param method ("pearson") which correlation should be computed, either "pearson" or "spearman" (has to briefly transform the ff's to vectors for rank)
##' @return pearson correlation between x and y
##' @author Mark Heron
cor_ff <- function(x,y, method="pearson") {
  
  if(method == "pearson") {
    centered_x = x - mean(x, na.rm=TRUE)
    centered_y = y - mean(y, na.rm=TRUE)
  } else if(method == "spearman") {
    ranked_x = as.ff(rank(as.ram(x), na.last="keep")) # as.ram has to be used or rank fails!
    ranked_y = as.ff(rank(as.ram(y), na.last="keep"))
    return(cor_ff(ranked_x, ranked_y, method="pearson"))
  } else {
    stop("only 'pearson' and 'spearman' correlations implemented in cor_ff")
  }
  return( cov_ff(centered_x, centered_y, centered=TRUE) / (sd_ff(centered_x, centered=TRUE)*sd_ff(centered_y, centered=TRUE)) )  
}
