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



##' Computes a smeared vector, basically a running_sum.
##' 
##' Add's the vector on its self several times, so good for short smears of long vectors.
##' 
##' @export
##' @title smear_ff
##' @param x ff vector to smear
##' @param from where to start
##' @param to where to go
##' @return smeared ff vector
##' @author Mark Heron
smear_ff <- function(x, from=0, to=1) {
  # the from, to parameter meaning seem reversed compared to my intuition (at least at the moment)
  # so maybe I should change it ...
  
  if(from == 0 & to == 0) {
    return(x)
  }
  
  smooth_len <- to-from+1
  times <- floor(log2(smooth_len))-1
  
  smeared <- c(x[], rep(0,to)) #to
  len <- length(smeared)
  for(iter in 0:times) {
    i <- 2^iter
    smeared <- smeared + c(rep(0,i),smeared[1:(len-i)])
  }
  if((2^(times+1)) < smooth_len) {
    len <- length(x)
    for(i in (-to-1+((2^(times+1)+1):smooth_len)) ) {
      smeared <- smeared + c(rep(0,to-min(-i,0)),x[max(-i+1,1):min(len-i,len)],rep(0,max(-i,0)))
    }
  }
  if(to == 0) {
    return(as.ff(smeared))
  } else {
    return(as.ff(smeared[-(1:to)])) #from
  }
}



##' Computest the covariance between two ff vectors.
##'
##' @export
##' @title cov_ff
##' @param x a ff vector
##' @param y a second ff vector
##' @param meaned boolean, is the mean already = 0 ?
##' @return covariance between x and y
##' @author Mark Heron
cov_ff <- function(x,y,meaned=FALSE) {
  
  if(meaned) {
    return( sum( x*y ) )
  } else {
    return( sum( (x- mean(x, na.rm=TRUE))*(y-mean(y, na.rm=TRUE)) ) )
  }
}



##' Computes the standard deviation of a vector.
##'
##' @export
##' @title sd_ff
##' @param x a ff vector
##' @param meaned boolean, is the mean already = 0 ?
##' @return standard deviation of x
##' @author Mark Heron
sd_ff <- function(x, meaned=FALSE) {
  
  if(meaned) {
    return( sqrt(sum( x^2 ) ) )
  } else {
    return( sqrt(sum( (x-mean(x, na.rm=TRUE))^2 )) )
  }  
}



##' Computes the correlation between two vectors.
##'
##' @export
##' @title cor_ff
##' @param x a ff vector
##' @param ya second ff vector
##' @return pearson correlation between x and y
##' @author Mark Heron
cor_ff <- function(x,y) {
  
  meaned_x = x - mean(x, na.rm=TRUE)
  meaned_y = y - mean(y, na.rm=TRUE)
  return( cov_ff(meaned_x, meaned_y, meaned=TRUE) / (sd_ff(meaned_x, meaned=TRUE)*sd_ff(meaned_y, meaned=TRUE)) )  
}
