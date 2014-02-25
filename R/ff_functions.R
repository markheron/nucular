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
    return( sum( x*y ) )
  } else {
    return( sum( (x- mean(x, na.rm=TRUE))*(y-mean(y, na.rm=TRUE)) ) )
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
    return( sqrt(sum( x^2 ) ) )
  } else {
    return( sqrt(sum( (x-mean(x, na.rm=TRUE))^2 )) )
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
