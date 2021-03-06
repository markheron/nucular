##' function variations to run on ff vector lists
##' 
##' Copies of several functions that are originally for numeric vectors, but in this *_ff_list variations work on lists of ff vectors
##' 
##' @name ff_list_functions
##' @author Mark Heron
##' 
##' @import ffbase
NULL



##' hist_ff_list
##'
##' Plots a histogram for the values in \code{column} for every ff_list element.
##' 
##' @export
##' @param ff_list list of ff matricies
##' @param column of the ff_list elements to plot the histogram for
##' @param ... parameters to be passed on to \code{\link{hist}}
##' 
hist_ff_list <- function(ff_list, column, ...) {
  hist( Reduce(c, lapply(ff_list, function (x) x[,column])), ...)
  gc()
}

# hist_ff_list <- function(ff_list, column, ...) {
#   hist( Reduce(c, lapply(ff_list, function (x) as.ff(x[,column]))), ...)
#   gc()
# }
# use whithin??


##' array_list_to_ff_list
##'
##' Creates a copy of array_list where each element is saved as an ff object.
##' 
##' @export
##' @param array_list to convert to a ff_list
##' @return list of ff objects
##' 
array_list_to_ff_list <- function(array_list) {
  
  ff_list <- list()
  for(name in names(array_list)) {
    if( !is.null(array_list[[name]]) & all(dim(array_list[[name]]) > 0) ) {
      ff_list[[name]] <- ff::as.ff(as.array(array_list[[name]]))
    } else {
      warning(paste("The array_list element", name, "is empty and skipped!"))
    }
  }
  return(ff_list)
}


##' df_list_to_ffdf_list
##'
##' Creates a copy of df_list where each element is saved as an ffdf object.
##' 
##' @export
##' @param df_list  to convert to a ffdf_list
##' @return list of ffdf objects
##' 
df_list_to_ffdf_list <- function(df_list) {
  
  ff_list <- list()
  for(name in names(df_list)) {
    ff_list[[name]] <- ff::as.ffdf(df_list[[name]])
  }
  return(ff_list)
}


##' filter_column_range_ff_list
##'
##' Filters the rows of the ff matrices in the list based on \code{column}, \code{min} and \code{max}.
##' 
##' @export
##' @param ff_list to filter
##' @param column based on which to filter
##' @param min allowed column value
##' @param max allowed column value
##' @return list of ff's filtered by min/max on column
##' 
filter_column_range_ff_list <- function(ff_list, column, min, max) {
  
  filtered_ff_list <- lapply(ff_list, function (x) {
    tmp <- (x[,column] >= min) & (x[,column] <= max)
    if(sum(tmp) > 0) {
      ff::as.ff(as.matrix(x[tmp ,]))
    } else {
      ff::as.ram(x)[0,] # slight hack since ff objects can't be empty
    }
  } )
  return(filtered_ff_list)
}


##' convertSparse2occ_ff_list
##'
##' Converts sparse representation of dyad positions to nucleosome occupancy ff vectors, for a list of data.
##' 
##' @export
##' @param sparse list of sparse representation matricies
##' @param lengths list of chromosome lengths, \code{names} must match those of \code{sparse}
##' @return list of ff vectors with genomic occuopancy
##' 
convertSparse2occ_ff_list <- function(sparse, lengths) {
  
  occ <- lapply(pRon::convertSparse2Complete_ff(sparse, lengths), function (x) maRs::smear_ff(x, -73,73))
  return(occ)
}



##' cov_ff_list
##' 
##' Computest the covariance between two ff vector lists.
##' 
##' @export
##' @param x list of ff vectors
##' @param y list of ff vectors
##' @param centered (boolean) is the mean_list of x and y == 0?
##' @return covariance between x and y
##' 
cov_ff_list <- function(x,y,centered=FALSE) {
  
  if(centered) {
    return( sum( mapply(function (a,b) sum(a*b, na.rm=TRUE) ,x,y), na.rm=TRUE) / (sum(mapply(function (a,b) sum(!is.na(a) & !is.na(b), na.rm=TRUE) ,x,y))-1) )
  } else {
    return( cov_ff_list( center_list(x), center_list(y), centered=TRUE) )
  }
}



##' sd_ff_list
##'
##' Computes the standard deviation of a ff vector list.
##' 
##' @export
##' @param x list of ff vectors
##' @param centered  boolean, is the mean_list of x == 0?
##' @return standard deviation of x
##' 
sd_ff_list <- function(x, centered=FALSE) {
  
  if(centered) {
    return( sqrt(sum( unlist(lapply( x, function (a) sum(a^2, na.rm=TRUE) )) , na.rm=TRUE) / (sum(unlist(lapply( x, function (a) sum(!is.na(a)))))-1) ) )
  } else {
    return( sd_ff_list( center_list(x), centered=TRUE) )
  }  
}



##' cor_ff_list
##' 
##' Computes the correlation between two lists of ff vectors.
##' 
##' @export
##' @param x list of ff vectors
##' @param y list of ff vectors
##' @return pearson correlation between x and y
##' 
cor_ff_list <- function(x,y) {
  
  not_na_pos <- mapply( function (a,b) !( is.na(a) | is.na(b) ), x, y, SIMPLIFY=FALSE)
  
  centered_x = center_list(mapply('[', x, not_na_pos, SIMPLIFY=FALSE))
  centered_y = center_list(mapply('[', y, not_na_pos, SIMPLIFY=FALSE))
  return( cov_ff_list(centered_x, centered_y, centered=TRUE) / (sd_ff_list(centered_x, centered=TRUE)*sd_ff_list(centered_y, centered=TRUE)) )  
}


##' center_ff_list
##' 
##' use center_list

##' scale_ff_list
##' 
##' use scale_list

#' intra_scaling_ff_list
#' 
#' use intre_scaling_list



##' make_prob_ff_list
##'
##' Scales all list elements so their sum is 1, i.e. if there are no negative values they can be used as probabilities.
##' 
##' @param x (list of ff_vectors)
##' @return (list of ff_vectors) probability version of x (sum == 1)
##' 
make_prob_ff_list <- function(x) {
  
  x_sum <- sum_list(x)
  return( lapply(x, function(a) a / x_sum ))
}



##' likelihood_ff_list
##' 
##' Computes a simple likelihood of a list of prediction vectors given a matching list of measurement vectors.
##' 
##' @export
##' @param predictions (list of ff vectors)
##' @param measurements (list of ff vectors)
##' @return (numeric) likelihood
##' 
likelihood_ff_list <- function(predictions, measurements) {
  
  # replace zero predictions with minimum predictions, so that the likelihood isn't -Inf because of them
  tmp_pred <- list()
  for(i in names(predictions)) {
    tmp_pred[[i]] <- ff::clone.ff(predictions[[i]])
    if( sum(tmp_pred[[i]] == 0, na.rm=TRUE) > 0) {
      tmp_pred[[i]][ (tmp_pred[[i]] == 0) & !(is.na(tmp_pred[[i]])) ] <- ff::as.ff(rep(min(tmp_pred[[i]][tmp_pred[[i]] > 0]), sum(tmp_pred[[i]] == 0, na.rm=TRUE)))
    }
  }
  
  return( sum( unlist(lapply( mapply( '*', lapply(make_prob_ff_list(tmp_pred), log), make_prob_ff_list(measurements), SIMPLIFY=FALSE) , sum, na.rm=TRUE)), na.rm=TRUE))
}



##' mae_ff_list
##' 
##' Computes the mean absolute error between two ff_lists.
##' 
##' @export
##' @param predictions (list of ff vectors)
##' @param measurements (list of ff vectors)
##' @return (numeric) mean absolute error
##' 
mae_ff_list <- function(predictions, measurements) {
  
  prob_pred <- make_prob_ff_list(predictions)
  prob_meas <- make_prob_ff_list(measurements)
  
  return( mean_list( lapply( mapply( '-', prob_pred, prob_meas, SIMPLIFY=FALSE), abs) ) )
}



##' rmse_ff_list
##' 
##' Computes the root mean square error between two ff_lists.
##' 
##' @export
##' @param predictions (list of ff vectors)
##' @param measurements (list of ff vectors)
##' @return (numeric) root mean square error
##' 
rmse_ff_list <- function(predictions, measurements) {
  
  prob_pred <- make_prob_ff_list(predictions)
  prob_meas <- make_prob_ff_list(measurements)
  
  return( sqrt( mean_list( mapply( function (a,b) {(a-b)^2} , prob_pred, prob_meas, SIMPLIFY=FALSE) ) ) )
}


