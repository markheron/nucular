##' function variations to run on ff vector lists
##' 
##' Copies of several functions that are originally for numeric vectors, but in this *_ff_list variations work on lists of ff vectors
##' 
##' @name ff_list_functions
##' @author Mark Heron
##' 
##' @import ff
##' @import ffbase
NULL



##' hist_ff_list
##'
##' Plots a histogram for the values in column for every ff_list element.
##' @export
##' @param ff_list list with ff objects
##' @param column to plot the histogram for
##' @param ... parameters that can be passed on the hist
##' @author Mark Heron
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
##' Creates a list where each element of array_list is saved in as an ff object
##' @export
##' @param array_list to convert to a ff_list
##' @return list of ff objects
##' @author Mark Heron
array_list_to_ff_list <- function(array_list) {
  
  ff_list <- list()
  for(name in names(array_list)) {
    ff_list[[name]] <- as.ff(as.matrix(array_list[[name]]))
  }
  return(ff_list)
}


##' df_list_to_ffdf_list
##'
##' Creates a list where each element of df_list is saved in as an ffdf object
##' @export
##' @param df_list  to convert to a ffdf_list
##' @return list of ffdf objects
##' @author Mark Heron
df_list_to_ffdf_list <- function(df_list) {
  
  ff_list <- list()
  for(name in names(df_list)) {
    ff_list[[name]] <- as.ffdf(df_list[[name]])
  }
  return(ff_list)
}


##' filter_column_range_ff_list
##'
##' Filters each ff matrix in the list based on the column, min and max.
##' @export
##' @param ff_list to filter
##' @param column based on which to filter
##' @param min allowed column value
##' @param max allowed column value
##' @return list of ff's filtered by min/max on column
##' @author Mark Heron
filter_column_range_ff_list <- function(ff_list, column, min, max) {
  
  filtered_ff_list <- lapply(ff_list, function (x) as.ff(as.matrix(x[ (x[,column] >= min) & (x[,column] <= max) ,])))
  return(filtered_ff_list)
}



##' convertSparse2occ_ff_list
##'
##' Converts sparse representation of dyad positions to nucleosome occupancy as list of ff vectors.
##' @export
##' @param sparse list of sparse representation matricies
##' @param lengths list of chromosome lengths, names must match those of sparse
##' @return list of ff vectors with genomic occuopancy
##' @author Mark Heron
convertSparse2occ_ff_list <- function(sparse, lengths) {
  
  occ <- lapply(convertSparse2Complete_ff(sparse, lengths), function (x) smear_ff(x, -73,73))
  return(occ)
}



##' cov_ff_list
##' 
##' Computest the covariance between two ff vectors.
##' @export
##' @param x a list of ff vectors
##' @param y a second list of ff vectors
##' @param centered boolean, is the mean of x and y already = 0 ?
##' @return covariance between x and y
##' @author Mark Heron
cov_ff_list <- function(x,y,centered=FALSE) {
  
  if(centered) {
    return( sum( mapply(function (a,b) sum(a*b, na.rm=TRUE) ,x,y), na.rm=TRUE) )
  } else {
    return( sum( mapply(function (a,b) sum(a*b, na.rm=TRUE) , center_list(x), center_list(y)), na.rm=TRUE))
  }
}



##' sd_ff_list
##'
##' Computes the standard deviation of a vector.
##' @export
##' @param x a list of ff vectors
##' @param centered  boolean, is the mean already = 0 ?
##' @return standard deviation of x
##' @author Mark Heron
sd_ff_list <- function(x, centered=FALSE) {
  
  if(centered) {
    return( sqrt(sum( unlist(lapply( x, function (a) sum(a^2, na.rm=TRUE) )) , na.rm=TRUE)) )
  } else {
    return( sqrt(sum( unlist(lapply( center_list(x), function (a) sum(a^2, na.rm=TRUE) )) , na.rm=TRUE)) )
  }  
}



##' cor_ff_list
##' 
##' Computes the correlation between two vectors.
##' @export
##' @param x a list of ff vectors
##' @param y a second list of ff vectors
##' @return pearson correlation between x and y
##' @author Mark Heron
cor_ff_list <- function(x,y) {
  
  centered_x = center_list(x)
  centered_y = center_list(y)
  return( cov_ff_list(centered_x, centered_y, centered=TRUE) / (sd_ff_list(centered_x, centered=TRUE)*sd_ff_list(centered_y, centered=TRUE)) )  
}
