##' function variations to run on lists of vectors, handling them as one vector
##' 
##' Copies of several functions that are originally for numeric vectors, but in this _list variations work on list of vectors
##' 
##' @name ff_list_list_functions
##' @author Mark Heron
##' 
##' @import ff
##' @import ffbase
NULL

##' cor_ff_list_list
##'
##' Computes the pairwise correlations between the list elements and saves them in a matrix
##' @export
##' @param data_list of the ff_list objects to compute the correlations between.
##' @return centered x (mean = 0)
##' @author Mark Heron
cor_ff_list_list <- function(data_list) {
  
  cor_matrix <- matrix(0, nrow=length(data_list), ncol=length(data_list))
  for(i in 1:(length(data_list)-1)) {
    for(j in (i+1):length(data_list)) {
      cor_matrix[i,j] <- cor_ff_list(data_list[[i]], data_list[[j]])
    }
  }
  cor_matrix <- cor_matrix + t(cor_matrix) + diag(1,length(data_list))
  
  colnames(cor_matrix) <- names(data_list)
  rownames(cor_matrix) <- names(data_list)
  return(cor_matrix)
}


##' convertSparse2occ_ff_list_list
##'
##' Individually converts a list of different sparse representation of dyad positions to nucleosome occupancy as list of ff vectors.
##' @export
##' @param data_list list of different sparse representation matrix lists
##' @param lengths list of chromosome lengths, names must match those of data_list subelements
##' @return list of lists of ff vectors with genomic occuopancy
##' @author Mark Heron
convertSparse2occ_ff_list_list <- function(data_list, lengths) {
  
  occ_list <- list()
  for(data_name in names(data_list)) {
    
    occ_list[[data_name]] <- convertSparse2occ_ff_list(data_list[[data_name]][names(lengths)], lengths)
  }
  return(occ_list)
}