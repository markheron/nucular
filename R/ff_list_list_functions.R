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


##' compare_ff_list_list
##'
##' Applies a comparison function pairwise between all data_list elements and saves the results in a matrix
##' @export
##' @param data_list of the ff_list objects to compare
##' @param comparison_function comparison function to use, e.g. cor_ff_list
##' @return comparison matrix
##' @author Mark Heron
compare_ff_list_list <- function(data_list, comparison_function) {
  
  result_matrix <- matrix(0, nrow=length(data_list), ncol=length(data_list))
  for(i in seq_along(data_list)) {
    for(j in (i+1):length(data_list)) {
      
      common_elements <- intersect(names(data_list[[i]]), names(data_list[[j]]))
      result_matrix[i,j] <- comparison_function(data_list[[i]][common_elements], data_list[[j]][common_elements])
      
      if(length(common_elements) < max(length(data_list[[i]]),length(data_list[[j]])) ) {
        warning("Elements of two list don't match!")
      }
    }
  }
  result_matrix <- result_matrix + t(result_matrix) + diag(1,length(data_list))
  
  colnames(result_matrix) <- names(data_list)
  rownames(result_matrix) <- names(data_list)
  return(result_matrix)
}



##' cor_ff_list_list
##'
##' Computes the pairwise correlations between the list elements and saves them in a matrix
##' @export
##' @param data_list of the ff_list objects to compute the correlations between.
##' @return correlation matrix
##' @author Mark Heron
cor_ff_list_list <- function(data_list) {
  
  return( compare_ff_list_list(data_list, cor_ff_list) )
}



##' compare_ff_list_list_vs_ff_list_list
##'
##' Applies a comparison function pairwise between the list elements of two different lists and saves the results in a matrix
##' @export
##' @param first_list of the ff_list objects to compute the correlations from.
##' @param second_list of the ff_list objects to compute the correlations to.
##' @param comparison_function comparison function to use, e.g. cor_ff_list
##' @return comparison matrix
##' @author Mark Heron
compare_ff_list_list_vs_ff_list_list <- function(first_list, second_list, comparison_function) {
  
  result_matrix <- matrix(0, nrow=length(first_list), ncol=length(second_list))
  for(i in 1:(length(first_list))) {
    for(j in 1:length(second_list)) {
      common_elements <- intersect(names(first_list[[i]]), names(second_list[[j]]))
      result_matrix[i,j] <- comparison_function(first_list[[i]][common_elements], second_list[[j]][common_elements])
      
      if(length(common_elements) < max(length(first_list[[i]]),length(second_list[[j]])) ) {
        warning("Elements of two list don't match!")
      }
    }
  }
  
  colnames(result_matrix) <- names(second_list)
  rownames(result_matrix) <- names(first_list)
  return(result_matrix)
}



##' cor_ff_list_list_vs_ff_list_list
##'
##' Computes the pairwise correlations between the list elements of two different lists and saves them in a matrix
##' @export
##' @param first_list of the ff_list objects to compute the correlations from.
##' @param second_list of the ff_list objects to compute the correlations to.
##' @return correlation matrix
##' @author Mark Heron
cor_ff_list_list_vs_ff_list_list <- function(first_list, second_list) {
  
  return(compare_ff_list_list_vs_ff_list_list(first_list, second_list, cor_ff_list))
  
}


##' likelihood_ff_list_list_vs_ff_list_list
##'
##' Computes the pairwise likelihood between the list elements of two different lists and saves them in a matrix
##' @export
##' @param first_list of the ff_list objects
##' @param second_list of the ff_list objects
##' @return likelihood matrix
##' @author Mark Heron
likelihood_ff_list_list_vs_ff_list_list <- function(first_list, second_list) {
  
  return(compare_ff_list_list_vs_ff_list_list(first_list, second_list, likelihood_ff_list))
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
