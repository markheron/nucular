##' Utility functions for analysing nucleosome data
##' 
##' a longer text explaining it
##' 
##' @name nucular-package
##' @author Mark Heron
##' @docType package
##' 
##' @import ff
##' @import ffbase
##' @import maRs
##' @import pRon
##' @import Biostrings
NULL




##' convertRaw2Rdata
##'
##' Converts raw data tables from the galaxy pipeline to lists of matricies that is saved as an Rdata for futher use.
##' @export
##' @param name of the data file
##' @param folder in which the file is to be found
##' @param extension of the data file, standard is ".tabular"
##' @author Mark Heron
convertRaw2Rdata <- function(name, folder="", extension=".tabular") {
  
  raw_table <- read.table.ffdf(file=paste0(folder,name,extension), sep="\t", comment.char="", header=FALSE, col.names=c("chr", "start", "stop", "length", NA), colClasses=c("factor", "integer", "integer", "integer","NULL"))
  
  chr_list_table <- list()
  
  for( chr_name in levels(raw_table[,1])) {
    chr_list_table[[chr_name]] <- subset(raw_table,chr == chr_name )[-1]
  }
  chr_table <- lapply(chr_list_table, as.ram)
  save(chr_table, file=paste0(name,"_chr_table.Rdata"), compress=TRUE)
  return(chr_table)
}




##' plotGenomicCutouts
##'
##' Create oligonucleotide frequency profile figure from genome positions
##' @export
##' @param pos ff_list of the nucleosome fragments
##' @param strand information at the moment ignored
##' @param size in either direction of the dyad position for which to count and plot the oligo frequencies
##' @param order of the oligonucleotides
##' @param genome_folder where to find the genome fasta files
##' @param chromosomes which chromosomes to use (must be contained in both genome_folder and pos)
##' @param sample if >0 only use the first sample positions per chromosome
##' @return olinucleotide frequency matrix, rows are the oligonucleotides, cols the positions around the dyad
##' @author Mark Heron
plotGenomicCutouts <- function(pos, strand, size, order, genome_folder, chromosomes, sample=0) {
  
  fasta_genome <- read_genome_fasta(genome_folder)
  
  freqs <- matrix(0, nrow=length(oligo_names(order+1)), size*2+1)
  
  for(chr in chromosomes) {
    
    freqs_chr <- matrix(0, nrow=length(oligo_names(order+1)), size*2+1)
    
    good <- (1:nrow(pos[[chr]]))[(pos[[chr]][,1] > size & pos[[chr]][,1] < width(fasta_genome[[paste0("chr",chr)]]) - size)]
    if(sample > 0) {
      good <- good[1:sample]
    }    
    
    chr_num <- fasta2num(fasta_genome[[paste0("chr",chr)]], order+1)
    chr_num[is.na(chr_num)] <- 0
    pos_chr <- pos[[chr]][good,]
    
    for(i in 1:nrow(freqs_chr)) {
      
      chr_bool <- chr_num == i
      for(j in 1:ncol(freqs_chr)) {
        freqs_chr[i,j] <- sum( chr_bool[pos_chr[,1]-size+j] * pos_chr[,2])
      }        
    }
    
    freqs <- freqs + freqs_chr/100000
  }
  
  freqs <- t(apply(freqs, 1, function (x) x/colSums(freqs)))
  rownames(freqs)  <- oligo_names(order+1) 
  
  
  plotOligoFreqs(freqs, x_pos=-size:size)
  
  invisible(freqs)
}




##' cor_nucs
##'
##' Calculates the pearson correlation between two nucleosome occupancy profiles given their sparse representation.
##' @export
##' @param data_list1 first ff_list of nucleosome positions
##' @param data_list2 second ff_list of nucleosome positions
##' @param lengths list of chromosome lengths (names must match subset of the data_list's )
##' @return correlation pearson correlation between the nucleosome occupancies
##' @author Mark Heron
cor_nucs <- function(data_list1, data_list2, lengths) {
  
  occ1 <- convertSparse2occ_ff_list(data_list1[names(lengths)], lengths)
  occ2 <- convertSparse2occ_ff_list(data_list2[names(lengths)], lengths)
  
  return(cor_ff_list(occ1, occ2))
}


##' cor_nucs_using_single_vecs
##'
##' Calculates the pearson correlation between two nucleosome occupancy profiles given their sparse representation.
##' Same as cor_nucs but uses an older non optimized method, being kept for testing purposis.
##' @param data_list1 first ff_list of nucleosome positions
##' @param data_list2 second ff_list of nucleosome positions
##' @param lengths list of chromosome lengths (names must match subset of the data_list's)
##' @return correlation pearson correlation between the nucleosome occupancies
##' @author Mark Heron
cor_nucs_using_single_vecs <- function(data_list1, data_list2, lengths) {
  
  complete1 <- convertSparse2Complete_ff(data_list1[names(lengths)], lengths)
  complete2 <- convertSparse2Complete_ff(data_list2[names(lengths)], lengths)
  
  occ1 <- smear_ff(complete1[[names(lengths)[1]]], -73,73)
  occ2 <- smear_ff(complete2[[names(lengths)[1]]], -73,73)
  for(chr in names(lengths)[-1]) {
    occ1 <- c(occ1, smear_ff(complete1[[chr]], -73,73))
    occ2 <- c(occ2, smear_ff(complete2[[chr]], -73,73))
  }
  return(cor_ff(occ1, occ2))
}



##' cor_nucs_multiple
##'
##' Computes all pairwise pearson correlations between the nucleosome occupancies.
##' @export
##' @param data_list list of ff_list of nucleosome dyad positions
##' @param lengths list of chromosome lengths (names must match subset of the data_list's)
##' @return top triangle matrix of the pearson correlations between nucleosome occupancies.
##' @author Mark Heron
cor_nucs_multiple <- function(data_list, lengths) {
  
  occ_list <- convertSparse2occ_ff_list_list(data_list, lengths)
  
  return(cor_ff_list_list(occ_list))
}





