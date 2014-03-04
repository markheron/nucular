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
##' @import pron
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



##' num2freq
##'
##' Calculates the mere frequencies given a oligonucleotide-coded representation matrix
##' @export
##' @param seqnum oligonucleotide-coded representation of fasta sequences
##' @param mer_length length of the oligonucleotide-coding used (order+1)
##' @return matrix of oligonucleotide frequencies
##' @author Mark Heron
num2freq <- function(seqnum, mer_length) {
  
  mers <- create_mers(mer_length)
  freqs <- sapply(1:length(mers), function (mere) rowMeans(seqnum==mere, na.rm=TRUE))
  colnames(freqs) <- mers
  return(t(freqs))
}


##' num2weightedfreq
##'
##' Calculates the weighted mere frequencies given a oligonucleotide-coded representation matrix
##' @export
##' @param seqnum oligonucleotide-coded representation of fasta sequences
##' @param weights of the diferent fasta sequences
##' @param mer_length length of the oligonucleotide-coding used (order+1)
##' @return matrix of weighted oligonucleotide frequencies
##' @author Mark Heron
num2weightedfreq <- function(seqnum, weights, mer_length) {
  
  mers <- create_mers(mer_length)
  freqs <- sapply(1:length(mers), function (mere) colSums(apply((seqnum==mere), 1, function (x) x*weights), na.rm=TRUE))
  colnames(freqs) <- mers
  return(t(freqs))
}



##' cut_out_fasta_single
##' 
##' Cut's out single region from a single DNAString.
##' @export
##' @param fasta DNAString to cut out from
##' @param start of region to cut out
##' @param end of region to cut out
##' @param strand cut out forward or reverse complement
##' @return DNAString of cut out fragment
##' @author Mark Heron
cut_out_fasta_single <- function(fasta, start, end, strand) {
  
  fasta_seq <- subseq(fasta, start, end)
  if(strand == "-") {
    fasta_seq <- reverseComplement(fasta_seq)
  }
  return(fasta_seq)
}


##' cut_out_fasta_multiple_from_one_chr
##'
##' Cut's out multiple regions from a single DNAString (i.e. typically chromosome).
##' @export
##' @param pos vector of central position for cuting out the region
##' @param strand vector of "+"/"-" to specify if the forward or the reverse complement sequence shall be used
##' @param size region size before and after the positions to cut out
##' @param order what oligonucleotide order will later be used (to extend the cut region for symetry)
##' @param chr_fasta DNAString of the chromosome from which to cut out the regions
##' @return DNAStringSet of cut out fragment
##' @author Mark Heron
cut_out_fasta_multiple_from_one_chr <- function(pos, strand, size, order, chr_fasta) {
  
  seqs <- DNAStringSet(rep("",length(pos)))
  
  for(i in 1:length(seqs)) {
    
    # seqs[i] <- cut_out_fasta_single(chr_fasta, start=pos[i]-size -((strand=="-")*order), end=pos[i]+size +((strand=="+")*order), strand )
    seqs <- subseq(rep(chr_fasta, length(pos)), start=pos-size -((strand=="-")*order), end=pos+size +((strand=="+")*order))
    
  }
  names(seqs) <- ""
  return(seqs)
}




##' cut_out_seqnums_from_one_chr
##'
##' Cut's out multiple regions from a single oligonucleotide-coded respresentation vector (i.e. typically chromosome).
##' @export
##' @param pos vector of central position for cuting out the region
##' @param strand vector of "+"/"-" to specify if the forward or the reverse complement sequence shall be used
##' @param size region size before and after the positions to cut out
##' @param order what oligonucleotide order was used to code the fasta sequence
##' @param chr_num oligonucleotide-coded respresentation of the chromosome from which to cut out the regions
##' @return matrix of cut out fragment
##' @author Mark Heron
cut_out_seqnums_from_one_chr <- function(pos, strand, size, order, chr_num) {
  
  for(i in 1:length(seqs)) {
    
    # seqs[i] <- cut_out_fasta_single(chr_fasta, start=pos[i]-size -((strand=="-")*order), end=pos[i]+size +((strand=="+")*order), strand )
    seqs <-  vapply(pos, function (i) chr_num[i+(-size -((strand=="-")*order)):(size +((strand=="+")*order))], FUN.VALUE=rep(0, length((-size -((strand=="-")*order)):(size +((strand=="+")*order)))))
    
  }
  return(seqs)
}



##' fasta2sparse
##'
##' Creates sparse coding matrix of a fasta sequence.
##' @export
##' @param fasta DNAString
##' @param mer_length oligonucleotide length to code the sequence in
##' @return sparse representation matrix of the DNAString
##' @author Mark Heron
fasta2sparse <- function(fasta, mer_length) {
  
  fasta_num <- fasta2num(fasta, mer_length)
  
  fasta_sparse <- matrix(FALSE, nrow=length(create_mers(mer_length)) ,ncol=length(fasta_num))
  
  for(i in 1:nrow(fasta_sparse)) {
    fasta_sparse[i,] <- fasta_num==i
  }
  chr_sparse[is.na(chr_sparse)] <- FALSE
  
  return( fasta_sparse )
}




##' fasta2sparse_ff
##'
##' Creates sparse coding ff object of a fasta sequence.
##' @export
##' @param fasta DNAString
##' @param mer_length oligonucleotide length to code the sequence in
##' @return sparse ff.matrix representation of the DNAString
##' @author Mark Heron
fasta2sparse_ff <- function(fasta, mer_length) {
  
  fasta_num <- fasta2num(fasta, mer_length)
  
  fasta_sparse <- ff(FALSE, vmode="boolean", dim=c(length(create_mers(mer_length)), length(fasta_num)))
  
  for(i in 1:nrow(fasta_sparse)) {
    fasta_sparse[i,] <- fasta_num==i
  }
  
  return( fasta_sparse )
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
  
  freqs <- matrix(0, nrow=length(create_mers(order+1)), size*2+1)
  
  for(chr in chromosomes) {
    
    freqs_chr <- matrix(0, nrow=length(create_mers(order+1)), size*2+1)
    
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
  rownames(freqs)  <- create_mers(order+1) 
  
  
  plotMereFreqs(freqs, x_pos=-size:size)
  
  invisible(freqs)
}



##' convertSparse2Complete_ff
##'
##' Converts sparse representation to complete representation as ff vectors of genomic tags.
##' @export
##' @param sparse list of sparse representation matricies
##' @param lengths list of chromosome lengths, names must match those of sparse
##' @return list of ff vectors with complete representation of genomic tags
##' @author Mark Heron
convertSparse2Complete_ff <- function(sparse, lengths) {
  
  complete <- list()
  
  for(chr in names(sparse)) {
    complete[[chr]] <- ff(0, length=lengths[[chr]])
    complete[[chr]][sparse[[chr]][,1]] <- sparse[[chr]][,2]
  }
  invisible(complete)
}




##' get_dyad_pos
##'
##' Extracts the dyad positions in the chosen way from a table with mapped fragment start and ends
##' @export
##' @param data_list list of mapped fragments for each chromosome
##' @param dyad_base based on what position should the dyad position be calculated
##' @param offset offset if the dyad position isn't calculated from the center
##' @return list of ff matricies with dyad positions and intensity
##' @author Mark Heron
get_dyad_pos <- function(data_list, dyad_base="center", offset=73) {
  
  dyad_pos <- list()
  
  for(chr_name in names(data_list)) {
    
    if(dyad_base == "center") {
      tmp <- as.ffdf(as.data.frame(table( floor((data_list[[chr_name]][,1] + data_list[[chr_name]][,2])/2) )))
    } else if(dyad_base == "start") {
      tmp <- as.ffdf(as.data.frame(table(data_list[[chr_name]][,1]+offset)))
    } else if(dyad_base == "end") {
      tmp <- as.ffdf(as.data.frame(table(data_list[[chr_name]][,2]-offset)))
    } else if(dyad_base == "dinucleosome") {
      tmp <- as.ffdf(as.data.frame(table( c(data_list[[chr_name]][,1]+offset, data_list[[chr_name]][,2]-offset))))
    }
    
    dyad_pos[[chr_name]] <- as.ff(matrix(c(as.numeric(as.character(tmp[,1])), tmp[,2]), ncol=2))
  }
  return(dyad_pos)
}



##' adjust_X_chr
##'
##' Adjusts lower X chromosome counts due cells being male or a mixture of male/female.
##' @export
##' @param ff_list list of ff objects each representing data of one chromosome
##' @param X_chr name of the X chromosome in ff_list
##' @param Xfactor factor by which the counts should be adjusted (2 for male celllines, 4/3 for male/female mixtures)
##' @return adjusted ff_list
##' @author Mark Heron
adjust_X_chr <- function(ff_list, X_chr, Xfactor) {
  
  ff_list[[X_chr]][,2] <- ff_list[[X_chr]][,2]*Xfactor
  return(ff_list)
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





