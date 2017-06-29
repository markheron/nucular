##' nucleosomal nucleotide frequencies Automatic Report
##' 
##' Generates a set of figures for nucleosomal data, showing the fragment length distribution, occupancy, nucleotide frequencies around the mapped ends or dyad positions and correlations between different data sets
##' 
##' This package provides basic functions to handle nucleosome data, nuccontrol and nucppp actually do the reporting, so it might need a renaming.
##' 
##' @name nucular-package
##' @author Mark Heron
##' @docType package
##' 
##' @import ff
##' @import ffbase
##' @import Biostrings
##' @import maRs
##' @import pRon
NULL




##' convertRaw2Rdata
##'
##' Converts raw data tables from the galaxy pipeline to lists of matricies that is saved as an Rdata for futher use.
##' 
##' @export
##' @param name of the data file
##' @param folder in which the file is to be found,  default: ""
##' @param extension of the data file, default: ".tabular"
##' 
convertRaw2Rdata <- function(name, folder="", extension=".tabular") {
  
  warning("convertRaw2Rdata assumes the data has a 1-based genome index, not 0-based!")
  
  raw_table <- read.table.ffdf(file=paste0(folder,name,extension), sep="\t", comment.char="", header=FALSE, col.names=c("chr", "start", "stop", "length", NA), colClasses=c("factor", "integer", "integer", "integer","NULL"))
  
  chr_list_table <- list()
  
  for( chr_name in levels(raw_table[,1])) {
    chr_list_table[[chr_name]] <- subset(raw_table, chr == chr_name )[-1]
  }
  chr_table <- lapply(chr_list_table, as.ram)
  save(chr_table, file=paste0(name,"_chr_table.Rdata"), compress=TRUE)
  return(chr_table)
}




##' plotGenomicCutouts
##'
##' Create oligonucleotide frequency profile figure from genome positions
##' 
##' @export
##' @param pos ff_list of the nucleosome fragments
##' @param strand list of "+"/"-" vectors
##' @param size in either direction of the dyad position for which to count and plot the oligo frequencies
##' @param order of the oligonucleotides
##' @param genome character string of the folder where the genome fasta files are located, or a DNAStringSet of the genome
##' @param chromosomes which chromosomes to use (must be contained in both genome and pos)
##' @param sample if >0 only use the first sample positions per chromosome
##' @return olinucleotide frequency matrix, rows are the oligonucleotides, cols the positions around the dyad
##' 
plotGenomicCutouts <- function(pos, strand, size, order, genome, chromosomes, sample=0) {
    
  comp_oli_pos <- complementary_oligo_positions(order+1)
  
  if( is.character(genome)) {
    fasta_genome <- read_genome_fasta(genome)
    names(fasta_genome) <- sub("chr", "", names(fasta_genome))
  } else if (class(genome) == "DNAStringSet"){
    fasta_genome <- genome
  } else {
    stop("genome is neither a string nor a DNAStringSet !")
  }
  
  freqs <- matrix(0, nrow=length(oligo_names(order+1)), size*2+1)
  
  for(chr in chromosomes) {
    
    freqs_chr <- matrix(0, nrow=length(oligo_names(order+1)), size*2+1)
    
    good <- (1:nrow(pos[[chr]]))[( pos[[chr]][,1] > size+order & pos[[chr]][,1] < length(fasta_genome[[chr]]) - size )]
    if(sample > 0) {
      good <- good[1:sample]
    }    
    
    chr_num <- fasta2num(fasta_genome[chr], order+1, method="memoryLimited")
    chr_num[is.na(chr_num)] <- 0
    pos_chr <- pos[[chr]][good,]
    
    if(setequal(names(strand), names(pos))) {
      strand_chr <- strand[[chr]][good,]
    } else {
      strand_chr <- "+"
    }
    
    for(i in 1:nrow(freqs_chr)) {
      
      chr_bool <- chr_num == i
      for(j in 1:ncol(freqs_chr)) {
        # plus strand
        freqs_chr[i,j] <- freqs_chr[i,j] + sum( chr_bool[pos_chr[strand_chr=="+",1]-(size+1)+j] * pos_chr[strand_chr=="+",2])
        # minus strand
        if( any(strand_chr=="-") ) {
          freqs_chr[comp_oli_pos[i],j] <- freqs_chr[comp_oli_pos[i],j] + sum( chr_bool[pos_chr[strand_chr=="-",1]+(size+1)-j-order] * pos_chr[strand_chr=="-",2])
        }
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
##' 
##' @export
##' @param data_list1 first ff_list of nucleosome positions
##' @param data_list2 second ff_list of nucleosome positions
##' @param lengths list of chromosome lengths (names must match subset of the data_list's )
##' @return correlation pearson correlation between the nucleosome occupancies
##'
cor_nucs <- function(data_list1, data_list2, lengths) {
  
  occ1 <- convertSparse2occ_ff_list(data_list1[names(lengths)], lengths)
  occ2 <- convertSparse2occ_ff_list(data_list2[names(lengths)], lengths)
  
  return(cor_ff_list(occ1, occ2))
}


##' cor_nucs_using_single_vecs
##'
##' Calculates the pearson correlation between two nucleosome occupancy profiles given their sparse representation.
##' Same as cor_nucs but uses an older non-optimized method, being kept for testing purposis.
##' 
##' @param data_list1 first ff_list of nucleosome positions
##' @param data_list2 second ff_list of nucleosome positions
##' @param lengths list of chromosome lengths (\code{names} must match subset of the data_list's \code{names})
##' @return correlation pearson correlation between the nucleosome occupancies
##' 
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
##' 
##' @export
##' @param data_list list of ff_list of nucleosome dyad positions
##' @param lengths list of chromosome lengths (\code{names} must match subset of the data_list's \code{names})
##' @return top triangle matrix of the pearson correlations between nucleosome occupancies
##' 
cor_nucs_multiple <- function(data_list, lengths) {
  
  occ_list <- convertSparse2occ_ff_list_list(data_list, lengths)
  
  return(cor_ff_list_list(occ_list))
}




##' get_dyad_pos
##'
##' Extracts the dyad positions in the chosen way from a table with mapped fragment start and ends
##' 
##' @export
##' @param data_list list of mapped fragments for each chromosome
##' @param dyad_base one of "center" (default), "start". "end" and "dinucleosome"; Based on which feature should the dyad position be calculated ("dinucleosome" == "start"+"end")
##' @param offset (integer) offset if the dyad position isn't calculated from the center (\code{-offset} is used if \code{dyad_base=="end"})
##' @return list of ff matricies with dyad positions and intensities
##' 
get_dyad_pos <- function(data_list, dyad_base="center", offset=73) {
  
  dyad_pos <- list()
  
  to_ffdf_table <- function (x) {
    return( as.ffdf(as.data.frame(table(x))) )
  }
  
  for(chr_name in names(data_list)) {
    
    if(dyad_base == "center" | dyad_base == "center") {
      tmp <- to_ffdf_table( floor((data_list[[chr_name]][,1] + data_list[[chr_name]][,2])/2) )
    } else if(dyad_base == "start") {
      tmp <- to_ffdf_table( data_list[[chr_name]][,1]+offset )
    } else if(dyad_base == "end") {
      tmp <- to_ffdf_table( data_list[[chr_name]][,2]-offset )
    } else if(dyad_base == "dinucleosome") {
      tmp <- to_ffdf_table( c(data_list[[chr_name]][,1]+offset, data_list[[chr_name]][,2]-offset) )
    }
    
    dyad_pos[[chr_name]] <- as.ff(matrix(c(as.numeric(as.character(tmp[,1])), tmp[,2]), ncol=2))
  }
  return(dyad_pos)
}



##' adjust_X_chr
##'
##' Adjusts lower X chromosome counts due to cells being male or a mixture of male/female.
##' 
##' @export
##' @param ff_list list of ff objects each representing data of one chromosome. IMPORTANT: the data in this list will be changed!
##' @param X_chr (character) name of the X chromosome in \code{ff_list}
##' @param Xfactor (numeric) factor by which the counts should be adjusted (2 for male cell lines, 4/3 for male/female mixtures, i.e. embryo's or unspecific adult flies)
##' @return reference to the ff_list, which is now adjusted. IMPORTANT: this is the same ff_list as the input, no copy is created!
##'
adjust_X_chr <- function(ff_list, X_chr, Xfactor) {
  
  ff_list[[X_chr]][,2] <- ff_list[[X_chr]][,2]*Xfactor
  return(ff_list)
}



##' mask_repeats
##'
##' Masks repeat regions read from repeatMasker output files.
##' 
##' @export
##' @param occ_single (list of ff vectors) each ff representing one chromosome
##' @param mask_file_folder (character) name of the folder where the repeatMasker files for each chromosome can be found
##' @param repeats_to_mask (character vector) of types of repeats to mask
##' @param nuc_half_width (integer) amount to additionally mask to the left and right of the repeat
##' @return the masked occupancy ff vector list (cloned so the original is \emph{not} changed)
##' 
mask_repeats <- function(occ_single, mask_file_folder, repeats_to_mask, nuc_half_width=73) {
  
  occ_masked <- list()
  
  for(chr in names(occ_single)) {
    col_classes <- c(rep("NULL", 6), "integer", "integer", "NULL", rep("factor",3), rep("NULL", 5))
    
    repeats <- read.table(paste0(mask_file_folder,"/chr", chr, "_rmsk.txt"), header=FALSE, colClasses=col_classes)
    
    repeats <- repeats[ repeats[,5] %in% repeats_to_mask,] 
    
    repeats[,1] <- pmax(1, repeats[,1]-nuc_half_width)
    repeats[,2] <- pmin(length(occ_single[[chr]]), repeats[,2]+nuc_half_width)
    
    occ_masked[[chr]] <- clone(occ_single[[chr]], pattern="ff")
    for(i in 1:nrow(repeats)) {
      occ_masked[[chr]][repeats[i,1]:repeats[i,2]] <- NA
    }
    
  }
  return(occ_masked)
}



##' mask_Ns_in_ff_list
##'
##' Masks an ff_list based on N's in the genomic sequence.
##' 
##' @export
##' @param ff_list to be masked
##' @param genome_fasta (DNAStringSet) of the genome
##' @return masked ff_list (cloned so the original is \emph{not} changed)
##' 
mask_Ns_in_ff_list <- function(ff_list, genome_fasta) {
  
  masked_ff_list <- list()
  for(name in names(ff_list)) {
    masked_ff_list[[name]] <- clone(ff_list[[name]])
    
    fasta_num <- fasta2num(genome_fasta[name], 1)
    masked_ff_list[[name]][is.na(fasta_num)] <- NA
  }
  return(masked_ff_list)
}



##' uniquify_ff_list
##'
##' Removes duplicate fragments in a list of ff objects representing all reads for a chromosome.
##' 
##' @export
##' @param ff_list list of ff matrices each representing data of one chromosome
##' @return demultiplexed 
##' 
uniquify_ff_list <- function(ff_list) {
  if (!requireNamespace("plyr", quietly = TRUE)) {
    stop("plyr needed for this function to work. Please install it.", call. = FALSE)
  }
  uniquified_ff_list <- lapply(ff_list, function(x) as.ff(as.matrix( plyr::count(x)[, 1:3])))
  return(uniquified_ff_list)
}
