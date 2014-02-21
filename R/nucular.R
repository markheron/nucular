##' Utility functions for analysing nucleosome data
##' 
##' a longer text explaining it
##' 
##' @name utility_functions
##' @author Mark Heron
##' 
##' @import ff
##' @import ffbase
##' @import pron
##' @import Biostrings
NULL


library(ff)
library(ffbase)
library(pron)
library(Biostrings)


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title table_to_ff
##' @param table_list 
##' @return list of ff objects
##' @author Mark Heron
table_to_ff <- function(table_list) {
  
  ff_list <- list()
  for(name in names(table_list)) {
    ff_list[[name]] <- as.ffdf(table_list[[name]])
  }
  return(ff_list)
}



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title filter_column_range_ff
##' @param table_list 
##' @param column 
##' @param min 
##' @param max 
##' @return list of ff's filtered by min/max on column
##' @author Mark Heron
filter_column_range_ff <- function(table_list, column, min, max) {
  
  table_list <- lapply(table_list, function (x) x[ (x[,column] >= min) & (x[,column] <= max) ,])
  return(table_list)
}



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title hist_ff_list
##' @param ff_list 
##' @param column 
##' @param ...
##' @author Mark Heron
hist_ff_list <- function(ff_list, column, ...) {
  hist( Reduce(c, lapply(ff_list, function (x) x[,column])), ...)
  gc()
}



##' cut's out single region from a single DNAString
##'
##' .. content for \details{} ..
##' @export
##' @title cut_out_fasta_single
##' @param fasta 
##' @param start 
##' @param end 
##' @param strand 
##' @return DNAString of cut out fragment
##' @author Mark Heron
cut_out_fasta_single <- function(fasta, start, end, strand) {
  
  fasta_seq <- subseq(fasta, start, end)
  if(strand == "-") {
    fasta_seq <- reverseComplement(fasta_seq)
  }
  return(fasta_seq)
}




##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title num2freq
##' @param seqnum 
##' @param mer_length 
##' @return matrix of oligonucleotide frequencies
##' @author Mark Heron
num2freq <- function(seqnum, mer_length) {
  
  mers <- create_mers(mer_length)
  freqs <- sapply(1:length(mers), function (mere) rowMeans(seqnum==mere, na.rm=TRUE))
  colnames(freqs) <- mers
  return(t(freqs))
}



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title num2weightedfreq
##' @param seqnum 
##' @param weights 
##' @param mer_length 
##' @return matrix of weighted oligonucleotide frequencies
##' @author Mark Heron
num2weightedfreq <- function(seqnum, weights, mer_length) {
  
  mers <- create_mers(mer_length)
  freqs <- sapply(1:length(mers), function (mere) colSums(apply((seqnum==mere), 1, function (x) x*weights), na.rm=TRUE))
  colnames(freqs) <- mers
  return(t(freqs))
}



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title cut_out_fasta_multiple_from_one_chr
##' @param pos 
##' @param strand 
##' @param size 
##' @param order 
##' @param chr_fasta 
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




##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title cut_out_seqnums_from_one_chr
##' @param pos 
##' @param strand 
##' @param size 
##' @param order 
##' @param chr_num 
##' @return matrix of cut out fragment
##' @author Mark Heron
cut_out_seqnums_from_one_chr <- function(pos, strand, size, order, chr_num) {
  
  for(i in 1:length(seqs)) {
    
    # seqs[i] <- cut_out_fasta_single(chr_fasta, start=pos[i]-size -((strand=="-")*order), end=pos[i]+size +((strand=="+")*order), strand )
    seqs <-  vapply(pos, function (i) chr_num[i+(-region -((strand=="-")*order)):(region +((strand=="+")*order))], FUN.VALUE=rep(0, length((-region -((strand=="-")*order)):(region +((strand=="+")*order)))))
    
  }
  return(seqs)
}



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title fasta2sparse
##' @param fasta 
##' @param mer_length 
##' @return matrix of sparse representation of the DNAString
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




##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title fasta2sparse_ff
##' @param fasta 
##' @param mer_length 
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



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title plotGenomicCutouts
##' @param pos 
##' @param strand 
##' @param size 
##' @param order 
##' @param genome_folder 
##' @param chromosomes 
##' @return frequency matrix
##' @author Mark Heron
plotGenomicCutouts <- function(pos, strand, size, order, genome_folder, chromosomes) {
  
  fasta_genome <- read_genome_fasta(genome_folder)
  
  freqs <- matrix(0, nrow=length(create_mers(order+1)), size*2+1)
  
  for(chr in chromosomes) { #[[1]]) {
    
    freqs_chr <- matrix(0, nrow=length(create_mers(order+1)), size*2+1)
    
    good <- (1:nrow(pos[[chr]]))[(pos[[chr]][,1] > size & pos[[chr]][,1] < width(fasta_genome[[paste0("chr",chr)]]) - size)] #[1:100000]
    
    #   chr_num <- fasta2num(genome[[paste0("chr",chr)]], mer_length=order+1 )
    #   seqnum <- t(cut_out_seqnums_from_one_chr(pos[[chr]][good,1], strand, size, order, chr_num))
    
    #     fasta <- cut_out_fasta_multiple_from_one_chr(pos[[chr]][good,1], "+", size, order, fasta_genome[[paste0("chr",chr)]])
    #   
    #     seqnum <- fasta2num(fasta, mer_length=order+1)
    #     
    #     freqs <- freqs + num2weightedfreq(seqnum, pos[[chr]][good,2], order+1)
    
    
    chr_num <- fasta2num(fasta_genome[[paste0("chr",chr)]], order+1)
    chr_num[is.na(chr_num)] <- 0
    pos_chr <- pos[[chr]][good,]
    
    for(i in 1:nrow(freqs_chr)) {
      
      chr_bool <- chr_num == i
      for(j in 1:ncol(freqs_chr)) {
        freqs_chr[i,j] <- sum( chr_bool[pos_chr[,1]-size+j] * pos_chr[,2])
      }        
    }
    
    
    #     chr_sparse <- fasta2sparse_ff(fasta_genome[[paste0("chr",chr)]], order+1)
    # 
    #     pos_chr <- pos[[chr]][good,]
    #     chr_length <- width(fasta_genome[[1]])-size
    #     frag_length <- 100000
    #     frag_start <- 1
    
    #     frag_nr <- 1+(pos_chr[,1] -1*size) %/% frag_length
    #     for(j in 1:(chr_length %/% frag_length)) {
    #       
    #      if(j == 2) {browser()}
    #       
    #       chr_sparse_frag <- chr_sparse[,frag_start+(0: min(frag_length+4*size, chr_length-frag_start-1))]
    #       for(p in pos_chr[frag_nr == j,]) {
    #         freqs_chr <- freqs_chr + chr_sparse_frag[,p[1]-frag_start +(-size):(size)]*p[2]
    #       }
    #       frag_start <- frag_start + frag_length
    #     }
    
    
    # chr_sparse_frag <- chr_sparse[,frag_start+(0: min(frag_length+4*size, chr_length-frag_start-1))]
    #     for(i in 1:nrow(pos_chr)) {
    #       pos_now <- pos_chr[i,1]
    #       if(pos_now >= (frag_start+frag_length+2*size)) {
    #         while(pos_now >= (frag_start+frag_length+2*size)) {
    #           frag_start <- frag_start+frag_length
    #         }
    #         chr_sparse_frag <- chr_sparse[,frag_start+(0: min(frag_length+4*size, chr_length-frag_start-1))]
    #       }
    #       freqs_chr <- freqs_chr + chr_sparse_frag[,pos_now-frag_start +(-size):(size)]*pos[[chr]][i,2]
    #     }
    
    #     i <- 1
    #     max_i <- nrow(pos_chr)
    #     for(frag_start in seq(1,chr_length, by=frag_length)) {
    #       pos_now <- pos_chr[i,1]
    #       chr_sparse_frag <- chr_sparse[,frag_start+(0: min(frag_length+4*size, chr_length-frag_start-1))]
    #       
    #       while(pos_now <= (frag_start+frag_length+2*size) & i < max_i) {
    #         freqs_chr <- freqs_chr + chr_sparse_frag[,pos_now-frag_start +(-size):(size)]*pos[[chr]][i,2]
    #         i <- i+1
    #         pos_now <- pos_chr[i,1]
    #       }
    #     }
    
    #     for(i in good) {
    #       freqs_chr <- freqs_chr + chr_sparse[,pos[[chr]][i,1]+(-size):(size)]*pos[[chr]][i,2]
    #     }
    
    freqs <- freqs + freqs_chr/100000
  }
  
  freqs <- t(apply(freqs, 1, function (x) x/colSums(freqs)))
  rownames(freqs)  <- create_mers(order+1) 
  
  
  plotMereFreqs(freqs, x_pos=-size:size)
  
  invisible(freqs)
}



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title convertSparse2Complete_ff
##' @param sparse 
##' @param lengths
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



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title smear_ff
##' @param x 
##' @param from 
##' @param to 
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



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title convertSparse2occ_ff
##' @param sparse 
##' @param lengths 
##' @return list of ff vectors with genomic occuopancy
##' @author Mark Heron
convertSparse2occ_ff <- function(sparse, lengths) {
  
  occ <- lapply(convertSparse2Complete_ff(sparse, lengths), function (x) smear_ff(x, -73,73))
  return(occ)
}



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title get_dyad_pos
##' @param data_list 
##' @param dyad_base 
##' @param offset 
##' @return list of matrices with dyad positions and intensity
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



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title cov_ff
##' @param x 
##' @param y 
##' @param meaned 
##' @return covariance
##' @author Mark Heron
cov_ff <- function(x,y,meaned=FALSE) {
  
  if(meaned) {
    return( sum( x*y ) )
  } else {
    return( sum( (x- mean(x, na.rm=TRUE))*(y-mean(y, na.rm=TRUE)) ) )
  }
}



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title sd_ff
##' @param x 
##' @param meaned 
##' @return standard deviation
##' @author Mark Heron
sd_ff <- function(x, meaned=FALSE) {
  
  if(meaned) {
    return( sqrt(sum( x^2 ) ) )
  } else {
    return( sqrt(sum( (x-mean(x, na.rm=TRUE))^2 )) )
  }  
}



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title cor_ff
##' @param x 
##' @param y 
##' @return pearson correlation
##' @author Mark Heron
cor_ff <- function(x,y) {
  
  meaned_x = x - mean(x, na.rm=TRUE)
  meaned_y = y - mean(y, na.rm=TRUE)
  return( cov_ff(meaned_x, meaned_y, meaned=TRUE) / (sd_ff(meaned_x, meaned=TRUE)*sd_ff(meaned_y, meaned=TRUE)) )  
}




##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title mean_list
##' @param x 
##' @return mean
##' @author Mark Heron
mean_list <- function(x) {
  
  x_sum <- sum(unlist(lapply(x, sum)))
  x_len <- sum(unlist(lapply(x, length)))
  
  return(x_sum/x_len)
}




##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title center_list
##' @param x 
##' @return list of x's after removing their overall mean
##' @author Mark Heron
center_list <- function(x) {
  return( lapply(x, function(a) a - mean_list(x)) )
}




##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title cov_ff_list
##' @param x 
##' @param y 
##' @param meaned 
##' @return covariance
##' @author Mark Heron
cov_ff_list <- function(x,y,meaned=FALSE) {
  
  if(meaned) {
    return( sum( mapply(function (a,b) sum(a*b) ,x,y)) )
  } else {
    return( sum( mapply(function (a,b) sum(a*b) , center_list(x), center_list(y))))
  }
}



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title sd_ff_list
##' @param x 
##' @param meaned 
##' @return standard deviation
##' @author Mark Heron
sd_ff_list <- function(x, meaned=FALSE) {
  
  if(meaned) {
    return( sqrt(sum( unlist(lapply( x, function (a) sum(a^2) )) )) )
  } else {
    return( sqrt(sum( unlist(lapply( center_list(x), function (a) sum(a^2) )) )) )
  }  
}



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title cor_ff_list
##' @param x 
##' @param y 
##' @return correlation
##' @author Mark Heron
cor_ff_list <- function(x,y) {
  
  meaned_x = center_list(x)
  meaned_y = center_list(y)
  return( cov_ff_list(meaned_x, meaned_y, meaned=TRUE) / (sd_ff_list(meaned_x, meaned=TRUE)*sd_ff_list(meaned_y, meaned=TRUE)) )  
}



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title cor_nucs
##' @param data_list1 
##' @param data_list2 
##' @param lengths 
##' @return correlation
##' @author Mark Heron
cor_nucs <- function(data_list1, data_list2, lengths) {
  
  occ1 <- convertSparse2occ_ff(data_list1[names(chr_lengths)], lengths)
  occ2 <- convertSparse2occ_ff(data_list2[names(chr_lengths)], lengths)
  
  return(cor_ff_list(occ1, occ2))
}



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title cor_nucs_using_single_vecs
##' @param data_list1 
##' @param data_list2 
##' @param lengths 
##' @return correlation
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



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @export
##' @title cor_nucs_multiple
##' @param data_list 
##' @param lengths 
##' @return correlation matrix
##' @author Mark Heron
cor_nucs_multiple <- function(data_list, lengths) {
  
  occ_list <- list()
  for(data_name in names(data_list)) {
    
    occ_list[[data_name]] <- ff(0, length=0)
    complete <- convertSparse2Complete_ff(data_list[[data_name]][names(lengths)], lengths)
    
    occ_list[[data_name]] <- smear_ff(complete[[chr]][[names(lengths)[1]]], -73,73)
    for(chr in names(lengths)[-1]) {
      occ_list[[data_name]] <- c(occ_list[[data_name]], smear_ff(complete[[chr]], -73,73))
    }
  }
  
  cor_matrix <- matrix(0, nrow=length(occ_list), ncol=length(occ_list))
  for(i in 1:(length(occ_list)-1)) {
    for(j in (i+1):length(occ_list)) {
      cor_matrix[i,j] <- cor_ff(occ_list[[i]], occ_list[[j]])
    }
  }
  return(cor_matrix)
}





