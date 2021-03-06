% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nucular.R
\name{plotGenomicCutouts}
\alias{plotGenomicCutouts}
\title{plotGenomicCutouts}
\usage{
plotGenomicCutouts(pos, strand, size, order, genome, chromosomes,
  sample = 0)
}
\arguments{
\item{pos}{ff_list of the nucleosome fragments}

\item{strand}{list of "+"/"-" vectors}

\item{size}{in either direction of the dyad position for which to count and plot the oligo frequencies}

\item{order}{of the oligonucleotides}

\item{genome}{character string of the folder where the genome fasta files are located, or a DNAStringSet of the genome}

\item{chromosomes}{which chromosomes to use (must be contained in both genome and pos)}

\item{sample}{if >0 only use the first sample positions per chromosome}
}
\value{
olinucleotide frequency matrix, rows are the oligonucleotides, cols the positions around the dyad
}
\description{
Create oligonucleotide frequency profile figure from genome positions
}
