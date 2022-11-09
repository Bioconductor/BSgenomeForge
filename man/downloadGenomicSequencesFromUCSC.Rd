\name{downloadGenomicSequencesFromUCSC}

\alias{downloadGenomicSequencesFromUCSC}

\title{The downloadGenomicSequencesFromUCSC function}

\description{
   A utility function that finds and downloads the 2bit file from
   UCSC, for a given genome.
 }

\usage{
downloadGenomicSequencesFromUCSC(genome, goldenPath.url=getOption("UCSC.goldenPath.url"), 
                                 destdir=".", method, quiet=FALSE)
 }

\arguments{
  \item{genome}{
    This is the name of the UCSC genome sequence to be downloaded.
    It is used to form the download URL.
  }
  
  \item{goldenPath.url}{
    a ‘character’ string set to getOption("UCSC.goldenPath.url)
    by default. getOption("UCSC.goldenPath.url) return the 
    goldenPath URL, http://hgdownload.cse.ucsc.edu/goldenPath. 
    This argument can be changed as needed. 
  }
  
  \item{destdir}{
    a 'character' string indicating the path where 
    the downloaded file is to be stored
    it is set to "." by default.   
  }
  
  \item{method, quiet}{
    Passed to the internal call to \code{download.file()}. See \code{?download.file}
    in the \pkg{utils} package for more information.
  }  
 }

\details{
  This function is intended for Bioconductor users who want 
  to download the 2bit genomic sequence file of a UCSC genome 
  specified by the 'genome' argument.
 }
  
  \seealso{
  The \code{\link[utils]{download.file}} function in the \pkg{utils} package, that
  \code{downloadGenomicSequencesFromUCSC} uses internally to download the 2bit file.
 }
  
\examples{
## Download a genome.2bit file 
## Downloading a genome.2bit file with only the genome passed
downloadGenomicSequencesFromUCSC("sacCer1")
   
## Downloading a genome.2bit file with genome and destdir arguments passed
downloadGenomicSequencesFromUCSC("sacCer1", destdir=tempdir())
   
## How to import the file in R as a DNAStringSet object. 
filepath <- downloadGenomicSequencesFromUCSC("sacCer1", destdir=tempdir())
genome_sequence <- import(filepath) 
genome_sequence
}
