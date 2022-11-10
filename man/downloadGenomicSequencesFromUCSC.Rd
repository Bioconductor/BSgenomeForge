\name{downloadGenomicSequencesFromUCSC}

\alias{downloadGenomicSequencesFromUCSC}

\title{Download genomic sequences from UCSC}

\description{
   A utility function that finds and downloads the 2bit file from
   UCSC, for a given genome.
}

\usage{
downloadGenomicSequencesFromUCSC(
            genome,
            goldenPath.url=getOption("UCSC.goldenPath.url"),                   
            destdir=".", method, quiet=FALSE)
}

\arguments{
  \item{genome}{
    This is the name of the UCSC genome sequence to be downloaded.
    It is used to form the download URL.
  }
  
  \item{goldenPath.url}{
    A 'character' string set to \code{getOption("UCSC.goldenPath.url")}
    by default. \code{getOption("UCSC.goldenPath.url")} returns the 
    goldenPath URL, \url{http://hgdownload.cse.ucsc.edu/goldenPath}.  
  }
  
  \item{destdir}{
    A 'character' string indicating the path where 
    the downloaded file is to be stored.
    It is set to \code{"."} by default.   
  }
  
  \item{method, quiet}{
    Passed to the internal call to \code{download.file()}. See 
    \code{?download.file} in the \pkg{utils} package for more information.
  }  
}

\details{
  This function is intended for Bioconductor users who want 
  to download the 2bit genomic sequence file of a UCSC genome 
  specified by the \code{genome} argument.
}

\value{
  An (invisible) \code{destfile} argument
}
  
\author{
  Emmanuel Chigozie Elendu (Simplecodez)
}
  
\seealso{
  The \code{\link[utils]{download.file}} function in the \pkg{utils}
  package, that \code{downloadGenomicSequencesFromUCSC} uses internally
  to download the 2bit file.
}
  
\examples{
## Download the 2bit file for UCSC genome sacCer1:
downloadGenomicSequencesFromUCSC("sacCer1")
   
## Use the 'destdir' argument to specify the directory where the
## file is to be downloaded:
downloadGenomicSequencesFromUCSC("sacCer1", destdir=tempdir())
   
## How to import the file in R as a DNAStringSet object. 
filepath <- downloadGenomicSequencesFromUCSC("sacCer1", destdir=tempdir())
genome_sequence <- import(filepath) 
genome_sequence
}
