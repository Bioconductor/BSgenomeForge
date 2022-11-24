\name{downloadGenomicSequencesFromUCSC}

\alias{downloadGenomicSequencesFromUCSC}

\title{Download genomic sequences from UCSC}

\description{
  A utility function to download the 2bit file that contains the genomic
  sequences of a given UCSC genome.
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
    A string set to \code{getOption("UCSC.goldenPath.url")}
    by default. \code{getOption("UCSC.goldenPath.url")} returns the
    goldenPath URL, \url{http://hgdownload.cse.ucsc.edu/goldenPath}.
  }
  \item{destdir}{
    A single string containing the path to the directory where the 2bit file
    is to be downloaded. This directory must already exist.
    Note that, by default, the file will be downloaded to the current
    directory (\code{"."}).
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
  The path to the downloaded file as an invisible string.
}

\author{
  Emmanuel Chigozie Elendu (Simplecodez)
}

\seealso{
  \itemize{
    \item The \code{\link[utils]{download.file}} function in the \pkg{utils}
          package that \code{downloadGenomicSequencesFromUCSC} uses internally
          to download the 2bit file.

    \item The \code{\link{downloadGenomicSequencesFromNCBI}} function to
          download genomic sequences from NCBI.
  }
}

\examples{
## Download the 2bit file for UCSC genome sacCer1:
downloadGenomicSequencesFromUCSC("sacCer1")

## Use the 'destdir' argument to specify the directory where to
## download the file:
downloadGenomicSequencesFromUCSC("sacCer1", destdir=tempdir())

## Download and import the file in R as a DNAStringSet object:
filepath <- downloadGenomicSequencesFromUCSC("sacCer1", destdir=tempdir())
genomic_sequences <- import(filepath)
genomic_sequences
}

\keyword{utilities}
