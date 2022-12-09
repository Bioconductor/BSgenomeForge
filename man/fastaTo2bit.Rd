\name{fastaTo2bit}

\alias{fastaTo2bit}

\title{Convert files from FASTA format to 2bit}

\description{
  This is a utility function to convert a downloaded FASTA file to 2bit format
}

\usage{
fastaTo2bit <- function(origfile, destfile)
}

\arguments{
  \item{origfile}{
    A single string containing the path to the intended FASTA file 
    for conversion (e.g.felCat9.fa).
  }
  \item{destfile}{
    A converted 2bit file within the specified directory.
  }
}

\details{
  This function is intended for Bioconductor users who want
  to convert a downloaded file from FASTA to a 2bit format.
}

\value{
  A converted 2bit file within the specified directory.
}

\author{
  Atuhurira Kirabo Kakopo
}

