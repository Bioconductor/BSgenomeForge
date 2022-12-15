\name{fastaTo2bit}

\alias{fastaTo2bit}

\title{Convert files from FASTA format to 2bit}

\description{
  This is a utility function to convert a FASTA file to 2bit format.
}

\usage{
fastaTo2bit(origfile, destfile, assembly_accession = NA)
}

\arguments{
  \item{origfile}{
    A single string containing the path to the FASTA file 
    to be read (e.g. \code{"felCat9.fasta"})
  }
  \item{destfile}{
    A single string containing the path to the 2bit file to be written (e.g. \code{"felCat9.2bit"}).
  }
  \item{assembly_accession}{
    A single string containing a GenBank assembly accession (e.g.
    \code{"GCA_009729545.1"}) or a RefSeq assembly accession (e.g.
    \code{"GCF_009729545.1"}).
  }
}

\details{
  This function is intended for Bioconductor users who want
  to convert a downloaded file from FASTA to a 2bit format.
}

\value{
  An invisible NULL.
}

\author{
  Atuhurira Kirabo Kakopo
}

\seealso{
  \itemize{
    \item The \code{\link{downloadGenomicSequencesFromNCBI}} function to
          download genomic sequences from NCBI.
          
    \item The \code{\link{readDNAStringSet}} function in the \pkg{Biostrings} 
          package that \code{fastaTo2bit} uses internally to import 
          the FASTA file.
          
    \item The \code{\link{export.2bit}} function in the \pkg{rtracklayer} 
          package that \code{fastaTo2bit} uses internally to export 
          the 2bit file.
  }
}

\examples{
## See https://www.ncbi.nlm.nih.gov/assembly/GCF_009729545.1 for
## the landing page of the ASM972954v1 assembly (Acidianus infernus).

## Download the compressed FASTA file for NCBI assembly ASM972954v1 
## and store in the tempdir folder (see 
## https://www.ncbi.nlm.nih.gov/assembly/GCF_009729545.1/):
fasta_path <- downloadGenomicSequencesFromNCBI("GCF_009729545.1", 
                                                destdir=tempdir())

## Convert the FASTA file to 2bit format and specify assembly accession:
twobitpath <- file.path(tempdir(), "ASM972954v1_sorted.2bit")
fastaTo2bit(fasta_path, twobitpath, assembly_accession = "GCF_009729545.1")
}

## Sanity check: try to import the 2bit file:
rtracklayer::import.2bit(twobit_path)
}

\keyword{utilities}
