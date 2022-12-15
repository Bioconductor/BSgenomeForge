\name{fastaTo2bit}

\alias{fastaTo2bit}

\title{Convert files from FASTA to 2bit}

\description{
  \code{fastaTo2bit} is a utility function to convert a FASTA file
  to the 2bit format..
}

\usage{
fastaTo2bit(origfile, destfile, assembly_accession = NA)
}

\arguments{
  \item{origfile}{
    A single string containing the path to the FASTA file 
    (possibly compressed) to read, e.g. \code{"felCat9.fa"},
    \code{"felCat9.fa.gz"}, or \code{"path/to/felCat9.fa.gz"}.
  }
  \item{destfile}{
    A single string containing the path to the 2bit file to be
    written, e.g. \code{"felCat9.2bit"} or \code{"path/to/felCat9.2bit"}.
  }
  \item{assembly_accession}{
    A single string containing a GenBank assembly accession (e.g.
    \code{"GCA_009729545.1"}) or a RefSeq assembly accession (e.g.
    \code{"GCF_009729545.1"}).
  }
}

\details{
  This function is intended for Bioconductor users who want
  to convert a FASTA file to the 2bit format.
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

## Download the FASTA file containing the genomic sequences for
## this assembly to the tempdir() folder:
fasta_path <- downloadGenomicSequencesFromNCBI("GCF_009729545.1", 
                                                destdir=tempdir())

## Convert the FASTA file to the 2bit format:
twobitpath <- file.path(tempdir(), "ASM972954v1_sorted.2bit")
fastaTo2bit(fasta_path, twobitpath, assembly_accession = "GCF_009729545.1")
}

## Sanity check: try to import the 2bit file:
rtracklayer::import.2bit(twobit_path)
}

\keyword{utilities}
