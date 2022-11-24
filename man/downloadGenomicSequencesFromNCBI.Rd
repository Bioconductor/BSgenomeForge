\name{downloadGenomicSequencesFromNCBI}

\alias{downloadGenomicSequencesFromNCBI}

\title{Download genomic sequences from NCBI}

\description{
  A utility function to download the compressed FASTA file that contains
  the genomic sequences of a given NCBI assembly.
}

\usage{
downloadGenomicSequencesFromNCBI(assembly_accession, assembly_name=NA,
                                 destdir=".", method, quiet=FALSE)
}

\arguments{
  \item{assembly_accession}{
    A single string containing a GenBank assembly accession (e.g.
    \code{"GCA_000001405.15"}) or a RefSeq assembly accession (e.g.
    \code{"GCF_000001405.26"}).
  }
  \item{assembly_name}{
    A single string or NA.
  }
  \item{destdir}{
    A single string containing the path to the directory where the
    compressed FASTA file is to be downloaded. This directory must already
    exist.
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
  to download the compressed FASTA file from NCBI for a given assembly
  specified by the \code{assembly_accession} argument.
}

\value{
  The path to the downloaded file as an invisible string.
}

\author{
  Prisca Chidimma Maduka
}

\seealso{
  \itemize{
    \item The \code{\link[utils]{download.file}} function in the \pkg{utils}
          package that \code{downloadGenomicSequencesFromNCBI} uses internally
          to download the compressed FASTA file.

    \item The \code{\link{downloadGenomicSequencesFromUCSC}} function to
          download genomic sequences from UCSC.
  }
}

\examples{
## Download the compressed FASTA file for NCBI assembly ASM972954v1 (see
## https://www.ncbi.nlm.nih.gov/assembly/GCF_009729545.1/):
downloadGenomicSequencesFromNCBI("GCF_009729545.1")

## Use the 'destdir' argument to specify the directory where to
## download the file:
downloadGenomicSequencesFromNCBI("GCF_009729545.1", destdir=tempdir())

## Download and import the file in R as a DNAStringSet object:
filepath <- downloadGenomicSequencesFromNCBI("GCF_009729545.1",
                                             destdir=tempdir())
genomic_sequences <- readDNAStringSet(filepath)
genomic_sequences
}

\keyword{utilities}
