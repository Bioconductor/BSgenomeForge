\name{downloadGenomicSequencesFromNCBI}

\alias{downloadGenomicSequencesFromNCBI}

\title{Download genomic sequences from NCBI}

\description{
  This is a utility function that finds and downloads
  the compressed FASTA file from NCBI, for a given assembly.
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
    A 'character' string indicating the path where
    the downloaded file is to be stored.
    It is set to \code{"."} by default.
  }
  \item{method, quiet}{
    Passed to the internal call to \code{download.file()}. See \code{?download.file}
    in the \pkg{utils} package for more information.
  }
}

\details{
  This function is intended for Bioconductor users who want
  to download the fasta file from NCBI for a given assembly
  specified by the \code{assembly_accession} argument.
}

\value{
  A 'character' string indicating the file path of the downloaded fasta file.
}

\author{
  Prisca Chidimma Maduka
}

\seealso{
  The \code{\link[utils]{download.file}} function in the \pkg{utils}
  package, that \code{downloadGenomicSequencesFromNCBI} uses internally
  to download the fasta file.
}

\examples{
## Download the fna file for NCBI genome ASM972954v1:
downloadGenomicSequencesFromNCBI("GCF_009729545.1")

## Use the 'destdir' argument to specify the directory where the
## file is to be downloaded:
downloadGenomicSequencesFromNCBI("GCF_009729545.1", destdir=tempdir())

## How to import the file in R as a DNAStringSet object:
file_path <- downloadGenomicSequencesFromNCBI("GCF_009729545.1", destdir=tempdir())
genome_sequence <- readDNAStringSet(file_path)
genome_sequence
}
