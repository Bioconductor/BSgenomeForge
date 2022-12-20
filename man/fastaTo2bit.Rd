\name{fastaTo2bit}

\alias{fastaTo2bit}

\title{Convert files from FASTA to 2bit}

\description{
  \code{fastaTo2bit} is a utility function to convert a FASTA file
  to the 2bit format.
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
    \code{"GCF_009729545.1"}). If missing, the function does not
    perform sequence renaming or reordering.
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
## Most assemblies at NCBI can be accessed using either their GenBank or
## RefSeq assembly accession. For example assembly ASM972954v1 (for
## Acidianus infernus) can be accessed either with GCA_009729545.1 (GenBank
## assembly accession) or GCF_009729545.1 (RefSeq assembly accession).
## See https://www.ncbi.nlm.nih.gov/assembly/GCA_009729545.1
## or https://www.ncbi.nlm.nih.gov/assembly/GCF_009729545.1 for
## the landing page of this assembly.

## ---------------------------------------------------------------------
## USING FASTA FILE FROM **GenBank** ASSEMBLY
## ---------------------------------------------------------------------

## Download the FASTA file containing the genomic sequences for
## the ASM972954v1 assembly to the tempdir() folder:
fasta_path <- downloadGenomicSequencesFromNCBI("GCF_009729545.1",
                                               destdir=tempdir())

## Use fastaTo2bit() to convert the file to 2bit. We're using the
## function in its simplest form here so there won't be any sequence
## renaming or reordering:
twobitpath1 <- tempfile(fileext=".2bit")
fastaTo2bit(fasta_path, twobitpath1)

## Take a look at the sequence names in the resulting 2bit file:
names(import.2bit(twobitpath1))

## Use fastaTo2bit() again to convert the file to 2bit. This time
## we want the function to rename and reorder the sequences as in
## getChromInfoFromNCBI("GCF_009729545.1"), so we supply "GCF_009729545.1":
twobitpath2 <- tempfile(fileext=".2bit")
fastaTo2bit(fasta_path, twobitpath2, assembly_accession="GCF_009729545.1")

## Take a look at the sequence names in the resulting 2bit file:
names(import.2bit(twobitpath2))

## ---------------------------------------------------------------------
## USING FASTA FILE FROM **RefSeq** ASSEMBLY
## ---------------------------------------------------------------------

#(same as above but using GCA_009729545.1 instead of GCF_009729545.1)
fasta_path <- downloadGenomicSequencesFromNCBI("GCA_009729545.1",
                                               destdir=tempdir())
twobitpath <- tempfile(fileext=".2bit")
fastaTo2bit(fasta_path, twobitpath, assembly_accession="GCA_009729545.1")
names(import.2bit(twobitpath))
}

\keyword{utilities}
