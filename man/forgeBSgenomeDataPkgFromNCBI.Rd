\name{forgeBSgenomeDataPkgFRomNCBI}

\alias{forgeBSgenomeDataPkgFromNCBI}

\title{Create a BSgenome data package from an NCBI assembly}

\description{
  The \code{forgeBSgenomeDataPkgFRomNCBI} function allows the user to
  create a BSgenome data package from an NCBI assembly.
}

\usage{
forgeBSgenomeDataPkgFromNCBI(assembly_accession, organism, genome,
                             pkg_maintainer, pkg_author=NA,
                             pkg_version="1.0.0",
                             pkg_license="Artistic-2.0",
                             circ_seqs=NULL,
                             destdir=".")
}

\arguments{
  \item{assembly_accession}{
    A single string containing a GenBank assembly accession (e.g.
    \code{"GCA_009729545.1"}) or a RefSeq assembly accession (e.g.
    \code{"GCF_009729545.1"}).
  }
  \item{organism}{
    The full name of the organism (e.g. \code{"Homo sapiens"},
    \code{"Felis catus"}, \code{"Acidianus infernus"}, etc...).
  }
  \item{genome}{
    The name of the NCBI assembly (e.g. \code{"GRCh38.p14"},
    \code{"felCat9.1_X"}, \code{"ASM972954v1"}, etc...).
  }
  \item{pkg_maintainer}{
    A single string containing the name and email address of the package
    maintainer (e.g \code{"Jane Doe, <janedoe@gmail.com>"}).
  }
  \item{pkg_author}{
    A single string containing the name of the package author. When
    unspecified, this takes the value of \code{pkg_maintainer}.
  }
  \item{pkg_version}{
    The version of the package. Set to \code{"1.0.0"} by default.
  }
  \item{pkg_license}{
    The license of the package. This must be the name of a software license
    used for free and open-source packages. Set to \code{"Artistic-2.0"}
    by default.
  }
    \item{circ_seqs}{
    The names of the circular sequences found in the genome. This must be
    specified if the organism is not registered in \code{GenomeInfoDb}, and
    supplied as a character vector. If the organism does not possess circular
    sequences, set this to character(0).
  }
  \item{destdir}{
    A single string containing the path to the directory where the
    compressed FASTA file is to be downloaded. This directory must
    already exist. Note that, by default, the file will be downloaded
    to the current directory (\code{"."}).
  }
}

\details{
  This function is intended for Bioconductor users who want
  to forge a BSgenome data package from an NCBI assembly.
}

\value{
  The path to the created package as an invisible string.
}

\author{
  Atuhurira Kirabo Kakopo
}

\seealso{
  \itemize{
    \item The \code{\link{downloadGenomicSequencesFromNCBI}} function that
          \code{forgeBSgenomeDataPkgFromNCBI} uses internally to download
          the genomic sequences from NCBI.

    \item The \code{\link{fastaTo2bit}} function that
          \code{forgeBSgenomeDataPkgFromNCBI} uses internally to convert
          the file downloaded by \code{downloadGenomicSequencesFromNCBI}
          from FASTA to 2bit.

  }
}

\examples{
## Create a BSgenome data package for NCBI assembly GCA_009729545.1
## (organism Acidianus infernus):
forgeBSgenomeDataPkgFromNCBI(assembly_accession="GCA_009729545.1",
                             organism="Acidianus infernus",
                             genome="ASM972954v1",
                             pkg_maintainer="Jane Doe <janedoe@gmail.com>",
                             circ_seqs=character(0))
}

\keyword{utilities}