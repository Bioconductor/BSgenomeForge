\name{forgeBSgenomeDataPkgFRomUCSC}

\alias{forgeBSgenomeDataPkgFromUCSC}

\title{Create a BSgenome data package from a UCSC genome}

\description{
  The \code{forgeBSgenomeDataPkgFromUCSC} function allows the user to
  create a BSgenome data package from a UCSC genome.
}

\usage{
forgeBSgenomeDataPkgFromUCSC(genome, organism,
             pkg_maintainer, pkg_author=NA,
             pkg_version="1.0.0", pkg_license="Artistic-2.0",
             circ_seqs=NULL,
             goldenPath.url=getOption("UCSC.goldenPath.url"),
             destdir=".")
}

\arguments{
  \item{genome}{
    A single string specifying the name of a UCSC genome
    (e.g. \code{"mm39"} or \code{"sacCer3"}).
  }
  \item{organism}{
    The full name of the organism e.g. \code{"Mus musculus"},
    \code{"Saccharomyces cerevisiae"}, etc...
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
    NULL (the default), or a character vector containing the names of the
    circular sequences in the UCSC genome. This only needs to be specified if
    the genome is \emph{not} registered in the \pkg{GenomeInfoDb} package
    (if the genome is registered then its circular sequences are known so
    there's no need to specify \code{circ_seqs}).

    Notes:
    \itemize{
      \item You can use \code{registered_UCSC_genomes()} to get the list
            of UCSC genomes that are registered in the \pkg{GenomeInfoDb}
            package.
      \item If the genome is not registered and does not have
            circular sequences, then \code{circ_seqs} must be set
            to \code{character(0)}.
    }
  }
  \item{goldenPath.url}{
    A single string specifying the URL to the UCSC goldenPath location
    where the genomic sequences and chromosome sizes are expected to be found.
  }
  \item{destdir}{
    A single string containing the path to the directory where the
    BSgenome data package is to be created. This directory must already
    exist. Note that, by default, the package will be created in the
    current directory (\code{"."}).
  }
}

\details{
  This function is intended for Bioconductor users who want to forge a BSgenome
  data package from a UCSC genome. It typically makes use of the
  \code{downloadGenomicSequencesFromUCSC} utility function to download
  the 2bit file that contains the genomic sequences of the genome,
  and stores it in the working directory. However, if the file already
  exists in the working directory, then it is used and not downloaded again.
}

\value{
  The path to the created package as an invisible string.
}

\author{
  Hervé Pagès
}

\seealso{
  \itemize{
    \item The \code{\link[GenomeInfoDb]{registered_UCSC_genomes}} and
          \code{\link[GenomeInfoDb]{getChromInfoFromUCSC}} functions
          defined in the \pkg{GenomeInfoDb} package.

    \item The \code{\link{downloadGenomicSequencesFromUCSC}} function that
          \code{forgeBSgenomeDataPkgFromUCSC} uses internally to download
          the genomic sequences from UCSC.

    \item The \code{\link{forgeBSgenomeDataPkgFromNCBI}} function for
          creating a BSgenome data package from an NCBI assembly.
  }
}

\examples{
## Create a BSgenome data package for UCSC genome wuhCor1 (SARS-CoV-2
## assembly, see https://genome.ucsc.edu/cgi-bin/hgGateway?db=wuhCor1):
forgeBSgenomeDataPkgFromUCSC(
    genome="wuhCor1",
    organism="Severe acute respiratory syndrome coronavirus 2",
    pkg_maintainer="Jane Doe <janedoe@gmail.com>",
    destdir=tempdir()
)
}

\keyword{utilities}
