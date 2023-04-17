\name{BSgenomeForge-package}

\alias{BSgenomeForge-package}
\alias{BSgenomeForge}

\docType{package}

\title{The BSgenomeForge package}

\description{
  A package that simplifies the process of forging a BSgenome data package,
  by allowing the user to use one function to create the package.
  }

\details{
  \pkg{BSgenomeForge} provides two major functions, the
  \code{forgeBSgenomeDataPkgFromNCBI} function and
  \code{forgeBSgenomeDataPkgFromUCSC} function which allow one to forge
  a BSgenome data package from a NCBI assembly or UCSC genome respectively.

  For an overview of the functionality provided by the package, please see the
  vignette:
  \code{vignette("UsingBSgenomeForge", package="BSgenomeForge")}
  }

\author{
  Atuhurira Kirabo Kakopo,
  Hervé Pagès

  Maintainer: Hervé Pagès
  }

\seealso{
  \itemize{
    \item The \code{\link{forgeBSgenomeDataPkgFromNCBI}} function for
          creating a BSgenome data package from a NCBI assembly.

    \item The \code{\link{forgeBSgenomeDataPkgFromUCSC}} function for
          creating a BSgenome data package from a UCSC genome.
  }
}

\examples{
## ---------------------------------------------------------------------
## EXAMPLE 1
## ---------------------------------------------------------------------

## Create a BSgenome data package for NCBI assembly GCF_000857545.1
## (organism Torque teno virus 1):
forgeBSgenomeDataPkgFromNCBI(assembly_accession="GCF_000857545.1",
                             organism="Torque teno virus 1",
                             pkg_maintainer="Jane Doe <janedoe@gmail.com>",
                             circ_seqs="NC_002076.2",
                             destdir=tempdir())

## ---------------------------------------------------------------------
## EXAMPLE 2
## ---------------------------------------------------------------------

## Create a BSgenome data package for UCSC genome wuhCor1 (SARS-CoV-2
## assembly, see https://genome.ucsc.edu/cgi-bin/hgGateway?db=wuhCor1):
forgeBSgenomeDataPkgFromUCSC(
    genome="wuhCor1",
    organism="Severe acute respiratory syndrome coronavirus 2",
    pkg_maintainer="Jane Doe <janedoe@gmail.com>",
    destdir=tempdir()
)
}
