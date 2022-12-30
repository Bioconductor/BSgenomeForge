\name{forgeBSgenomeDataPkgFRomNCBI}

\alias{forgeBSgenomeDataPkgFromNCBI}

\title{Make a BSgenome package from a NCBI assembly}

\description{
  \code{forgeBSgenomeDataPkgFRomNCBI} is a function used to make a BSgenome
   package from an NCBI assembly.
}

\usage{
forgeBSgenomeDataPkgFromNCBI(assembly_accession, organism, genome,
                             pkg_maintainer, pkg_author=NA,
                             pkg_version="0.01",
                             license="Artistic-2.0",
                             destdir=".")
}

\arguments{
  \item{assembly_accession}{
    A single string containing a GenBank assembly accession (e.g.
    \code{"GCA_009729545.1"}) or a RefSeq assembly accession (e.g.
    \code{"GCF_009729545.1"}).
  }
  \item{organism}{
    The full name of the organism (e.g. "Homo sapiens", "Felis catus", 
    "Acidianus infernus", etc...).
  }
  \item{genome}{
    The name of the NCBI assembly (e.g. "GRCh38.p14", "felCat9.1_X", 
    "ASM972954v1", etc...).
  }
  \item{pkg_maintainer}{
    A single string containing the name and email address of the package
    maintainer (e.g Jane Doe, <janedoe@gmail.com>).
  }
  \item{pkg_author}{
    A single string containing the name of the package author. When unspecified
    , this takes the value of pkg_author.
  }
  \item{pkg_version}{
    The current version of the created BSgenome package.
  }
  \item{license}{
    A single string containing the software license used for free and 
    open-source packages. The default version is Artistic-2.0.
  }
  \item{destdir}{
    A single string containing the path to the directory where the 
    compressed FASTA file is to be downloaded. This directory must
    already exist. Note that, by default, the file will be downloaded
    to the current directory (".").
  }
}

\details{
  This function is intended for Bioconductor users who want
  to forge a BSgenome data package from an NCBI assembly.
}

\value{
  A logical TRUE if the package is created, a logical FALSE if it is not.
}

\author{
  Atuhurira Kirabo Kakopo
}

\seealso{
  \itemize{
    \item The \code{\link{downloadGenomicSequencesFromNCBI}} function that
          downloads genomic sequences from NCBI.
          
    \item The \code{\link{fastaTo2bit}} function that converts a FASTA file
          to the 2bit format.

  }
}

\examples{
## Make a BSgenome package for NCBI assembly GCA_009729545.1, organism 
## Acidianus infernus:
forgeBSgenomeDataPkgFromNCBI("GCA_009729545.1", "Acidianus infernus",
 "ASM972954v1", "Jane Doe <janedoe@gmail.com>")
}

\keyword{utilities}
