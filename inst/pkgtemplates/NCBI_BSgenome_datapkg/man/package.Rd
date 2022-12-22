\name{@PKGNAME@}
\docType{package}

\alias{@PKGNAME@-package}
\alias{@PKGNAME@}
\alias{@BSGENOMEOBJNAME@}

\title{@PKGTITLE@}

\description{
  @PKGDESCRIPTION@
}

\author{@AUTHOR@}

\seealso{
  \itemize{
    \item \link[BSgenome]{BSgenome} objects and the
          \code{\link[BSgenome]{available.genomes}} function
          in the \pkg{BSgenome} software package.
    \item \link[Biostrings]{DNAString} objects in the \pkg{Biostrings}
          package.
    \item The BSgenomeForge vignette (\code{vignette("BSgenomeForge")})
          in the \pkg{BSgenome} software package for how to make a BSgenome
          data package.
  }
}

\examples{
@PKGNAME@
genome <- @PKGNAME@
head(seqlengths(genome))

## ---------------------------------------------------------------------
## Genome-wide motif searching
## ---------------------------------------------------------------------
## See the GenomeSearching vignette in the BSgenome software
## package for some examples of genome-wide motif searching using
## Biostrings and the BSgenome data packages:
if (interactive())
    vignette("GenomeSearching", package="BSgenome")
}

\keyword{package}
\keyword{data}
