\name{@PKGNAME@}
\docType{package}

\alias{@PKGNAME@-package}
\alias{@PKGNAME@}

\title{@PKGTITLE@}

\description{
  @PKGDESCRIPTION@
}

\details{
  @PKGDETAILS@
}

\note{
  The masks in this BSgenome data package were made from the following
  source data files:
  \preformatted{
@SRCDATAFILES@
  }

  See \code{?\link[@REFPKGNAME@]{@REFPKGNAME@}} in the
  \pkg{@REFPKGNAME@} package for information about how the sequences
  were obtained.
}

\author{@PKGAUTHOR@}

\seealso{
  \itemize{
    \item \link[@REFPKGNAME@]{@REFPKGNAME@} in the \pkg{@REFPKGNAME@} package
          for information about how the sequences were obtained.

    \item \link[BSgenome]{BSgenome} objects in the \pkg{BSgenome}
          software package.

    \item The \code{\link[Seqinfo]{seqinfo}} getter and \link[Seqinfo]{Seqinfo}
          objects in the \pkg{Seqinfo} package.

    \item The \code{\link[GenomeInfoDb]{seqlevelsStyle}} getter and
          setter in the \pkg{GenomeInfoDb} package.

    \item \link[Biostrings]{MaskedDNAString} objects in the \pkg{Biostrings}
          package.

    \item The \code{\link[BSgenome]{available.genomes}} function
          in the \pkg{BSgenome} software package.

    \item The "Advanced BSgenomeForge usage" vignette in the
          \pkg{BSgenomeForge} package for how to create a BSgenome data
          package with masked sequences. Install the package and access the
          vignette with \code{browseVignettes(package="BSgenomeForge")}.
  }
}

\examples{
@PKGNAME@
mbsg <- @PKGNAME@
head(seqlengths(mbsg))
seqinfo(mbsg)

@PKGEXAMPLES@

if ("AGAPS" \%in\% masknames(mbsg)) {

  ## Check that the assembly gaps contain only Ns:
  checkOnlyNsInGaps <- function(seq)
  {
    ## Replace all masks by the inverted AGAPS mask
    masks(seq) <- gaps(masks(seq)["AGAPS"])
    unique_letters <- uniqueLetters(seq)
    if (any(unique_letters != "N"))
        stop("assembly gaps contain more than just Ns")
  }

  ## A message will be printed each time a sequence is removed
  ## from the cache:
  options(verbose=TRUE)

  for (seqname in seqnames(mbsg)) {
    cat("Checking sequence", seqname, "... ")
    seq <- mbsg[[seqname]]
    checkOnlyNsInGaps(seq)
    cat("OK\n")
  }
}

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
