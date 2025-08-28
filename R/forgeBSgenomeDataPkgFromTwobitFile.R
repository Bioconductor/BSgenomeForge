### =========================================================================
### forgeBSgenomeDataPkgFromTwobitFile()
### -------------------------------------------------------------------------
###
### Create a BSgenome data package from a 2bit file.
###


.get_seqnames_from_twobit_file <- function(filepath)
{
    if (!isSingleString(filepath) || filepath == "")
        stop(wmsg("'filepath' must be a single (non-empty) string"))
    if (!file.exists(filepath))
        stop(wmsg(filepath, ": no such file or directory"))
    twobitfile <- rtracklayer::TwoBitFile(filepath)
    ans <- try(seqlevels(twobitfile), silent=TRUE)
    if (inherits(ans, "try-error"))
        stop(wmsg("Failed to read file ", filepath, ". Is this a 2bit file?"))
    ans
}

.check_seqnames <- function(seqnames, filepath, file_seqnames)
{
    ## Same checks as check_circ_seqs() plus the "cannot be empty" check.
    if (is.null(seqnames))
        return()
    if (!is.character(seqnames))
        stop(wmsg("'seqnames' must be NULL or a character vector"))
    if (length(seqnames) == 0L)
        stop(wmsg("'seqnames' cannot be empty"))
    ## "primary key" constraints (see GenomeInfoDb:::is_primary_key)
    if (anyNA(seqnames))
        stop(wmsg("'seqnames' cannot contain NA's"))
    if (!all(nzchar(seqnames)))
        stop(wmsg("'seqnames' cannot contain empty strings"))
    if (anyDuplicated(seqnames))
        stop(wmsg("'seqnames' cannot contain duplicate values"))

    bad_seqnames <- setdiff(seqnames, file_seqnames)
    if (length(bad_seqnames) != 0L) {
        errmsg <- c("sequence name(s) in 'seqnames' not in file ",
                    filepath, ": ", paste(bad_seqnames, collapse=", "))
        stop(wmsg(errmsg))
    }
}

.normarg_circ_seqs <- function(circ_seqs, filepath, file_seqnames,
                               seqnames=NULL)
{
    check_circ_seqs(circ_seqs)  # shallow check

    if (is.null(seqnames)) {
        effective_seqnames <- file_seqnames
        where <- c("file ", filepath)
    } else {
        effective_seqnames <- seqnames
        where <- "'seqnames'"
    }
    if (is.null(circ_seqs)) {
        circ_flags <- tolower(effective_seqnames) %in%
                      tolower(DEFAULT_CIRC_SEQS)
        circ_seqs <- effective_seqnames[circ_flags]
        if (length(circ_seqs) == 0L) {
            msg <- c("based on their names we didn't see any circular ",
                     "sequences in ", where, " (this is just a guess!)")
        } else {
            msg <- c("the following circular sequences were guessed based ",
                     "on their names: ", paste(circ_seqs, collapse=", "))
        }
        warning(wmsg("you didn't specify 'circ_seqs' --> ", msg))
        return(circ_seqs)
    }

    ## Deep 'circ_seqs' check.
    bad_circ_seqs <- setdiff(circ_seqs, effective_seqnames)
    if (length(bad_circ_seqs) != 0L) {
        errmsg <- c("sequence name(s) in 'circ_seqs' not in ", where, ": ",
                    paste(bad_circ_seqs, collapse=", "))
        stop(wmsg(errmsg))
    }
    circ_seqs
}

.make_pkgtitle <- function(organism, provider, genome)
{
    paste0("Full genomic sequences for ", organism,
           " (genome assembly ", genome, " from ", provider, ")")
}

.make_pkgdesc <- function(organism, provider, genome)
{
    paste0("Full genomic sequences for ", organism, " as ",
           "provided by ", provider, " (genome assembly ", genome, "). ",
           "The sequences are stored in DNAString objects.")
}

.sort_twobit_file1 <- function(origfile, destfile, seqnames)
{
    dna <- rtracklayer::import.2bit(origfile)
    m <- match(seqnames, names(dna))
    ## 'seqnames' should have been checked already so this should never happen,
    ## but just in case...
    bad_idx <- which(is.na(m))
    if (length(bad_idx) != 0L) {
        errmsg <- c("sequence name(s) in 'seqnames' not in file ",
                    origfile, ": ", paste(seqnames[bad_idx], collapse=", "))
        stop(wmsg(errmsg))
    }
    dna <- dna[m]
    rtracklayer::export.2bit(dna, destfile)
    destfile
}

create_2bit_BSgenome_datapkg <-
    function(twobit_path, pkgname, BSgenome_objname,
             pkg_title, pkg_desc, pkg_version,
             pkg_author, pkg_maintainer, pkg_license,
             organism, provider, genome, organism_biocview, forge_function,
             circ_seqs, destdir=".", move_twobit_file=FALSE)
{
    stopifnot(isSingleString(twobit_path),
              isSingleString(pkgname),
              isSingleString(destdir),
              isTRUEorFALSE(move_twobit_file))
    origdir <- system.file("pkgtemplates", "2bit_BSgenome_datapkg",
                           package="BSgenomeForge", mustWork=TRUE)
    symValues <- list(BSGENOMEOBJNAME=BSgenome_objname,
                      PKGTITLE=pkg_title,
                      PKGDESCRIPTION=pkg_desc,
                      PKGVERSION=pkg_version,
                      PKGAUTHOR=pkg_author,
                      PKGMAINTAINER=pkg_maintainer,
                      PKGLICENSE=pkg_license,
                      ORGANISM=organism,
                      PROVIDER=provider,
                      GENOME=genome,
                      ORGANISMBIOCVIEW=organism_biocview,
                      CIRCSEQS=build_Rexpr_as_string(circ_seqs),
                      FORGEFUN=forge_function)
    stopifnot(all(vapply(symValues, isSingleString, logical(1))))
    pkg_dir <- createPackage(pkgname, destdir, origdir, symValues,
                             unlink=TRUE, quiet=FALSE)[[1L]]
    dir.create(file.path(pkg_dir, "inst", "extdata"), recursive = TRUE)
    to <- file.path(pkg_dir, "inst", "extdata", "single_sequences.2bit")
    if (move_twobit_file) {
        file.rename(twobit_path, to)
    } else {
        stopifnot(file.copy(twobit_path, to))
    }
    invisible(pkg_dir)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### forgeBSgenomeDataPkgFromTwobitFile()
###

forgeBSgenomeDataPkgFromTwobitFile <- function(filepath,
    organism, provider, genome,
    pkg_maintainer, pkg_author=NA,
    pkg_version="1.0.0",
    pkg_license="Artistic-2.0",
    seqnames=NULL,
    circ_seqs=NULL,  # must be a subset of 'seqnames' when 'seqnames' not NULL
    destdir=".")
{
    if (!isSingleString(organism) || organism == "")
        stop(wmsg("'organism' must be a single (non-empty) string"))
    if (!isSingleString(provider) || provider == "")
        stop(wmsg("'provider' must be a single (non-empty) string"))
    if (!isSingleString(genome) || genome == "")
        stop(wmsg("'genome' must be a single (non-empty) string"))
    check_pkg_maintainer(pkg_maintainer)
    if (identical(pkg_author, NA)) {
        pkg_author <- pkg_maintainer
    } else if (!isSingleString(pkg_author) || pkg_author == "") {
        stop(wmsg("'pkg_author' must be a single (non-empty) string"))
    }
    if (!isSingleString(pkg_version) || pkg_version == "")
        stop(wmsg("'pkg_version' must be a single (non-empty) string"))
    if (!isSingleString(pkg_license) || pkg_license == "")
        stop(wmsg("'pkg_license' must be a single (non-empty) string"))
    if (!isSingleString(destdir) || destdir == "")
        stop(wmsg("'destdir' must be a single (non-empty) string"))

    file_seqnames <- .get_seqnames_from_twobit_file(filepath)
    .check_seqnames(seqnames, filepath, file_seqnames)
    circ_seqs <- .normarg_circ_seqs(circ_seqs, filepath, file_seqnames,
                                    seqnames=seqnames)

    if (!is.null(seqnames)) {
        sorted_twobit_file <- file.path(tempdir(), "single_sequences.2bit")
        filepath <- .sort_twobit_file1(filepath, sorted_twobit_file, seqnames)
    }

    organism <- format_organism(organism)
    abbr_organism <- abbreviate_organism_name(organism)
    pkgname <- make_pkgname(abbr_organism, provider, genome)
    pkg_title <- .make_pkgtitle(organism, provider, genome)
    pkg_desc <- .make_pkgdesc(organism, provider, genome)
    biocview <- organism2biocview(organism)
    forge_function <- "forgeBSgenomeDataPkgFromTwobitFile"

    create_2bit_BSgenome_datapkg(filepath, pkgname,
                                 BSgenome_objname=abbr_organism,
                                 pkg_title=pkg_title,
                                 pkg_desc=pkg_desc,
                                 pkg_version=pkg_version,
                                 pkg_author=pkg_author,
                                 pkg_maintainer=pkg_maintainer,
                                 pkg_license=pkg_license,
                                 organism=organism,
                                 provider=provider,
                                 genome=genome,
                                 organism_biocview=biocview,
                                 forge_function=forge_function,
                                 circ_seqs=circ_seqs,
                                 destdir=destdir,
                                 move_twobit_file=!is.null(seqnames))
}

