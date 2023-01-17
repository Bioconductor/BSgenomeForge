### =========================================================================
### forgeBSgenomeDataPkgFRomNCBI()
### -------------------------------------------------------------------------
###
### Create a BSgenome data package from an NCBI assembly.
###


.format_organism <- function(organism)
{
    ## Remove leading and trailing whitespaces (like strip() in Python):
    organism <- gsub("^\\s*|\\s*$", "", organism)
    ## Replace multiple whitespaces with single space:
    organism <- gsub("\\s+", " ", organism)
    first_letter <- substring(organism, 1, 1)
    other_letters <- substring(organism, 2)
    paste0(toupper(first_letter), tolower(other_letters))
}

.abbreviate_organism_name <- function(organism)
{
    first_letter <- substring(organism, 1, 1)
    parts <- strsplit(organism, " +")[[1]]
    last_part <- tail(parts, 1)
    paste0(first_letter, last_part)
}

.create_pkgname <- function(abbr_organism, genome)
{
    assembly_name <- gsub("[^0-9a-zA-Z.]", "", genome)
    paste0("BSgenome.", abbr_organism, ".NCBI.", assembly_name)
}

.create_pkgtitle <- function(organism, genome)
{
    paste0("Full genomic sequences for ", organism,
           " (NCBI assembly ", genome, ")")
}

.create_pkgdesc <- function(organism, genome, assembly_accession)
{
    paste0("Full genomic sequences for ", organism,
           " as provided by NCBI (assembly ", genome,
           ", assembly accession ", assembly_accession, "). ",
           "The sequences are stored in DNAString objects.")
}

.check_pkg_maintainer <- function(pkg_maintainer)
{
    pattern <- "\\<[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,}\\>"
    if (grepl(pattern, pkg_maintainer, ignore.case=TRUE))
        return(pkg_maintainer)
    stop(wmsg("Please enter a valid email address"))
}

.create_organism_biocview <- function(organism)
{
    parts <- strsplit(organism, " +")[[1]]
    first_part <- head(parts, 1)
    last_part <- tail(parts, 1)
    paste0(first_part, "_", last_part)
}

.get_all_seqnames_in_one_string <- function(assembly_accession)
{
    seqinfo <- getChromInfoFromNCBI(assembly_accession)
    seqnames <- seqinfo$SequenceName
    paste0('c', '(', paste0('"', seqnames, '"', collapse=","), ')')
}

.move_seq_file <- function(twobit_file, pkg_dir)
{
    to <- file.path(pkg_dir, "inst", "extdata", basename(twobit_file))
    file.rename(twobit_file, to)
}

forgeBSgenomeDataPkgFromNCBI <- function(assembly_accession, organism, genome,
                                         pkg_maintainer, pkg_author=NA,
                                         pkg_version="1.0.0",
                                         pkg_license="Artistic-2.0",
                                         destdir=".")
{
    if (!isSingleString(organism) || organism == "")
        stop(wmsg("'organism' must be a single (non-empty) string"))
    if (!isSingleString(genome) || genome == "")
        stop(wmsg("'genome' must be a single (non-empty) string"))
    if (!isSingleString(pkg_maintainer) || pkg_maintainer == "")
        stop(wmsg("'pkg_maintainer' must be a single (non-empty) string"))
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

    ## Download file and convert from FASTA to 2bit.
    fasta_file <- downloadGenomicSequencesFromNCBI(assembly_accession)
    twobit_file <- file.path(tempdir(), "single_sequences.2bit")
    fastaTo2bit(fasta_file, twobit_file, assembly_accession)

    organism <- .format_organism(organism)
    abbr_organism <- .abbreviate_organism_name(organism)
    pkgname <- .create_pkgname(abbr_organism, genome)
    pkgtitle <- .create_pkgtitle(organism, genome)
    pkgdesc <- .create_pkgdesc(organism, genome, assembly_accession)
    pkg_maintainer <- .check_pkg_maintainer(pkg_maintainer)
    organism_biocview <- .create_organism_biocview(organism)
    seqnames <- .get_all_seqnames_in_one_string(assembly_accession)

    symValues <- list(BSGENOMEOBJNAME=abbr_organism,
                      PKGTITLE=pkgtitle,
                      PKGDESCRIPTION=pkgdesc,
                      PKGVERSION=pkg_version,
                      AUTHOR=pkg_author,
                      MAINTAINER=pkg_maintainer,
                      LICENSE=pkg_license,
                      ORGANISM=organism,
                      GENOME=genome,
                      ORGANISMBIOCVIEW=organism_biocview,
                      SEQNAMES=seqnames,
                      CIRCSEQS="character(0)")

    origdir <- system.file("pkgtemplates", "NCBI_BSgenome_datapkg",
                           package="BSgenomeForge")
    pkg_dir <- unlist(createPackage(pkgname, destdir, origdir, symValues,
                                    unlink=TRUE, quiet=FALSE),
                      use.names=FALSE)

    .move_seq_file(twobit_file, pkg_dir)
    invisible(pkg_dir)
}
