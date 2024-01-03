### =========================================================================
### forgeBSgenomeDataPkgFromNCBI()
### -------------------------------------------------------------------------
###
### Create a BSgenome data package from an NCBI assembly.
###


.extract_assembly_name_from_ftp_file_prefix <- function(ftp_file_prefix)
{
    ## 'ftp_file_prefix' is a string that is expected to contain at least
    ## two underscores. We want to extract everything that comes after the
    ## second underscore.
    sub("^[^_]*_[^_]*_", "", ftp_file_prefix)
}

.fetch_assembly_name_from_NCBI <- function(assembly_accession)
{
    ## 'ftp_dir' will be a length 2 character vector containing URL
    ## to FTP dir and prefix of file names located in FTP dir.
    ftp_dir <- find_NCBI_assembly_ftp_dir(assembly_accession)
    .extract_assembly_name_from_ftp_file_prefix(ftp_dir[2])
}

.make_pkgname_for_NCBI_datapkg <- function(abbr_organism, assembly_name)
{
    part4 <- gsub("[^0-9a-zA-Z.]", "", assembly_name)
    paste0("BSgenome.", abbr_organism, ".NCBI.", part4)
}

.make_pkgtitle_for_NCBI_datapkg <- function(organism, assembly_name)
{
    paste0("Full genomic sequences for ", organism,
           " (NCBI assembly ", assembly_name, ")")
}

.make_pkgdesc_for_NCBI_datapkg <- function(organism, assembly_name,
                                           assembly_accession)
{
    paste0("Full genomic sequences for ", organism, " as ",
           "provided by NCBI (assembly ", assembly_name, ", assembly ",
           "accession ", assembly_accession, "). ",
           "The sequences are stored in DNAString objects.")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extract_NCBI_assembly_info()
###

### Returns a named list of length 4 with components "assembly_accession",
### "organism", "circ_seqs", and "assembly_name". This is a subset of
### the "NCBI_assembly_info" attribute found on the 'chrominfo' data frame
### returned by getChromInfoFromNCBI() when the assembly is registered.
.extract_NCBI_assembly_info <- function(assembly_accession, chrominfo,
                                        organism=NULL, circ_seqs=NULL)
{
    seqnames <- chrominfo[ , "SequenceName"]

    ## If the requested assembly is registered, then 'chrominfo' has
    ## a "NCBI_assembly_info" attribute.
    NCBI_assembly_info <- attr(chrominfo, "NCBI_assembly_info")
    if (is.null(NCBI_assembly_info)) {
        ## NCBI assembly is **not** registered.
        if (is.null(organism))
            stop(wmsg("\"", assembly_accession, "\" is not a registered NCBI ",
                      "assembly --> argument 'organism' must be supplied"))
        organism <- format_organism(organism)
        is_assembled <- chrominfo[ , "SequenceRole"] %in% "assembled-molecule"
        circ_seqs <- get_circ_seqs_for_unregistered_assembly_or_genome(
                         assembly_accession,
                         seqnames,
                         is_assembled,
                         circ_seqs,
                         what="assembly")
        ## Obtain assembly name from NCBI FTP repository.
        assembly <- .fetch_assembly_name_from_NCBI(assembly_accession)
    } else {
        ## NCBI assembly is registered.
        if (!is.null(organism))
            warning(wmsg("\"", assembly_accession, "\" is a registered NCBI ",
                         "assembly for organism \"",
                         NCBI_assembly_info$organism, "\" ",
                         "--> ignoring supplied 'organism' argument"),
                    immediate.=TRUE)
        organism <- NCBI_assembly_info$organism
        is_circular <- chrominfo[ , "circular"]
        circ_seqs <- get_circ_seqs_for_registered_assembly_or_genome(
                         assembly_accession,
                         seqnames,
                         is_circular,
                         circ_seqs,
                         what="assembly")
        assembly_accession <- NCBI_assembly_info$assembly_accession
        assembly <- NCBI_assembly_info$assembly
    }
    list(assembly_accession=assembly_accession,
         organism=organism,
         circ_seqs=circ_seqs,
         assembly=assembly)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### forgeBSgenomeDataPkgFromNCBI()
###

forgeBSgenomeDataPkgFromNCBI <- function(assembly_accession,
                                         pkg_maintainer, pkg_author=NA,
                                         pkg_version="1.0.0",
                                         pkg_license="Artistic-2.0",
                                         organism=NULL,
                                         circ_seqs=NULL,
                                         destdir=".")
{
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
    if (!(is.null(organism) || isSingleString(organism) && organism != ""))
        stop(wmsg("when suplied, 'organism' must be a single (non-empty) ",
                  "string"))
    check_circ_seqs(circ_seqs)
    if (!isSingleString(destdir) || destdir == "")
        stop(wmsg("'destdir' must be a single (non-empty) string"))

    ## Returns TRUE if 'assembly_accession' is a GenBank accession, or FALSE
    ## if it's a RefSeq accession, or an error if it's none.
    ## Make sure to do this before the call to getChromInfoFromNCBI() below.
    is_GCA <- is_GenBank_accession(assembly_accession)
    accession_type <- if (is_GCA) "GenBank" else "RefSeq"
    accession_col <- paste0(accession_type, "Accn")

    ## Retrieve chromosome information for specified NCBI assembly.
    chrominfo <- getChromInfoFromNCBI(assembly_accession)
    chrominfo <- drop_rows_with_NA_accns(chrominfo, accession_col)

    NCBI_assembly_info <- .extract_NCBI_assembly_info(assembly_accession,
                                                      chrominfo,
                                                      organism=organism,
                                                      circ_seqs=circ_seqs)
    assembly_accession <- NCBI_assembly_info$assembly_accession
    organism <- NCBI_assembly_info$organism
    circ_seqs <- NCBI_assembly_info$circ_seqs
    assembly_name <- NCBI_assembly_info$assembly
    seqnames <- chrominfo[ , "SequenceName"]

    ## Download genomic sequences and convert from FASTA to 2bit.
    file_url <- get_URL_to_genomic_sequences_from_NCBI(assembly_accession)
    fasta_file <- basename(file_url)
    if (file.exists(fasta_file)) {
        message(wmsg("File ", fasta_file, " is already in current ",
                     "directory so will be used."))
    } else {
        fasta_file <- downloadGenomicSequencesFromNCBI(assembly_accession)
    }
    sorted_twobit_file <- file.path(tempdir(), "single_sequences.2bit")
    fastaTo2bit(fasta_file, sorted_twobit_file, assembly_accession)

    abbr_organism <- abbreviate_organism_name(organism)
    pkgname <- .make_pkgname_for_NCBI_datapkg(abbr_organism, assembly_name)
    pkgtitle <- .make_pkgtitle_for_NCBI_datapkg(organism, assembly_name)
    pkgdesc <- .make_pkgdesc_for_NCBI_datapkg(organism, assembly_name,
                                              assembly_accession)
    check_pkg_maintainer(pkg_maintainer)
    biocview <- organism2biocview(organism)
    seqnames <- build_Rexpr_as_string(seqnames)
    circ_seqs <- build_Rexpr_as_string(circ_seqs)

    ## Create the package.
    origdir <- system.file("pkgtemplates", "NCBI_BSgenome_datapkg",
                           package="BSgenomeForge")
    symValues <- list(BSGENOMEOBJNAME=abbr_organism,
                      PKGTITLE=pkgtitle,
                      PKGDESCRIPTION=pkgdesc,
                      PKGVERSION=pkg_version,
                      PKGAUTHOR=pkg_author,
                      PKGMAINTAINER=pkg_maintainer,
                      PKGLICENSE=pkg_license,
                      ORGANISM=organism,
                      GENOME=assembly_name,
                      ORGANISMBIOCVIEW=biocview,
                      SEQNAMES=seqnames,
                      CIRCSEQS=circ_seqs)
    pkg_dir <- createPackage(pkgname, destdir, origdir, symValues,
                             unlink=TRUE, quiet=FALSE)[[1]]
    move_file_to_datapkg(sorted_twobit_file, pkg_dir)

    invisible(pkg_dir)
}

