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
### .get_circ_seqs_from_NCBI()
###

### 'circ_seqs' contains the circular sequences specified by the user.
.get_circ_seqs_from_NCBI <- function(assembly_accession, chrominfo,
                                     circ_seqs=NULL)
{
    check_circ_seqs(circ_seqs)
    seqnames <- chrominfo[ , "SequenceName"]
    NCBI_assemblies <- registered_NCBI_assemblies()[ , "assembly_accession"]
    is_registered <- assembly_accession %in% NCBI_assemblies
    if (is_registered) {
        ## NCBI assembly is registered.
        FUN <- get_circ_seqs_for_registered_assembly_or_genome
        is_xxx <- chrominfo[ , "circular"]
    } else {
        ## NCBI assembly is **not** registered.
        FUN <- get_circ_seqs_for_unregistered_assembly_or_genome 
        is_xxx <- chrominfo[ , "SequenceRole"] %in% "assembled-molecule"
    }
    FUN(assembly_accession, seqnames, is_xxx, circ_seqs, what="assembly")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### forgeBSgenomeDataPkgFromNCBI()
###

forgeBSgenomeDataPkgFromNCBI <- function(assembly_accession, organism,
                                         pkg_maintainer, pkg_author=NA,
                                         pkg_version="1.0.0",
                                         pkg_license="Artistic-2.0",
                                         circ_seqs=NULL,
                                         destdir=".")
{
    if (!isSingleString(organism) || organism == "")
        stop(wmsg("'organism' must be a single (non-empty) string"))
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

    ## Retrieve chromosome information for specified NCBI assembly.
    chrominfo <- getChromInfoFromNCBI(assembly_accession)
    circ_seqs <- .get_circ_seqs_from_NCBI(assembly_accession, chrominfo,
                                          circ_seqs)

    ## Obtain assembly name from NCBI FTP repository.
    assembly_name <- .fetch_assembly_name_from_NCBI(assembly_accession)

    ## Download genomic sequences and convert from FASTA to 2bit.
    file_url <- get_URL_to_genomic_sequences_from_NCBI(assembly_accession)
    fasta_file <- basename(file_url)
    if (file.exists(fasta_file)) {
        message("File ", fasta_file, " is already in current ",
                "directory so will be used.")
    } else {
        fasta_file <- downloadGenomicSequencesFromNCBI(assembly_accession)
    }
    sorted_twobit_file <- file.path(tempdir(), "single_sequences.2bit")
    fastaTo2bit(fasta_file, sorted_twobit_file, assembly_accession)

    organism <- format_organism(organism)
    abbr_organism <- abbreviate_organism_name(organism)
    pkgname <- .make_pkgname_for_NCBI_datapkg(abbr_organism, assembly_name)
    pkgtitle <- .make_pkgtitle_for_NCBI_datapkg(organism, assembly_name)
    pkgdesc <- .make_pkgdesc_for_NCBI_datapkg(organism, assembly_name,
                                              assembly_accession)
    check_pkg_maintainer(pkg_maintainer)
    biocview <- organism2biocview(organism)
    seqnames <- build_Rexpr_as_string(chrominfo[ , "SequenceName"])
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
