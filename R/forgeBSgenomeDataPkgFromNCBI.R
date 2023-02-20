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

.create_pkgname_for_NCBI_datapkg <- function(abbr_organism, genome)
{
    assembly_name <- gsub("[^0-9a-zA-Z.]", "", genome)
    paste0("BSgenome.", abbr_organism, ".NCBI.", assembly_name)
}

.create_pkgtitle_for_NCBI_datapkg <- function(organism, genome)
{
    paste0("Full genomic sequences for ", organism,
           " (NCBI assembly ", genome, ")")
}

.create_pkgdesc_for_NCBI_datapkg <- function(organism, genome,
                                             assembly_accession)
{
    paste0("Full genomic sequences for ", organism, " as ",
           "provided by NCBI (assembly ", genome, ", assembly ",
           "accession ", assembly_accession, "). ",
           "The sequences are stored in DNAString objects.")
}

.move_seq_file <- function(twobit_file, pkg_dir)
{
    to <- file.path(pkg_dir, "inst", "extdata", basename(twobit_file))
    file.rename(twobit_file, to)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .get_circ_seqs_from_NCBI()
###
### TODO: Refactor .get_circ_seqs_from_NCBI() and .get_circ_seqs_from_UCSC()
### so that they share as much code as possible.

.check_circ_seqs <- function(circ_seqs)
{
    if (is.null(circ_seqs))
        return(circ_seqs)
    if (!is.character(circ_seqs))
        stop(wmsg("'circ_seqs' must be NULL or a character vector"))
    if (anyNA(circ_seqs))
        stop(wmsg("'circ_seqs' cannot contain NA's"))
    if ("" %in% circ_seqs)
        stop(wmsg("'circ_seqs' cannot contain empty strings"))
    if (anyDuplicated(circ_seqs))
        stop(wmsg("'circ_seqs' cannot contain duplicate values"))
    circ_seqs
}

.assembly_has_no_assembled_molecules <- function(circ_seqs)
{
    if (length(circ_seqs) == 0)
        return(character(0))
    stop(wmsg("This assembly contains no assembled molecules ",
              "so cannot have circular sequences."))
}

### 'circ_seqs' contains the circular sequences specified by the user.
.get_circ_seqs_from_NCBI <- function(assembly_accession, chrominfo,
                                     circ_seqs=NULL)
{
    circ_seqs <- .check_circ_seqs(circ_seqs)
    NCBI_assemblies <- registered_NCBI_assemblies()
    if (assembly_accession %in% NCBI_assemblies[ , "assembly_accession"]) {
        ## NCBI assembly is registered.
        known_circ_seqs <- chrominfo[chrominfo[ , "circular"], "SequenceName"]
        if (is.null(circ_seqs))
            return(known_circ_seqs)
        if (setequal(circ_seqs, known_circ_seqs))
            return(circ_seqs)
        msg <- "This assembly is registered in the GenomeInfoDb package "
        if (length(known_circ_seqs) == 0) {
            msg <- c(msg, "and it has no known circular sequences.")
        } else {
            in1string <- paste0("\"", known_circ_seqs, "\"", collapse=", ")
            msg <- c(msg, "which means that its circular sequences are known ",
                          "so you are not required to specify them. However, ",
                          "if you do specify them, then they must match the ",
                          "known ones. The circular sequences for registered ",
                          "assembly ", assembly_accession, " are: ", in1string)
        }
        stop(wmsg(msg))
    } else {
        ## NCBI assembly is **not** registered.
        if (!("assembled-molecule" %in% chrominfo[ , "SequenceRole"]))
            return(.assembly_has_no_assembled_molecules(circ_seqs))
        if (is.null(circ_seqs))
            stop(wmsg("This assembly is not registered in the ",
                      "GenomeInfoDb package so I don't know what its ",
                      "circular sequences are (if any). Please provide ",
                      "their names in a character vector passed to ",
                      "the 'circ_seqs' argument (set 'circ_seqs' to ",
                      "character(0) if the assembly has no circular ",
                      "sequences)."))
        ## The sequence names in 'circ_seqs' must belong to the assembly.
        if (anyNA(match(circ_seqs, chrominfo[ , "SequenceName"])))
            stop(wmsg("'circ_seqs' contains sequence names that ",
                      "do not belong to the specified assembly (",
                      assembly_accession, ")"))
        ## 'is_assembled' will be a logical vector with 1 element per row
        ## in the 'chrominfo' data frame.
        is_assembled <- chrominfo[ , "SequenceRole"] %in% "assembled-molecule"
        assembled_molecules <- chrominfo[is_assembled, "SequenceName"]
        ## The sequence names in 'circ_seqs' must be names of assembled
        ## molecules.
        if (!all(circ_seqs %in% assembled_molecules))
            stop(wmsg("all the sequence names in 'circ_seqs' must be ",
                      "names of assembled molecules"))
        return(circ_seqs)
    }
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
    genome <- .fetch_assembly_name_from_NCBI(assembly_accession)

    ## Download genomic sequences and convert from FASTA to 2bit.
    file_url <- get_URL_to_genomic_sequences_from_NCBI(assembly_accession)
    fasta_file <- basename(file_url)
    if (file.exists(fasta_file)) {
        message(wmsg("The file ", fasta_file, " is already in the current ",
                     "directory so will be used."))
    } else {
        fasta_file <- downloadGenomicSequencesFromNCBI(assembly_accession)
    }
    twobit_file <- file.path(tempdir(), "single_sequences.2bit")
    fastaTo2bit(fasta_file, twobit_file, assembly_accession)

    organism <- format_organism(organism)
    abbr_organism <- abbreviate_organism_name(organism)
    pkgname <- .create_pkgname_for_NCBI_datapkg(abbr_organism, genome)
    pkgtitle <- .create_pkgtitle_for_NCBI_datapkg(organism, genome)
    pkgdesc <- .create_pkgdesc_for_NCBI_datapkg(organism, genome,
                                                assembly_accession)
    check_pkg_maintainer(pkg_maintainer)
    biocview <- organism2biocview(organism)
    seqnames <- build_Rexpr_as_string(chrominfo[ , "SequenceName"])
    circ_seqs <- build_Rexpr_as_string(circ_seqs)

    symValues <- list(BSGENOMEOBJNAME=abbr_organism,
                      PKGTITLE=pkgtitle,
                      PKGDESCRIPTION=pkgdesc,
                      PKGVERSION=pkg_version,
                      PKGAUTHOR=pkg_author,
                      PKGMAINTAINER=pkg_maintainer,
                      PKGLICENSE=pkg_license,
                      ORGANISM=organism,
                      GENOME=genome,
                      ORGANISMBIOCVIEW=biocview,
                      SEQNAMES=seqnames,
                      CIRCSEQS=circ_seqs)

    origdir <- system.file("pkgtemplates", "NCBI_BSgenome_datapkg",
                            package="BSgenomeForge")
    pkg_dir <- unlist(createPackage(pkgname, destdir, origdir, symValues,
                                    unlink=TRUE, quiet=FALSE),
                      use.names=FALSE)

    .move_seq_file(twobit_file, pkg_dir)
    invisible(pkg_dir)
}
