### =========================================================================
### forgeBSgenomeDataPkgFromNCBI()
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

.check_circ_seqs <- function(circ_seqs)
{
    # check if circ_seqs is null first - if so, return immediately
    if (is.null(circ_seqs))
        return(circ_seqs)
        # check is circ_seqs is a character vector & stop w/ error message if not
    if (!is.character(circ_seqs))
        stop(wmsg("'circ_seqs' must be NULL or a valid character vector"))
            # check for nas, empty strings and duplicates
    if (anyNA(circ_seqs))
        stop(wmsg("'circ_seqs' must contain valid, non empty character values"))
    if ("" %in% circ_seqs)
        stop(wmsg("'circ_seqs' must be a non-empty string"))
    if (anyDuplicated(circ_seqs))
        stop(wmsg("'circ_seqs' contains duplicate values"))
    return(circ_seqs)
}

.get_circseqs <- function(assembly_accession, circ_seqs=NULL)
{
    circ_seqs <- .check_circ_seqs(circ_seqs)
    NCBI_assemblies <- registered_NCBI_assemblies()
    seq_info <- getChromInfoFromNCBI(assembly_accession)
    ## if NCBI assembly is registered
    if (assembly_accession %in% NCBI_assemblies[ , "assembly_accession"]) {
        true_circ_seq <- seq_info[seq_info$circular == "TRUE", ]
        inferred_circ_seqs<- true_circ_seq$SequenceName

        if (is.null(circ_seqs))
            return(inferred_circ_seqs)
        if (!setequal(circ_seqs, inferred_circ_seqs))
            stop(wmsg("'circ_seqs' values do not match those of the assembly"))
            return(circ_seqs)

    } else {
        ## if NCBI assembly is not registered.
        if (is.null(circ_seqs))
            stop(wmsg("This assembly is not registered in
            GenomeInfoDb so I don't know what sequences in this assembly are
            circular. Please provide them in a character vector passed to
            the 'circ_seqs' argument (set 'circ_seqs' to 'character(0)'
                      if the assembly has no circular sequences)."))
        ## Check if circ_seqs match assembly sequence names
        if(anyNA(match(circ_seqs, seq_info$SequenceName)))
            stop(wmsg("Please enter valid circular sequence names"))
        ## Check if circ_seqs are names of assembled molecules
        ## Approach creates Warning: longer object length is not a multiple of
        ## shorter object length
        subset_seq_info <- seq_info[seq_info$SequenceName == circ_seqs, ]
        if (isEmpty(subset_seq_info)){
            return(circ_seqs)
        } else {
            if (! "assembled-molecule" %in% subset_seq_info[ , "SequenceRole"])
                stop(wmsg("'circ_seqs' must be an assembled molecule"))
        return(circ_seqs) }
    }
}

.build_Rexpr_as_string <- function(sequencename)
{
    if (length(sequencename) == 0)
        return("character(0)")
    paste0('c', '(', paste0('"', sequencename, '"', collapse=","), ')')
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
                                         circ_seqs=NULL,
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
#    if (anyNA(circ_seqs))
#        stop(wmsg("'circ_seqs' must contain valid, non empty character values"))
#    if (anyDuplicated(circ_seqs))
#        stop(wmsg("'circ_seqs' contains duplicate values"))
#    if (!(is.null(circ_seqs) || is.character(circ_seqs)))
#        stop(wmsg("'circ_seqs' must be NULL or a valid character vector"))
    if (!isSingleString(destdir) || destdir == "")
        stop(wmsg("'destdir' must be a single (non-empty) string"))

    ## Download file and convert from FASTA to 2bit.
    fasta_file <- downloadGenomicSequencesFromNCBI(assembly_accession)
    twobit_file <- file.path(tempdir(), "single_sequences.2bit")
    fastaTo2bit(fasta_file, twobit_file, assembly_accession)

    # retrieve the sequence names for supplied NCBI assembly
    seqinfo <- getChromInfoFromNCBI(assembly_accession)
    seqnames <- seqinfo$SequenceName

    organism <- .format_organism(organism)
    abbr_organism <- .abbreviate_organism_name(organism)
    pkgname <- .create_pkgname(abbr_organism, genome)
    pkgtitle <- .create_pkgtitle(organism, genome)
    pkgdesc <- .create_pkgdesc(organism, genome, assembly_accession)
    pkg_maintainer <- .check_pkg_maintainer(pkg_maintainer)
    organism_biocview <- .create_organism_biocview(organism)
    seqnames <- .build_Rexpr_as_string(seqnames)
    circ_seqs <- .check_circ_seqs(circ_seqs)
    circ_seqs <- .get_circseqs(assembly_accession, circ_seqs)
    circ_seqs <- .build_Rexpr_as_string(circ_seqs)

    symValues <- list(BSGENOMEOBJNAME=abbr_organism,
                      PKGTITLE=pkgtitle,
                      PKGDESCRIPTION=pkgdesc,
                      PKGVERSION=pkg_version,
                      PKGAUTHOR=pkg_author,
                      PKGMAINTAINER=pkg_maintainer,
                      PKGLICENSE=pkg_license,
                      ORGANISM=organism,
                      GENOME=genome,
                      ORGANISMBIOCVIEW=organism_biocview,
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
