### =========================================================================
### forgeBSgenomeDataPkgFromUCSC()
### -------------------------------------------------------------------------
###
### Create a BSgenome data package from a UCSC genome.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .get_circ_seqs_from_UCSC()
###
### TODO: Refactor .get_circ_seqs_from_UCSC() and .get_circ_seqs_from_NCBI()
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

.genome_has_no_assembled_molecules <- function(circ_seqs)
{
    if (length(circ_seqs) == 0)
        return(character(0))
    stop(wmsg("This UCSC genome contains no assembled molecules ",
              "so cannot have circular sequences."))
}

### 'circ_seqs' contains the circular sequences specified by the user.
.get_circ_seqs_from_UCSC <- function(genome, chrominfo, circ_seqs=NULL)
{
    circ_seqs <- .check_circ_seqs(circ_seqs)
    UCSC_genomes <- registered_UCSC_genomes()
    if (genome %in% UCSC_genomes[ , "genome"]) {
        ## UCSC genome is registered.
        known_circ_seqs <- chrominfo[chrominfo[ , "circular"], "SequenceName"]
        if (is.null(circ_seqs))
            return(known_circ_seqs)
        if (setequal(circ_seqs, known_circ_seqs))
            return(circ_seqs)
        msg <- "This UCSC genome is registered in the GenomeInfoDb package "
        if (length(known_circ_seqs) == 0) {
            msg <- c(msg, "and it has no known circular sequences.")
        } else {
            in1string <- paste0("\"", known_circ_seqs, "\"", collapse=", ")
            msg <- c(msg, "which means that its circular sequences are known ",
                          "so you are not required to specify them. However, ",
                          "if you do specify them, then they must match the ",
                          "known ones. The circular sequences for registered ",
                          "genome ", genome, " are: ", in1string)
        }
        stop(wmsg(msg))
    } else {
        ## UCSC genome is **not** registered.
        if (!("assembled-molecule" %in% chrominfo[ , "SequenceRole"]))
            return(.genome_has_no_assembled_molecules(circ_seqs))
        if (is.null(circ_seqs))
            stop(wmsg("This UCSC genome is not registered in the ",
                      "GenomeInfoDb package so I don't know what its ",
                      "circular sequences are (if any). Please provide ",
                      "their names in a character vector passed to ",
                      "the 'circ_seqs' argument (set 'circ_seqs' to ",
                      "character(0) if the genome has no circular ",
                      "sequences)."))
        ## The sequence names in 'circ_seqs' must belong to the genome.
        if (anyNA(match(circ_seqs, chrominfo[ , "SequenceName"])))
            stop(wmsg("'circ_seqs' contains sequence names that ",
                      "do not belong to the specified genome (",
                      genome, ")"))
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
### forgeBSgenomeDataPkgFromUCSC()
###

forgeBSgenomeDataPkgFromUCSC <- function(genome, organism,
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

    ## Retrieve chromosome information for specified UCSC genome.
    chrominfo <- getChromInfoFromUCSC(genome)
    circ_seqs <- .get_circ_seqs_from_UCSC(genome, chrominfo, circ_seqs)

    ## Download genomic sequences.
    twobit_file <- downloadGenomicSequencesFromUCSC(genome)

    organism <- format_organism(organism)
    abbr_organism <- abbreviate_organism_name(organism)

    stop("NOT READY YET!")
}

