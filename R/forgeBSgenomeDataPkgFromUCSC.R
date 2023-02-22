### =========================================================================
### forgeBSgenomeDataPkgFromUCSC()
### -------------------------------------------------------------------------
###
### Create a BSgenome data package from a UCSC genome.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .get_circ_seqs_from_UCSC()
###

### 'circ_seqs' contains the circular sequences specified by the user.
.get_circ_seqs_from_UCSC <- function(genome, chrominfo, circ_seqs=NULL)
{
    check_circ_seqs(circ_seqs)
    UCSC_genomes <- registered_UCSC_genomes()
    seqnames <- chrominfo[ , "chrom"]
    if (genome %in% UCSC_genomes[ , "genome"]) {
        ## UCSC genome is registered.
        FUN <- get_circ_seqs_for_registered_assembly_or_genome
        is_xxx <- chrominfo[ , "circular"]
    } else {
        ## UCSC genome is **not** registered.
        FUN <- get_circ_seqs_for_unregistered_assembly_or_genome
        is_xxx <- NULL
    }
    FUN(genome, seqnames, is_xxx, circ_seqs, what="UCSC genome")
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

