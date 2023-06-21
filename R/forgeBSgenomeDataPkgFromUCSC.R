### =========================================================================
### forgeBSgenomeDataPkgFromUCSC()
### -------------------------------------------------------------------------
###
### Create a BSgenome data package from a UCSC genome.
###


.make_pkgname_for_UCSC_datapkg <- function(abbr_organism, genome)
{
    assembly_name <- gsub("[^0-9a-zA-Z.]", "", genome)
    paste0("BSgenome.", abbr_organism, ".UCSC.", assembly_name)
}

.make_pkgtitle_for_UCSC_datapkg <- function(organism, genome)
{
    paste0("Full genomic sequences for ", organism,
           " (UCSC genome ", genome, ")")
}

.make_pkgdesc_for_UCSC_datapkg <- function(organism, genome)
{
    paste0("Full genomic sequences for ", organism, " as ",
           "provided by UCSC (genome ", genome, "). ",
           "The sequences are stored in DNAString objects.")
}

.sort_twobit_file <- function(origfile, destfile, chrominfo, genome)
{
    dna <- import.2bit(origfile)
    m <- match(chrominfo[ , "chrom"], names(dna))
    if (length(dna) != nrow(chrominfo) || anyNA(m))
        stop(wmsg("sequence names in file ", origfile, " are not the same ",
                  "as in 'getChromInfoFromUCSC(\"", genome, "\")'"))
    dna <- dna[m]
    if (!identical(width(dna), chrominfo[ , "size"]))
        stop(wmsg("sequence lengths in file ", origfile, " are not the same ",
                  "as in 'getChromInfoFromUCSC(\"", genome, "\")'"))
    export.2bit(dna, destfile)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .get_circ_seqs_from_UCSC()
###

### 'circ_seqs' contains the circular sequences specified by the user.
.get_circ_seqs_from_UCSC <- function(genome, is_registered, chrominfo,
                                     circ_seqs=NULL)
{
    check_circ_seqs(circ_seqs)
    seqnames <- chrominfo[ , "chrom"]
    if (is_registered) {
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
    goldenPath.url=getOption("UCSC.goldenPath.url"),
    destdir=".")
{
    if (!isSingleString(genome) || genome == "")
        stop(wmsg("'genome' must be a single (non-empty) string"))
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
    if (!isSingleString(goldenPath.url) || goldenPath.url == "")
        stop(wmsg("'goldenPath.url' must be a single (non-empty) string"))
    if (!isSingleString(destdir) || destdir == "")
        stop(wmsg("'destdir' must be a single (non-empty) string"))

    ## Is genome registered?
    UCSC_genomes <- registered_UCSC_genomes()[ , "genome"]
    is_registered <- genome %in% UCSC_genomes

    ## Retrieve chromosome information for specified UCSC genome.
    chrominfo <- getChromInfoFromUCSC(genome, goldenPath.url=goldenPath.url)
    circ_seqs <- .get_circ_seqs_from_UCSC(genome, is_registered, chrominfo,
                                          circ_seqs)

    ## Download genomic sequences.
    file_url <- get_URL_to_genomic_sequences_from_UCSC(genome,
                                              goldenPath.url=goldenPath.url)
    twobit_file <- basename(file_url)
    if (file.exists(twobit_file)) {
        message(wmsg("File ", twobit_file, " is already in current ",
                     "directory so will be used."))
    } else {
        twobit_file <- downloadGenomicSequencesFromUCSC(genome,
                                              goldenPath.url=goldenPath.url)
    }
    sorted_twobit_file <- file.path(tempdir(), "single_sequences.2bit")
    .sort_twobit_file(twobit_file, sorted_twobit_file, chrominfo, genome)

    organism <- format_organism(organism)
    abbr_organism <- abbreviate_organism_name(organism)
    pkgname <- .make_pkgname_for_UCSC_datapkg(abbr_organism, genome)
    pkgtitle <- .make_pkgtitle_for_UCSC_datapkg(organism, genome)
    pkgdesc <- .make_pkgdesc_for_UCSC_datapkg(organism, genome)
    check_pkg_maintainer(pkg_maintainer)
    biocview <- organism2biocview(organism)
    seqnames <- build_Rexpr_as_string(chrominfo[ , "chrom"])
    circ_seqs <- build_Rexpr_as_string(circ_seqs)

    ## Create the package.
    origdir <- system.file("pkgtemplates", "UCSC_BSgenome_datapkg",
                           package="BSgenomeForge")
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
    pkg_dir <- createPackage(pkgname, destdir, origdir, symValues,
                             unlink=TRUE, quiet=FALSE)[[1]]
    move_file_to_datapkg(sorted_twobit_file, pkg_dir)

    invisible(pkg_dir)
}

