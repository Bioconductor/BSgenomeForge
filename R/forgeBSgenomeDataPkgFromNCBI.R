.abbreviate_pkgname <- function(organism)
{
    ## Abbreviate package name
    f_name <- substring(organism,1,1)
    split_organism <- strsplit(organism, " +")[[1]]
    l_name <- tail(split_organism, 1)
    abbr_name <- paste0(f_name, l_name)
    abbr_name
}

.make_pkgname_part4 <- function(genome)
{
    if (grepl("[^0-9a-zA-Z.]", genome)) {
        genome_name <- gsub("[^0-9a-zA-Z.]", "", genome)
    } else {
        genome_name <- genome
    }
    genome_name
}

.get_seqnames <- function(assembly_accession)
{
    seqinfo <- getChromInfoFromNCBI(assembly_accession)
    seq_sub <- seqinfo$SequenceName
    seqnames <- paste0('c', '(', paste0('"', seq_sub, '"', collapse = ","), ')')
    seqnames
}

.create_pkgtitle <- function(organism, genome)
{
    pkgtitle <- paste0("Full genome sequences for ",
                       organism, " (NCBI version ", genome, ")")
    pkgtitle
}

.create_pkgname <- function(organism, genome)
{
    abbr_name <- .abbreviate_pkgname(organism)
    genome_name <- .make_pkgname_part4(genome)
    pkgname <- pkgname <- paste0("BSgenome.", abbr_name, ".NCBI.", genome_name)
    pkgname
}

.create_pkgdesc <- function(organism, genome)
{
    pkgdesc <- paste0("Full genome sequences for ",
                      organism,
                      "as provided by NCBI (",
                      genome, ") and stored in Biostrings objects.")
    pkgdesc
}

.create_org_score <- function(organism)
{
    org_score <- gsub(" ", "_", organism)
    org_score
}

.move_seq_file <- function(file_dir)
{
    new_dir <- file.path(file_dir, "inst/extdata")
    file.copy(from = file.path(tempdir(), "single_sequences.2bit"),
              to = file.path(new_dir, "single_sequences.2bit"))
}

forgeBSgenomeDataPkgFromNCBI <- function(assembly_accession, organism, genome,
                                         pkg_maintainer, pkg_author=NA,
                                         pkg_version="1.0.0",
                                         license="Artistic-2.0",
                                         destdir=".")
{
    if (!isSingleString(organism))
        stop(wmsg("'organism' must be a single string"))
    if (!isSingleString(genome))
        stop(wmsg("'genome' must be a single string"))
    if (is.na(pkg_author))
        pkg_author <- pkg_maintainer

    ## Download file and convert from fasta to 2bit
    origfile <- downloadGenomicSequencesFromNCBI(assembly_accession)
    twobitfile <- file.path(tempdir(), "single_sequences.2bit")
    fastaTo2bit(origfile, twobitfile, assembly_accession)

    abbr_name <- .abbreviate_pkgname(organism)
    pkgname <- .create_pkgname(organism, genome)
    pkgtitle <- .create_pkgtitle(organism, genome)
    pkgdesc <- .create_pkgdesc(organism, genome)
    org_score <- .create_org_score(organism)
    seqnames <- .get_seqnames(assembly_accession)

    symValues <- list(PKGNAME = pkgname,
                      BSGENOMEOBJNAME = abbr_name,
                      PKGTITLE = pkgtitle,
                      AUTHOR = pkg_author,
                      PKGVERSION = pkg_version,
                      MAINTAINER = pkg_maintainer,
                      PKGDESCRIPTION = pkgdesc,
                      LICENSE = license,
                      ORGANISM = organism,
                      GENOME = genome,
                      ORGANISMBIOCVIEW = org_score,
                      SEQNAMES = seqnames,
                      CIRCSEQS = "character(0)")

    originDir <- system.file("pkgtemplates", "NCBI_BSgenome_datapkg",
                             package = "BSgenomeForge")
    file_dir <- createPackage(pkgname, destdir, originDir, symValues,
                              unlink=TRUE, quiet=FALSE)

    .move_seq_file(file_dir)
}
