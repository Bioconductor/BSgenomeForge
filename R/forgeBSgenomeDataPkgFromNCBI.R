.check_organism_case <- function(organism)
{
    f_letter <- substring(organism,1,1)
    r_letter <- substring(organism,2)
    if (grepl("[[:upper:]]", f_letter) == TRUE & grepl("[[:upper:]]", r_letter) == FALSE){
        return(organism)
    }
    else{
        stop(wmsg("'organism' must be in sentence case"))
    }
}

.get_abbr_organism <- function(organism)
{
    ## Abbreviate organism name
    organism <- .check_organism_case(organism)
    f_name <- substring(organism,1,1)
    split_organism <- strsplit(organism, " +")[[1]]
    l_name <- tail(split_organism, 1)
    abbr_organism <- paste0(f_name, l_name)
    abbr_organism
}

.make_pkgname_part4 <- function(genome) gsub("[^0-9a-zA-Z.]", "", genome)

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

.create_pkgname <- function(abbr_organism, genome)
{
    genome_name <- .make_pkgname_part4(genome)
    pkgname <- paste0("BSgenome.", abbr_organism, ".NCBI.", genome_name)
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
    split_organism <- strsplit(organism, " +")[[1]]
    h_name <- head(split_organism, 1)
    l_name <- tail(split_organism, 1)
    org_score <- paste0(h_name, "_", l_name)
    org_score
}

.move_seq_file <- function(file_dir)
{
    current_dir <- getwd()
    new_dir <- file.path(file_dir, "inst", "extdata")
    file.rename(from = file.path(current_dir, "single_sequences.2bit"),
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
    if (!isSingleString(pkg_maintainer))
        stop(wmsg("'package maintainer' must be a single string"))
    if (is.na(organism) || organism == '')
        stop(wmsg("'organism' must be a valid string value"))
    if (is.na(genome) || genome == '')
        stop(wmsg("'genome' must be a valid string value"))
    if (!(grepl(pattern = "@", pkg_maintainer) & grepl(pattern = ".com", pkg_maintainer)))
        stop(wmsg("Please enter a valid email address"))
    if (is.na(pkg_author))
        pkg_author <- pkg_maintainer

    ## Download file and convert from fasta to 2bit
    origfile <- downloadGenomicSequencesFromNCBI(assembly_accession)
    fastaTo2bit(origfile, "single_sequences.2bit", assembly_accession)

    abbr_organism <- .get_abbr_organism(organism)
    pkgname <- .create_pkgname(abbr_organism, genome)
    pkgtitle <- .create_pkgtitle(organism, genome)
    pkgdesc <- .create_pkgdesc(organism, genome)
    org_score <- .create_org_score(organism)
    seqnames <- .get_seqnames(assembly_accession)

    symValues <- list(PKGNAME = pkgname,
                      BSGENOMEOBJNAME = abbr_organism,
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
    file_dir <- unlist(createPackage(pkgname, destdir, originDir, symValues,
                              unlink=TRUE, quiet=FALSE))

    .move_seq_file(file_dir)
    return(file_dir)
}
