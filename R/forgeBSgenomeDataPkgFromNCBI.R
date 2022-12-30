forgeBSgenomeDataPkgFromNCBI <- function(assembly_accession, organism, genome,
                                         pkg_maintainer, pkg_author=NA,
                                         pkg_version="0.01",
                                         license="Artistic-2.0",
                                         destdir=".")
{
    if (!isSingleString(organism))
        stop(wmsg("'organism' must be a single string"))
    if (!isSingleString(genome))
        stop(wmsg("'genome' must be a single string"))
    if(is.na(pkg_author))
        pkg_author <- pkg_maintainer

    ## Download file and convert from fasta to 2bit
    local_file <-downloadGenomicSequencesFromNCBI(assembly_accession)
    origfile <- local_file
    fastaTo2bit(origfile, "single_sequences.2bit", assembly_accession)

    ## Abbreviate package name
    f_name <- substring(organism,1,1)
    split_organism <- strsplit(organism, " +")[[1]]
    l_name <- tail(split_organism, 1)
    abbr_name <- paste0(f_name, l_name)
    pkgname <- paste0("BSgenome.", abbr_name, ".NCBI.", genome)

    pkgtitle <- paste0("Full genome sequences for ", organism, " (NCBI version ", genome, ")")
    pkgdesc <- paste0("Full genome sequences for ", organism, "as provided by NCBI (", genome, ") and stored in Biostrings objects.")

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
                      ORGANISMBIOCVIEW = genome,
                      SEQNAMES = "NA",
                      CIRCSEQS = "character(0)")

    originDir <- system.file("pkgtemplates", "NCBI_BSgenome_datapkg", package = "BSgenomeForge")
    file_dir <- createPackage(pkgname, destdir, originDir, symValues, unlink=TRUE, quiet=FALSE)

    current_dir <- getwd()
    new_dir <- file.path(file_dir, "inst/extdata")
    file.copy(from = file.path(current_dir, "single_sequences.2bit"),
              to = file.path(new_dir, "single_sequences.2bit"))
}
