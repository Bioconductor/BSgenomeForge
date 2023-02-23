### =========================================================================
### downloadGenomicSequencesFromUCSC()
### -------------------------------------------------------------------------
###
### A utility function to download the 2bit file that contains the genomic
### sequences of a given UCSC genome.
###

### NOT exported.
get_URL_to_genomic_sequences_from_UCSC <- function(genome,
            goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    file_name <- paste0(genome, ".2bit")
    paste(goldenPath.url, genome, "bigZips", file_name, sep="/")
}

downloadGenomicSequencesFromUCSC <- function(genome,
            goldenPath.url=getOption("UCSC.goldenPath.url"),
            destdir=".", method, quiet=FALSE)
{
    if (!isSingleString(genome) || genome == "")
        stop(wmsg("'genome' must be a single (non-empty) string"))
    if (!isSingleString(goldenPath.url) || goldenPath.url == "")
        stop(wmsg("'goldenPath.url' must be a single (non-empty) string"))
    if (!isSingleString(destdir))
        stop(wmsg("'destdir' must be a single string"))
    if (!dir.exists(destdir))
        stop(wmsg("'destdir' must be the path to an existing directory"))
    if (file.access(destdir, 2))
        stop(wmsg("you don't have write permission to 'destdir'"))
    if (!isTRUEorFALSE(quiet))
        stop(wmsg("'quiet' must be TRUE or FALSE"))

    file_url <- get_URL_to_genomic_sequences_from_UCSC(genome,
                                        goldenPath.url=goldenPath.url)
    destfile <- file.path(destdir, basename(file_url))
    download.file(file_url, destfile, method, quiet)
    invisible(destfile)
}

