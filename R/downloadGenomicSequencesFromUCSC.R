### =========================================================================
### downloadGenomicSequencesFromUCSC()
### -------------------------------------------------------------------------
###
### A utility function to download the 2bit file that contains the genomic
### sequences of a given UCSC genome.
###

downloadGenomicSequencesFromUCSC <- function(genome,
            goldenPath.url=getOption("UCSC.goldenPath.url"),
            destdir=".", method, quiet=FALSE)
{
    if (!isSingleString(genome))
        stop(wmsg("'genome' must be a single string"))
    if (!isSingleString(goldenPath.url))
        stop(wmsg("'goldenPath.url' must be a single string"))
    if (!isSingleString(destdir))
        stop(wmsg("'destdir' must be a single string"))
    if (!dir.exists(destdir))
        stop(wmsg("'destdir' must be the path to an existing directory"))
    if (file.access(destdir, 2))
        stop(wmsg("you don't have write permission to 'destdir'"))
    if (!isTRUEorFALSE(quiet))
        stop(wmsg("'quiet' must be TRUE or FALSE"))

    file_name <- paste0(genome, ".2bit")
    file_url <- paste0(goldenPath.url, "/", genome, "/bigZips/", file_name)
    destfile <- file.path(destdir, file_name)
    download.file(file_url, destfile, method, quiet)
    invisible(destfile)
}

