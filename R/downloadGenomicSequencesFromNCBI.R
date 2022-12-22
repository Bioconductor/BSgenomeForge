### =========================================================================
### downloadGenomicSequencesFromNCBI()
### -------------------------------------------------------------------------
###
### A utility function to download the compressed FASTA file that contains
### the genomic sequences of a given NCBI assembly.
###

downloadGenomicSequencesFromNCBI <- function(assembly_accession,
                                             assembly_name=NA,
                                             destdir=".", method, quiet=FALSE)
{
    if (!isSingleString(destdir))
        stop(wmsg("'destdir' must be a single string"))
    if (!dir.exists(destdir))
        stop(wmsg("'destdir' must be the path to an existing directory"))
    if (file.access(destdir, 2))
        stop(wmsg("you don't have write permission to 'destdir'"))
    if (!isTRUEorFALSE(quiet))
        stop(wmsg("'quiet' must be TRUE or FALSE"))

    ftp_dir <- find_NCBI_assembly_ftp_dir(assembly_accession, assembly_name)
    ftp_dir_url <- ftp_dir[[1L]]
    file_prefix <- ftp_dir[[2L]]
    file_name <- paste0(file_prefix, "_genomic.fna.gz")
    file_url <- paste0(ftp_dir_url, "/", file_name)
    destfile <- file.path(destdir, file_name)
    download.file(file_url, destfile, method, quiet)
    invisible(destfile)
}

