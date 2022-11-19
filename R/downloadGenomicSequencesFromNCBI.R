### The function helps to download NCBI genomic sequence of a particular organism
downloadGenomicSequencesFromNCBI <- function(assembly_accession,
                                             assembly_name=NA,
                                             destdir=".", method, quiet=FALSE)
{
    if (!isSingleString(destdir))
        stop(wmsg("'destdir' must be a single string"))
    if (!dir.exists(destdir))
        stop(wmsg("'destdir' must be the path to an existing directory"))
    if (file.access(destdir, 2))
        stop(wmsg("'destdir'  must be a path with permission"))
    if (!isTRUEorFALSE(quiet))
        stop(wmsg("'quiet' must be TRUE or FALSE"))

    ftp_dir <- find_NCBI_assembly_ftp_dir(assembly_accession, assembly_name)
    ftp_dir_url <- ftp_dir[[1]]
    file_prefix <- ftp_dir[[2]]
    filename <- paste0(file_prefix, "_", "genomic.fna.gz")
    full_url <- paste0(ftp_dir_url, "/", filename)
    destfile <- file.path(destdir, filename)
    download.file(full_url, destfile, method, quiet)
    invisible(destfile)
}
