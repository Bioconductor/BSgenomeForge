### The function helps to convert a file from FASTA to 2bit
fastaTo2bit <- function(origfile, destfile)
  
{
    if (!isSingleString(origfile))
        stop(wmsg("'origfile' must be a single string"))
    if (!file.exists(origfile))
        stop(wmsg("'origfile' must be the path to an existing file"))
    if (!isSingleString(destfile))
        stop(wmsg("'destfile' must be a single string"))
    dna <- readDNAStringSet(origfile)
    export.2bit(dna, destfile)
}

