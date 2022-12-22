### =========================================================================
### fastaTo2bit()
### -------------------------------------------------------------------------


.sort_and_rename <- function(dna, assembly_accession)
{
    current_Accn <- unlist(heads(strsplit(names(dna), " ", fixed=TRUE), n=1L))
    chrominfo <- getChromInfoFromNCBI(assembly_accession)

    ## Check if RefSeq or GenBank assembly accession.
    if (grepl("GCF", assembly_accession)) {
        expected_Accn <- chrominfo[ , "RefSeqAccn"]
    } else if (grepl("GCA", assembly_accession)) {
        expected_Accn <- chrominfo[ , "GenBankAccn"]
    }
    stopifnot(setequal(expected_Accn, current_Accn))

    ## Reorder the sequences.
    dna <- dna[match(expected_Accn, current_Accn)]

    ## Rename the sequences.
    names(dna) <- chrominfo[ , "SequenceName"]

    ## Check sequence lengths.
    expected_seqlengths <- chrominfo[ , "SequenceLength"]
    stopifnot(all(width(dna) == expected_seqlengths |
                  is.na(expected_seqlengths)))
    dna
}

fastaTo2bit <- function(origfile, destfile, assembly_accession=NA)
{
    if (!isSingleString(origfile))
        stop(wmsg("'origfile' must be a single string"))
    if (dir.exists(origfile))
        stop(wmsg("'origfile' must be path to a file, not a directory"))
    if (!file.exists(origfile))
        stop(wmsg("'origfile' must be the path to an existing file"))
    if (!isSingleString(destfile))
        stop(wmsg("'destfile' must be a single string"))
    if (!isSingleStringOrNA(assembly_accession))
        stop(wmsg("'assembly_accession' must be a single string or NA value"))

    dna <- readDNAStringSet(origfile)

    if (!is.na(assembly_accession))
        dna <- .sort_and_rename(dna, assembly_accession)

    ### Export file as 2bit
    export.2bit (dna, destfile)
}
