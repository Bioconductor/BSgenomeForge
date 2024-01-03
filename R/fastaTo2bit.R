### =========================================================================
### fastaTo2bit()
### -------------------------------------------------------------------------


.sort_and_rename_fasta_sequences <- function(dna, assembly_accession)
{
    dna_names <- names(dna)
    stopifnot(!is.null(dna_names))  # sanity check
    dna_accns <- unlist(heads(strsplit(dna_names, " ", fixed=TRUE), n=1L))
    stopifnot(!anyDuplicated(dna_accns))  # sanity check

    ## Returns TRUE if 'assembly_accession' is a GenBank accession, or FALSE
    ## if it's a RefSeq accession, or an error if it's none.
    ## Make sure to do this before the call to getChromInfoFromNCBI() below.
    is_GCA <- is_GenBank_accession(assembly_accession)
    accession_type <- if (is_GCA) "GenBank" else "RefSeq"
    accession_col <- paste0(accession_type, "Accn")

    ## Retrieve chromosome information for specified NCBI assembly.
    chrominfo <- getChromInfoFromNCBI(assembly_accession)
    chrominfo <- drop_rows_with_NA_accns(chrominfo, accession_col)
    if (length(dna) != nrow(chrominfo))
        stop(wmsg("Incompatible assembly report and FASTA file ",
                  "for assembly ", assembly_accession, ": ",
                  "the number of sequences in the FASTA file is ",
                  "not equal to the number of sequences in ",
                  "'getChromInfoFromNCBI(\"", assembly_accession, "\")' ",
                  "that have a non-NA ", accession_type, " accession"))
    chrominfo_accns <- chrominfo[ , accession_col]
    stopifnot(!anyDuplicated(chrominfo_accns))  # sanity check

    ## Reorder the sequences in 'dna' as in 'chrominfo'. Note that at this
    ## point 'chrominfo_accns' and 'dna_accns' are guaranteed to have the
    ## same length.
    m <- match(chrominfo_accns, dna_accns)
    if (anyNA(m))
        stop(wmsg("Incompatible assembly report and FASTA file ",
                  "for assembly ", assembly_accession, ": ",
                  "the non-NA ", accession_type, " accessions in ",
                  "'getChromInfoFromNCBI(\"", assembly_accession, "\")' ",
                  "cannot be mapped to the sequence names in the FASTA file"))
    dna <- dna[m]

    ## Rename the sequences.
    names(dna) <- chrominfo[ , "SequenceName"]

    ## Check sequence lengths.
    chrominfo_seqlengths <- chrominfo[ , "SequenceLength"]
    if (!all(lengths(dna) == chrominfo_seqlengths |
             is.na(chrominfo_seqlengths)))
        stop(wmsg("Incompatible assembly report and FASTA file ",
                  "for assembly ", assembly_accession, ": ",
                  "lengths of sequences in the FASTA file ",
                  "do not match lengths reported in ",
                  "'getChromInfoFromNCBI(\"", assembly_accession, "\")'"))
    dna
}

fastaTo2bit <- function(origfile, destfile, assembly_accession=NA)
{
    if (!isSingleString(origfile))
        stop(wmsg("'origfile' must be a single string"))
    if (dir.exists(origfile))
        stop(wmsg("'origfile' must be a path to a file, not to a directory"))
    if (!file.exists(origfile))
        stop(wmsg("'origfile' must be the path to an existing file"))
    if (!isSingleString(destfile))
        stop(wmsg("'destfile' must be a single string"))
    if (!isSingleStringOrNA(assembly_accession))
        stop(wmsg("'assembly_accession' must be a single string or NA value"))

    dna <- readDNAStringSet(origfile)
    if (!is.na(assembly_accession))
        dna <- .sort_and_rename_fasta_sequences(dna, assembly_accession)
    ## Export file as 2bit.
    export.2bit(dna, destfile)
}

