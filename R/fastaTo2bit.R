### =========================================================================
### fastaTo2bit()
### -------------------------------------------------------------------------


.is_GenBank_accession <- function(assembly_accession)
{
    stopifnot(isSingleString(assembly_accession))
    if (assembly_accession == "")
        stop(wmsg("'assembly_accession' cannot be an empty string"))
    ## Make sure that 'assembly_accession' looks either like a GenBank
    ## or like a RefSeq accession. We only look at its first 5 characters.
    if (grepl("^GCA_[0-9]", assembly_accession))
        return(TRUE)
    if (grepl("^GCF_[0-9]", assembly_accession))
        return(FALSE)
    stop(wmsg("malformed assembly accession: ", assembly_accession))
}

.extract_chrominfo_accns <- function(chrominfo, use_GenBankAccn)
{
    if (use_GenBankAccn) {
        chrominfo_accns <- chrominfo[ , "GenBankAccn"]
    } else {
        chrominfo_accns <- chrominfo[ , "RefSeqAccn"]
    }
    stopifnot(!anyDuplicated(chrominfo_accns))  # sanity check
    chrominfo_accns
}

.sort_and_rename_fasta_sequences <- function(dna, assembly_accession)
{
    dna_names <- names(dna)
    stopifnot(!is.null(dna_names))  # sanity check
    dna_accns <- unlist(heads(strsplit(dna_names, " ", fixed=TRUE), n=1L))
    stopifnot(!anyDuplicated(dna_accns))  # sanity check

    ## This also checks that 'assembly_accession' is either a GenBank
    ## or a RefSeq accession. Make sure to do this before the call to
    ## getChromInfoFromNCBI() below.
    is_GCA <- .is_GenBank_accession(assembly_accession)

    chrominfo <- getChromInfoFromNCBI(assembly_accession)
    if (length(dna) != nrow(chrominfo))
        stop(wmsg("number of sequences in FASTA file ",
                  "does not match number of sequences ",
                  "in 'getChromInfoFromNCBI(\"", assembly_accession, "\")'"))
    chrominfo_accns <- .extract_chrominfo_accns(chrominfo, is_GCA)

    ## Reorder the sequences in 'dna' as in 'chrominfo'. Note that at this
    ## point 'chrominfo_accns' and 'dna_accns' are guaranteed to have the
    ## same length.
    m <- match(chrominfo_accns, dna_accns)
    if (anyNA(m)) {
        what <- if (is_GCA) "GenBank" else "RefSeq"
        stop(wmsg("failed to map the ", what, " accessions ",
                  "in 'getChromInfoFromNCBI(\"", assembly_accession, "\")' ",
                  "to the sequence names in FASTA file"))
    }
    dna <- dna[m]

    ## Rename the sequences.
    names(dna) <- chrominfo[ , "SequenceName"]

    ## Check sequence lengths.
    chrominfo_seqlengths <- chrominfo[ , "SequenceLength"]
    if (!all(lengths(dna) == chrominfo_seqlengths |
             is.na(chrominfo_seqlengths)))
        stop(wmsg("lengths of sequences in FASTA file ",
                  "do not match lengths of sequences ",
                  "in 'getChromInfoFromNCBI(\"", assembly_accession, "\")'"))
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

