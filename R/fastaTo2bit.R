### The function helps to convert a file from FASTA to 2bit
fastaTo2bit <- function(origfile, destfile, assembly_accession = NA)
  
{
    if (!isSingleString(origfile))
        stop(wmsg("'origfile' must be a single string"))
    if (dir.exists(origfile))
        stop(wmsg("'origfile' must be path to a file, not a directory"))
    if (!file.exists(origfile))
        stop(wmsg("'origfile' must be the path to an existing file"))
    if (!isSingleString(destfile))
        stop(wmsg("'destfile' must be a single string"))
    dna <- readDNAStringSet(origfile)
    
    current_RefSeqAccn <- unlist(heads(strsplit(names(dna), " ", fixed=TRUE), n=1L)) 
    chrominfo <- getChromInfoFromNCBI(assembly_accession)
    expected_RefSeqAccn <- chrominfo[ , "RefSeqAccn"]
    stopifnot(setequal(expected_RefSeqAccn, current_RefSeqAccn))
    
    ### Reorder the sequences.
    dna <- dna[match(expected_RefSeqAccn, current_RefSeqAccn)]
    
    ### Rename the sequences.
    names(dna) <- chrominfo[ , "SequenceName"]
    
    ### Check sequence lengths.
    expected_seqlengths <- chrominfo[ , "SequenceLength"]
    stopifnot(all(width(dna) == expected_seqlengths | is.na(expected_seqlengths)))
    
    ### Export file as 2bit
    export.2bit (dna, destfile)
}


