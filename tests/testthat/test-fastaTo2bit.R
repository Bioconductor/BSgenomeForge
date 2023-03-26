
sort_and_rename_fasta_sequences <-
    BSgenomeForge:::.sort_and_rename_fasta_sequences

test_that("sort_and_rename_fasta_sequences() works with GCA_009729545.1", {
    ## In practice 'dna' will be a DNAStringSet object obtained with
    ## readDNAStringSet() but sort_and_rename_fasta_sequences() doesn't
    ## really care about that. It only cares about 'dna' being a list-like
    ## object with names, and about the lengths of the individual list objects.
    dna <- list(Rle("A", 2186600),
                Rle("A", 16994),
                Rle("A", 3763),
                Rle("A", 15057))
    names(dna) <- c("WFIY01000004.1 blah blah",
                    "WFIY01000003.1 blah blah",
                    "WFIY01000002.1 blah blah blah",
                    "WFIY01000001.1 blah blah blah blah")
    expected_seqnames <- c("contig_1", "contig_2", "contig_3", "contig_4")

    dna2 <- sort_and_rename_fasta_sequences(dna, "GCA_009729545.1")
    expect_identical(names(dna2), expected_seqnames)

    dna3 <- sort_and_rename_fasta_sequences(rev(dna), "GCA_009729545.1")
    expect_identical(names(dna3), expected_seqnames)

    names(dna) <- c("WFIY01000004.1",
                    "WFIY01000003.1",
                    "WFIY01000002.1",
                    "WFIY01000001.1")
    dna4 <- sort_and_rename_fasta_sequences(rev(dna), "GCA_009729545.1")
    expect_identical(names(dna4), expected_seqnames)

    regexp <-
        "number[\\s]+of[\\s]+sequences[\\s]+.*[\\s]+does[\\s]+not[\\s]+match"
    expect_error(sort_and_rename_fasta_sequences(dna[-2], "GCA_009729545.1"),
                 regexp, ignore.case=TRUE, perl=TRUE)

    names(dna) <- tolower(names(dna))
    regexp <- "failed[\\s]+to[\\s]+map"
    expect_error(sort_and_rename_fasta_sequences(dna, "GCA_009729545.1"),
                 regexp, ignore.case=TRUE, perl=TRUE)

    names(dna) <- toupper(names(dna))
    dna$WFIY01000004.1 <- Rle("A", 9)
    regexp <-
        "lengths[\\s]+of[\\s]+sequences[\\s]+.*[\\s]+do[\\s]+not[\\s]+match"
    expect_error(sort_and_rename_fasta_sequences(rev(dna), "GCA_009729545.1"),
                 regexp, ignore.case=TRUE, perl=TRUE)
})

