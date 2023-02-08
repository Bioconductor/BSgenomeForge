
test_that("get_circ_seqs works for an unregistered assembly without
          assembled-molecules", {
    ## NCBI assembly Acidianus infernus
    chrominfo <- getChromInfoFromNCBI("GCA_009729545.1")
    test1 <- BSgenomeForge:::.get_circ_seqs("GCA_009729545.1", chrominfo)
    expect_equal(test1, character(0))
    test2 <- BSgenomeForge:::.get_circ_seqs("GCA_009729545.1", chrominfo,
                                           character(0))
    expect_equal(test2, character(0))
    expect_error(BSgenomeForge:::.get_circ_seqs("GCA_009729545.1",
                                                chrominfo, c("MT")),
                 "This assembly|circular sequences")
})

test_that("get_circ_seqs works for an unregistered assembly with
          assembled-molecules", {
    ## NCBI assembly Vibrio cholerae
    chrominfo <- getChromInfoFromNCBI("GCA_008369605.1")
    test3 <- BSgenomeForge:::.get_circ_seqs("GCA_008369605.1", chrominfo,
                                            circ_seqs = c("1", "2"))
    expect_equal(test3, c("1", "2"))
    test4 <- BSgenomeForge:::.get_circ_seqs("GCA_008369605.1", chrominfo,
                                            circ_seqs = c("unnamed"))
    expect_equal(test4, "unnamed")
    expect_error(BSgenomeForge:::.get_circ_seqs("GCA_008369605.1", chrominfo),
                 "This assembly is not registered|has no circular sequences")
    expect_error(BSgenomeForge:::.get_circ_seqs("GCA_008369605.1", chrominfo,
                                                "MT"), "'circ_seqs' contains|the
                 specified assembly")
})

test_that("get_circ_seqs works for an registered assembly with
          assembled-molecules", {
    ## NCBI assembly Felis_catus_9.0
    chrominfo <- getChromInfoFromNCBI("GCF_000181335.3")
    test5 <- BSgenomeForge:::.get_circ_seqs("GCF_000181335.3", chrominfo)
    expect_equal(test5, "MT")
    test6 <- BSgenomeForge:::.get_circ_seqs("GCF_000181335.3",
                                           chrominfo, circ_seqs = "MT")
    expect_equal(test6, "MT")
    expect_error(BSgenomeForge:::.get_circ_seqs("GCF_000181335.3", chrominfo,
                                                circ_seqs = "chrA1"),
                 "This assembly is registered in the GenomeInfoDb|circular
                 sequences for registered assembly")
    expect_error(BSgenomeForge:::.get_circ_seqs("GCF_000181335.3", chrominfo,
                                                circ_seqs = c("chrA1", "MT")),
                 "This assembly is registered in the GenomeInfoDb|circular
                 sequences for registered assembly")
})

test_that("get_circ_seqs works for an registered assembly with
          assembled-molecules and more than one circular sequence", {
    ## NCBI assembly Vitis vinifera
    chrominfo <- getChromInfoFromNCBI("GCF_000003745.3")
    test7 <- BSgenomeForge:::.get_circ_seqs("GCF_000003745.3", chrominfo)
    expect_equal(test7, c("MT", "Pltd"))
    test8 <- BSgenomeForge:::.get_circ_seqs("GCF_000003745.3",
                                            chrominfo, circ_seqs = c("MT",
                                                                     "Pltd"))
    expect_equal(test8, c("MT", "Pltd"))
    expect_error(BSgenomeForge:::.get_circ_seqs("GCF_000003745.3", chrominfo,
                                                circ_seqs = "1"), "This assembly
                                                is registered in the
                 GenomeInfoDb|circular sequences for registered assembly")
})

test_that("get_circ_seqs works for an registered assembly with
          assembled-molecules", {
    ## NCBI assembly Betacoronavirus
    chrominfo <- getChromInfoFromNCBI("GCA_009948555.1")
    test9 <- BSgenomeForge:::.get_circ_seqs("GCA_009948555.1", chrominfo)
    expect_equal(test9, character(0))
    test10 <- BSgenomeForge:::.get_circ_seqs("GCA_009948555.1", chrominfo,
                                             circ_seqs = character(0))
    expect_equal(test10, character(0))
    expect_error(BSgenomeForge:::.get_circ_seqs("GCA_009948555.1", chrominfo,
                                                circ_seqs = "MT"),
                 "This assembly is registered|no known circular sequences.")
})
