
get_circ_seqs <- BSgenomeForge:::.get_circ_seqs

test_that("get_circ_seqs works for an unregistered assembly without
          assembled-molecules", {
    ## NCBI assembly Acidianus infernus (no circular sequences)
    chrominfo <- getChromInfoFromNCBI("GCA_009729545.1")
    circ_seqs <- get_circ_seqs("GCA_009729545.1", chrominfo)
    expect_equal(circ_seqs, character(0))
    circ_seqs <- get_circ_seqs("GCA_009729545.1", chrominfo, character(0))
    expect_equal(circ_seqs, character(0))
    regexp <- "cannot[\\s]+have[\\s]+circular[\\s]+sequences"
    expect_error(get_circ_seqs("GCA_009729545.1", chrominfo, c("MT")), regexp,
                 ignore.case=TRUE, perl=TRUE)
})

test_that("get_circ_seqs works for an unregistered assembly with
          assembled-molecules", {
    ## NCBI assembly Vibrio cholerae (2 circular sequences: 1 ,2)
    chrominfo <- getChromInfoFromNCBI("GCA_008369605.1")
    circ_seqs <- get_circ_seqs("GCA_008369605.1", chrominfo,
                               circ_seqs = c("1", "2"))
    expect_equal(circ_seqs, c("1", "2"))
    circ_seqs <- get_circ_seqs("GCA_008369605.1", chrominfo,
                               circ_seqs = c("unnamed"))
    expect_equal(circ_seqs, "unnamed")
    regexp <- "registered[\\s]+in[\\s]+the[\\s]+GenomeInfoDb[\\s]+package"
    expect_error(get_circ_seqs("GCA_008369605.1", chrominfo), regexp,
                 ignore.case=TRUE, perl=TRUE)
    regexp <- "'circ_seqs'[\\s]+contains[\\s]+sequence[\\s]+names[\\s]"
    expect_error(get_circ_seqs("GCA_008369605.1", chrominfo, "MT"), regexp,
                 ignore.case=TRUE, perl=TRUE)
})

test_that("get_circ_seqs works for an unregistered assembly with
          assembled-molecules and no circular sequences", {
    chrominfo <- getChromInfoFromNCBI("GCA_003692155.1")
    circ_seqs <- get_circ_seqs("GCA_003692155.1", chrominfo,
                               circ_seqs=character(0))
    expect_equal(circ_seqs, character(0))
    regexp <- "registered[\\s]+in[\\s]+the[\\s]+GenomeInfoDb[\\s]+package"
    expect_error(get_circ_seqs("GCA_003692155.1", chrominfo), regexp,
                 ignore.case=TRUE, perl=TRUE)
    regexp <- "'circ_seqs'[\\s]+contains[\\s]+sequence[\\s]+names[\\s]"
    expect_error(get_circ_seqs("GCA_003692155.1", chrominfo, "MT"), regexp,
                 ignore.case=TRUE, perl=TRUE)
})

test_that("get_circ_seqs works for a registered assembly with one
          assembled-molecule", {
    ## NCBI assembly Felis_catus_9.0 (1 circular sequence: MT)
    chrominfo <- getChromInfoFromNCBI("GCF_000181335.3")
    circ_seqs <- get_circ_seqs("GCF_000181335.3", chrominfo)
    expect_equal(circ_seqs, "MT")
    circ_seqs <- get_circ_seqs("GCF_000181335.3", chrominfo, circ_seqs = "MT")
    expect_equal(circ_seqs, "MT")
    regexp <- "circular[\\s]+sequences[\\s]+are[\\s]+known"
    expect_error(get_circ_seqs("GCF_000181335.3", chrominfo,
                               circ_seqs = "chrA1"), regexp, ignore.case=TRUE,
                 perl=TRUE)
    regexp <- "circular[\\s]+sequences[\\s]+are[\\s]+known"
    expect_error(get_circ_seqs("GCF_000181335.3", chrominfo,
                               circ_seqs = c("chrA1", "MT")), regexp,
                 ignore.case=TRUE, perl=TRUE)
})

test_that("get_circ_seqs works for a registered assembly with
          assembled-molecules and more than one circular sequence", {
    ## NCBI assembly Vitis vinifera (2 circular sequences: MT, Pltd)
    chrominfo <- getChromInfoFromNCBI("GCF_000003745.3")
    circ_seqs <- get_circ_seqs("GCF_000003745.3", chrominfo)
    expect_equal(circ_seqs, c("MT", "Pltd"))
    circ_seqs <- get_circ_seqs("GCF_000003745.3", chrominfo,
                               circ_seqs = c("MT", "Pltd"))
    expect_equal(circ_seqs, c("MT", "Pltd"))
    regexp <- "circular[\\s]+sequences[\\s]+are[\\s]+known"
    expect_error(get_circ_seqs("GCF_000003745.3", chrominfo, circ_seqs = "1"),
                 regexp, ignore.case=TRUE, perl=TRUE)
})

test_that("get_circ_seqs works for a registered assembly with no circular
          sequences", {
    ## NCBI assembly Betacoronavirus (no circular sequences)
    chrominfo <- getChromInfoFromNCBI("GCA_009948555.1")
    circ_seqs <- get_circ_seqs("GCA_009948555.1", chrominfo)
    expect_equal(circ_seqs, character(0))
    circ_seqs <- get_circ_seqs("GCA_009948555.1", chrominfo,
                               circ_seqs = character(0))
    expect_equal(circ_seqs, character(0))
    regexp <- "no[\\s]+known[\\s]+circular[\\s]+sequences"
    expect_error(get_circ_seqs("GCA_009948555.1", chrominfo, circ_seqs = "MT"),
                 regexp, ignore.case=TRUE, perl=TRUE)
    regexp <- "no[\\s]+known[\\s]+circular[\\s]+sequences"
    expect_error(get_circ_seqs("GCA_009948555.1", chrominfo,
                               circ_seqs = "MN996531.1"), regexp,
                 ignore.case=TRUE, perl=TRUE)
})

