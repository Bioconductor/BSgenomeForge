
.run_forgeBSgenomeDataPkgFromUCSC_in_extdata <- function(genome, organism)
{
    extdata_path <- system.file("extdata", package="BSgenomeForge",
                                mustWork=TRUE)
    old_wd <- setwd(extdata_path)
    on.exit(setwd(old_wd))
    forgeBSgenomeDataPkgFromUCSC(genome=genome,
                                 organism=organism,
                                 pkg_maintainer="Jane Doe <janedoe@gmail.com>",
                                 destdir=tempdir())
}

test_that("forgeBSgenomeDataPkgFromUCSC() works on \"eboVir3\"", {
    genome <- "eboVir3"
    organism <- "Ebola Virus"
    pkg_dir <- .run_forgeBSgenomeDataPkgFromUCSC_in_extdata(genome, organism)

    expect_true(dir.exists(pkg_dir))
    expect_identical(basename(pkg_dir), "BSgenome.Evirus.UCSC.eboVir3")
    path <- file.path(pkg_dir, "DESCRIPTION")
    expect_true(file.exists(path))
    path <- file.path(pkg_dir, "NAMESPACE")
    expect_true(file.exists(path))
    path <- file.path(pkg_dir, "R", "zzz.R")
    expect_true(file.exists(path))
    path <- file.path(pkg_dir, "inst", "extdata", "single_sequences.2bit")
    expect_true(file.exists(path))
})

