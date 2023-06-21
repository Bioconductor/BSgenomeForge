
.run_forgeBSgenomeDataPkgFromNCBI_in_extdata <-
    function(assembly_accession, organism, circ_seqs=NULL)
{
    extdata_path <- system.file("extdata", package="BSgenomeForge",
                                mustWork=TRUE)
    old_wd <- setwd(extdata_path)
    on.exit(setwd(old_wd))
    forgeBSgenomeDataPkgFromNCBI(assembly_accession=assembly_accession,
                                 pkg_maintainer="Jane Doe <janedoe@gmail.com>",
                                 organism=organism,
                                 circ_seqs=circ_seqs,
                                 destdir=tempdir())
}

test_that("forgeBSgenomeDataPkgFromNCBI() works on \"Torque teno virus 1\"", {
    assembly_accession="GCF_000857545.1"
    organism <- "Torque teno virus 1"
    circ_seqs <- "NC_002076.2"
    pkg_dir <- .run_forgeBSgenomeDataPkgFromNCBI_in_extdata(
                                 assembly_accession, organism, circ_seqs)

    expect_true(dir.exists(pkg_dir))
    expect_identical(basename(pkg_dir), "BSgenome.Tvirus1.NCBI.ViralProj15247")
    path <- file.path(pkg_dir, "DESCRIPTION")
    expect_true(file.exists(path))
    path <- file.path(pkg_dir, "NAMESPACE")
    expect_true(file.exists(path))
    path <- file.path(pkg_dir, "R", "zzz.R")
    expect_true(file.exists(path))
    path <- file.path(pkg_dir, "inst", "extdata", "single_sequences.2bit")
    expect_true(file.exists(path))
})

