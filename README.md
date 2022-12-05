[<img src="https://www.bioconductor.org/images/logo/jpg/bioconductor_logo_rgb.jpg" width="200" align="right"/>](https://bioconductor.org/)

WORK IN PROGRESS!

**BSgenomeForge** is an R/Bioconductor package that provides a set of tools to forge _BSgenome data packages_. These tools supersede the old tools from the **BSgenome** software package.

### About BSgenome data packages

_BSgenome data packages_ are one of the many types of annotation packages available in Bioconductor. They contain the genomic sequences, which comprise chromosome sequences and other DNA sequences, of a particular genome assembly for a given organism. For example **BSgenome.Hsapiens.UCSC.hg19** is a BSgenome data package that contains the genomic sequences of the `hg19` genome from UCSC. Users can easily and efficiently access the sequences, or portions of the sequences, stored in these packages, via a common API implemented in the **BSgenome** _software_ package.

Bioconductor currently provides more than 100 BSgenome data packages, for more than 30 organisms. Most of them contain the genomic sequences of UCSC genomes (i.e. genomes supported by the UCSC Genome Browser) or NCBI assemblies. The packages are used in various Bioconductor workflows, as well as in man page examples and vignettes of other Bioconductor packages, typically in conjunction with tools available in the **BSgenome** and **Biostrings** software packages. New BSgenome data packages get added on a regular basis, based on user demand.

### Useful links

- **BSgenome.Hsapiens.UCSC.hg19** landing page: https://bioconductor.org/packages/BSgenome.Hsapiens.UCSC.hg19

- List of BSgenome data packages available in Bioconductor: https://bioconductor.org/packages/release/BiocViews.html#___BSgenome

- **BSgenome** software package landing page (with link to the old "How to forge a BSgenome data package" vignette): https://bioconductor.org/packages/BSgenome

- **BSgenome** software package on GitHub: https://github.com/Bioconductor/BSgenome

- **Biostrings** software package landing page: https://bioconductor.org/packages/Biostrings

- **GenomeInfoDb** software package on GitHub: https://github.com/Bioconductor/GenomeInfoDb

